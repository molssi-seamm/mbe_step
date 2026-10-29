# -*- coding: utf-8 -*-

"""Non-graphical part of the MBE step in a SEAMM flowchart"""

import json
import logging
from pathlib import Path
import importlib.resources
import pprint  # noqa: F401
import time

import numpy as np

import mbe_step
import molsystem
import seamm
import seamm_mbe
from seamm_util import ureg, Q_  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

from . import counterpoise
from . import labels as labels_
from . import levels
from . import thermo
from .mbe_parameters import CRITERIA

# In addition to the normal logger, two logger-like printing facilities are
# defined: "job" and "printer". "job" send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# "printer" sends output to the file "step.out" in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("MBE")

# Add this module's properties to the standard properties
path = importlib.resources.files("mbe_step") / "data"
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)


def _aux_listed(P):
    """The fragments listed for explicit auxiliary basis sets, by level:
    'auxiliary basis fragments' (increments of two or more units) and
    'auxiliary basis monomers' (1-body terms) together. {str: set of str}"""
    listed = levels.parse_aux_fragments(P["auxiliary basis fragments"])
    for level, keys in _aux_monomers(P).items():
        listed.setdefault(level, set()).update(keys)
    return listed


def _aux_monomers(P):
    """The 1-body terms listed for explicit auxiliary basis sets, by level."""
    return levels.parse_aux_fragments(P.get("auxiliary basis monomers", ""))


def _distance_set(value):
    """Whether a distance parameter is set (> 0): a Quantity when run() reads it,
    a string such as '2.5 Å' or a $variable in the description."""
    if hasattr(value, "m_as"):
        return value.m_as("Å") > 0
    try:
        return float(str(value).split()[0]) > 0
    except (ValueError, IndexError):
        return True  # e.g. a variable, known only when the step runs


def _truthy(value):
    return value is True or (isinstance(value, str) and value.lower() == "yes")


def _triples_level(P):
    """The triples' own high level, or None when it is the high level's or there
    are no triples."""
    return _triples_override(P, "triples high level", "high level")


def _triples_low_level(P):
    """The triples' own molecular low level, or None when it is the molecular low
    level or there are no triples."""
    return _triples_override(P, "triples low level", "molecular low level")


def _shell_level(P):
    """The shell fragments' own high level, or None when it is the high level's
    or there are no ion shells."""
    return _shell_override(P, "shell high level", "high level")


def _shell_low_level(P):
    """The shell fragments' own molecular low level, or None when it is the
    molecular low level or there are no ion shells."""
    return _shell_override(P, "shell low level", "molecular low level")


def _shell_override(P, key, base):
    """The shell fragments' own level from parameter ``key``, or None when it is
    blank, "same as the ...", the same text as ``base``, or there are no ion
    shells."""
    if P["ion shells"] == "none":
        return None
    text = (P[key] or "").strip()
    if text == "" or text.startswith("same as the") or text == P[base]:
        return None
    return text


def _level_names(P):
    """The level names ("high", "high:3", ...) that the task options may use with
    these settings. The cell and cluster levels depend on the structures, known
    only when the step runs, so both are allowed here."""
    names = {"high", "molecular", "cell", "cluster"}
    if _triples_level(P):
        names.add("high:3")
    if _triples_low_level(P):
        names.add("molecular:3")
    if _shell_level(P):
        names.add("high:shell")
    if _shell_low_level(P):
        names.add("molecular:shell")
    if P["periodic low level"] != "none":
        names.add("periodic")
    return names


def _level_maps(P):
    """The seamm_mbe level maps ({order or "shell": level name}) for the high and
    the molecular low levels, or None for each without levels of their own."""
    high, low = {}, {}
    if _triples_level(P):
        high[3] = "high:3"
    if _shell_level(P):
        high[seamm_mbe.SHELL] = "high:shell"
    if _triples_low_level(P):
        low[3] = "molecular:3"
    if _shell_low_level(P):
        low[seamm_mbe.SHELL] = "molecular:shell"
    return high or None, low or None


def _triples_override(P, key, base):
    """The triples' own level from parameter ``key``, or None when it is blank,
    "same as the ...", the same text as ``base``, or there are no triples."""
    text = (P[key] or "").strip()
    if text == "" or text.startswith("same as the") or text == P[base]:
        return None
    if P["triple rule"] == "none":
        return None
    try:
        order = int(P["maximum order"])
    except (TypeError, ValueError):
        order = 3
    return text if order >= 3 else None


class Mbe(seamm.Node):
    """
    The non-graphical part of a MBE step in a flowchart.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    parameters : MbeParameters
        The control parameters for MBE.

    See Also
    --------
    TkMbe,
    Mbe, MbeParameters
    """

    def __init__(self, flowchart=None, title="MBE", extension=None, logger=logger):
        """A step for MBE in a SEAMM flowchart.

        You may wish to change the title above, which is the string displayed
        in the box representing the step in the flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.

        title: str
            The name displayed in the flowchart.
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating MBE {self}")

        super().__init__(
            flowchart=flowchart,
            title="MBE",
            extension=extension,
            module=__name__,
            logger=logger,
        )  # yapf: disable

        self._metadata = mbe_step.metadata
        self.parameters = mbe_step.MbeParameters()
        self.model = None
        self._grid = {"max_spacing": 0.0829, "padding": 7.5}
        self._counterpoise = False
        self._counterpoise_low = False

    @property
    def version(self):
        """The semantic version of this module."""
        return mbe_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return mbe_step.__git_revision__

    def description_text(self, P=None):
        """Create the text description of what this step will do.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control
            parameters.
        Returns
        -------
        str
            A description of the current step.
        """
        if not P:
            P = self.parameters.values_to_dict()

        which = seamm.standard_parameters.structure_selection_description(P)
        which = which[0].lower() + which[1:]
        if which.endswith(" will be used."):
            which = which[: -len(" will be used.")]

        order = int(P["maximum order"])
        what = {1: "monomers", 2: "monomers and pairs"}.get(order, "")
        if order >= 3:
            what = (
                "monomers, pairs and triples"
                if P["triple rule"] != "none"
                else "monomers and pairs"
            )
        high = P["high level"]
        if high == levels.CURRENT:
            high = "the current model chemistry"
        text = (
            f"Correct a low-level calculation of {which} with many-body "
            f"increments of [high - low] from its {what}, with {high} as the "
            f"high level and {P['molecular low level'] or '(not set)'} as the "
            "molecular low level"
        )
        if P["periodic low level"] != "none":
            text += f", and {P['periodic low level']} for the compact fragments"
        text += "."
        triples = _triples_level(P)
        if triples:
            text += (
                f" The triples' increments use {triples} as their high level, with "
                "their pairs and monomers at that level too."
            )
        triples_low = _triples_low_level(P)
        if triples_low:
            text += (
                f" The triples' increments use {triples_low} as their molecular low "
                "level, with their pairs and monomers at that level too."
            )
        shell = _shell_level(P)
        shell_low = _shell_low_level(P)
        if shell or shell_low:
            parts = []
            if shell:
                parts.append(f"{shell} as their high level")
            if shell_low:
                parts.append(f"{shell_low} as their molecular low level")
            text += (
                " The increments of the fragments holding an ion shell use "
                + " and ".join(parts)
                + ", with their sub-fragments at that level too."
            )
        task_options = levels.parse_task_options(P["task options"])
        if task_options:
            text += (
                " Task options: "
                + "; ".join(
                    f"{name} " + " ".join(f"{k}={v}" for k, v in options.items())
                    for name, options in task_options.items()
                )
                + "."
            )
        listed = _aux_listed(P)
        if listed:
            sets = levels.parse_aux_sets(P["auxiliary basis sets"])
            text += (
                " These fragments are computed, with their sub-fragments, with "
                "explicit auxiliary basis sets instead of AutoAux: "
                + "; ".join(
                    f"{', '.join(sorted(keys))} at the {level} level "
                    f"({sets.get(level, 'explicit')})"
                    for level, keys in listed.items()
                )
                + "."
            )
            if _aux_monomers(P):
                text += (
                    " The 1-body energies of those units are not on the AutoAux "
                    "scale of the others."
                )
        tiers = levels.parse_resources(P["fragment resources"])
        if tiers:
            text += (
                " Larger molecular fragments run with more resources: "
                + "; ".join(
                    f"from {a} atoms {r} cores of {m:.0f} MB"
                    + (f", {b} per job" if b else "")
                    for a, r, m, b in tiers
                )
                + "."
            )
        if order >= 2:
            if P["cutoffs"] == "single value":
                text += f" Pairs closer than {P['pair cutoff']} are selected"
            else:
                text += " Pairs are selected with cutoffs by molecule type"
            if order >= 3 and P["triple rule"] != "none":
                text += f", and {P['triple rule']} triples"
            text += f" ({P['distance criterion']})."
            if P["ion shells"] != "none":
                text += (
                    " Each Li+ and the molecules in its first shell are one unit of "
                    f"the expansion (shell cutoffs {P['ion shell cutoffs']}"
                )
                if int(P["shell members"]) > 0:
                    text += f", at most {P['shell members']} molecules"
                text += ")"
                if P["shell truncation"] == "pairs":
                    text += ", and a shell's pairs are selected but not its triples"
                text += "."
                if _distance_set(P["shell pair cap"]):
                    text += (
                        f" Pairs of two shells {P['shell pair cap']} or more apart "
                        "are not computed (capped)."
                    )
            if P["counterpoise"] == "pairwise":
                text += (
                    " Each pair is corrected for the basis-set superposition error "
                    "(pairwise counterpoise) at the high level."
                )
            elif P["counterpoise"] == "pairwise, both levels":
                text += (
                    " Each pair is corrected for the basis-set superposition error "
                    "(pairwise counterpoise) at the high level and the molecular low "
                    "level."
                )
        offsets = (P["energy offsets"] or "").strip()
        if thermo.is_automatic(offsets):
            text += (
                " The energies are put on the DfE0 scale (energy of formation at 0 "
                "K) with the atomic references at the high level from the "
                "thermochemistry database."
            )
        elif offsets.lower() not in ("", "none"):
            text += f" Energy offsets (eV per molecule): {offsets}."
        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    # ------------------------------------------------------------------ #
    # Running
    # ------------------------------------------------------------------ #

    def run(self):
        """Run the MBE step.

        Returns
        -------
        seamm.Node
            The next node object in the flowchart.
        """
        next_node = super().run(printer)
        context = seamm.flowchart_variables._data
        P = self.parameters.current_values_to_dict(context=context)

        printer.important(__(self.description_text(P), indent=self.indent))
        printer.important("")

        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        configurations = self.select_configurations(P)
        if len(configurations) == 0:
            raise ValueError("The MBE step found no configurations to label.")
        rules = self._rules(P)
        automatic = thermo.is_automatic(P["energy offsets"])
        offsets = None if automatic else labels_.parse_offsets(P["energy offsets"])
        frames = [self._frame(c, rules, P) for c in configurations]
        self._check_aux_fragments(P, frames)
        levels_ = self._levels(P, context, frames)
        if automatic:
            # Before any calculation: a level without atomic references is refused
            offsets = self._dfe0_offsets(levels_["high"], frames)
        self._grid = {
            "max_spacing": P["grid spacing"].m_as("Å"),
            "padding": P["box padding"].m_as("Å"),
        }
        self._counterpoise = P["counterpoise"] != "none"
        self._counterpoise_low = P["counterpoise"] == "pairwise, both levels"
        self.model = levels_["high"].level

        t0 = time.perf_counter()
        results = self._evaluate(levels_, frames)
        elapsed = time.perf_counter() - t0

        failed = []
        rows = []
        for frame in frames:
            try:
                labels = self._assemble(frame, levels_, results, offsets)
            except seamm_mbe.MissingFragmentsError as e:
                failed.append((frame, self._failure_text(frame, e, results)))
                continue
            self._store(frame, labels, levels_, P)
            rows.append((frame, labels))

        self._cite()
        self.analyze(rows=rows, levels=levels_, elapsed=elapsed)
        if failed:
            lines = "\n".join(
                f"    {frame['configuration'].name}: {text}" for frame, text in failed
            )
            raise RuntimeError(
                f"{len(failed)} of {len(frames)} configurations are incomplete and "
                f"were given no labels; the others are stored. A rerun of the job "
                f"retries the failed calculations.\n{lines}"
            )
        return next_node

    # ------------------------------------------------------------- the setup
    def _rules(self, P):
        """The seamm_mbe selection rules from the parameters."""
        order = int(P["maximum order"])
        cutoffs, rules = {}, {}
        table = P["cutoffs"] != "single value"
        if order >= 2:
            cutoffs[2] = (
                labels_.parse_table(P["pair cutoff table"], "pair cutoff")
                if table
                else P["pair cutoff"].m_as("Å")
            )
        if order >= 3:
            cutoffs[3] = (
                labels_.parse_table(P["triple cutoff table"], "triple cutoff")
                if table
                else P["triple cutoff"].m_as("Å")
            )
            rules[3] = P["triple rule"]
        shells = self._shell_rules(P) is not None
        cap = P["shell pair cap"].m_as("Å")
        return seamm_mbe.SelectionRules(
            max_order=order,
            criterion=CRITERIA[P["distance criterion"]],
            cutoffs=cutoffs,
            rules=rules,
            shell_max_order=2 if shells and P["shell truncation"] == "pairs" else None,
            shell_pair_cap=cap if shells and cap > 0 else None,
        )

    def _shell_rules(self, P):
        """The seamm_mbe ion-shell rules from the parameters, or None."""
        if P["ion shells"] == "none":
            return None
        members = int(P["shell members"])
        return seamm_mbe.IonShellRules(
            cutoffs=labels_.parse_shells(P["ion shell cutoffs"]),
            max_members=members if members > 0 else None,
        )

    def _frame(self, configuration, rules, P):
        """The molecules and fragments of one configuration."""
        if configuration.n_bonds == 0 and configuration.n_atoms > 1:
            configuration.perceive_bonds()
            printer.important(
                __(
                    f"Configuration {configuration.name} had no bonds, so they "
                    f"were perceived: {configuration.n_bonds} bonds.",
                    indent=4 * " ",
                )
            )
        system = seamm_mbe.System.from_configuration(configuration)
        # With ion shells the fragments are built from units (each shell one
        # unit); the real system keeps the molecules, for the energy offsets,
        # the charges and the molecular virial.
        shell_rules = self._shell_rules(P)
        units = (
            system if shell_rules is None else seamm_mbe.ion_shells(system, shell_rules)
        )
        fragments = seamm_mbe.enumerate_fragments(units, rules)
        if fragments.capped:
            printer.important(
                __(
                    f"Configuration {configuration.name}: {len(fragments.capped)} "
                    f"pairs of two shells are {rules.shell_pair_cap:.2f} Å or more "
                    "apart, so they are not computed (capped).",
                    indent=4 * " ",
                )
            )
        periodic = {}
        if P["periodic low level"] != "none":
            periodic[1] = _truthy(P["periodic monomers"])
            r1 = P["periodic pair cutoff"].m_as("Å")
            if r1 > 0:
                periodic[2] = r1
        seamm_mbe.assign_levels(fragments, periodic)
        self._check_frame(configuration, system, rules, P)
        high_levels, low_levels = _level_maps(P)
        prefix = f"c{configuration.id}-"
        variants = self._aux_variants(
            configuration, fragments, prefix, high_levels, low_levels, P
        )
        return {
            "configuration": configuration,
            "system": system,
            "fragments": fragments,
            "calculations": fragments.calculations(
                high_levels=high_levels, low_levels=low_levels, variants=variants
            ),
            "prefix": prefix,
            "variants": variants,
        }

    @staticmethod
    def _aux_variants(configuration, fragments, prefix, high_levels, low_levels, P):
        """The seamm_mbe level variants of a configuration for the fragments
        listed to use explicit auxiliary basis sets: {level: (level@aux,
        {fragment name})}."""
        listed = _aux_listed(P)
        monomers = _aux_monomers(P)
        variants = {}
        for level, keys in listed.items():
            names = {k[len(prefix) :] for k in keys if k.startswith(prefix)}
            for name in sorted(names):
                where = f"{prefix}{name} (configuration {configuration.name})"
                if name not in fragments:
                    raise ValueError(
                        f"The auxiliary basis fragment {where} is not a fragment."
                    )
                f = fragments[name]
                if not f.in_sum:
                    raise ValueError(
                        f"The auxiliary basis fragment {where} is not a selected "
                        "fragment, so it has no increment of its own."
                    )
                one_body = prefix + name in monomers.get(level, set())
                if f.order < 2 and not one_body:
                    raise ValueError(
                        f"The auxiliary basis fragment {where} is a monomer: the "
                        "1-body terms set the energy scale of every unit and stay "
                        "at their level, unless listed in 'Explicit auxiliary basis "
                        "for 1-body terms'."
                    )
                if f.order >= 2 and one_body:
                    raise ValueError(
                        f"The 1-body auxiliary basis entry {where} is a fragment of "
                        f"{f.order} units: list it in 'Explicit auxiliary basis for' "
                        "instead."
                    )
                uses = (
                    seamm_mbe.level_key(f, high_levels, "high", fragments.shells),
                    seamm_mbe.low_level_of(f, low_levels, fragments.shells),
                )
                if level not in uses:
                    raise ValueError(
                        f"The auxiliary basis fragment {where} is given for the "
                        f"{level} level, but its increment uses {uses[0]} and "
                        f"{uses[1]}."
                    )
            if names:
                variants[level] = (f"{level}@aux", names)
        return variants

    @staticmethod
    def _check_aux_fragments(P, frames):
        """Refuse auxiliary basis fragments that name no configuration's fragment,
        e.g. a mistyped configuration number."""
        listed = _aux_listed(P)
        prefixes = [frame["prefix"] for frame in frames]
        for level, keys in listed.items():
            for key in sorted(keys):
                if not any(key.startswith(prefix) for prefix in prefixes):
                    raise ValueError(
                        f"The auxiliary basis fragment {key} ({level}) is not in "
                        "any of the configurations: the keys are 'c<configuration "
                        "id>-<fragment>', e.g. c1-d01_09."
                    )

    def _check_frame(self, configuration, system, rules, P):
        """Refuse, before any calculation, what would fail or be wrong later:
        missing energy offsets, a configuration charge that disagrees with the
        molecules', and open-shell molecules inside larger fragments."""
        name = configuration.name
        if thermo.is_automatic(P["energy offsets"]):
            offsets = None  # from the database for every type, in run()
        else:
            offsets = labels_.parse_offsets(P["energy offsets"])
        if offsets is not None:
            missing = sorted(set(system.type_counts()) - set(offsets))
            if missing:
                raise ValueError(
                    f"Configuration {name} has molecules of type(s) {missing} with "
                    "no energy offset. Give an offset for every type, or 'none'."
                )
        total = sum(system.types[m.type].charge for m in system.molecules)
        given = int(getattr(configuration, "charge", 0) or 0)
        if given != 0 and given != total:
            raise ValueError(
                f"Configuration {name} has charge {given:+d}, but its molecules' "
                f"charges add up to {total:+d}."
            )
        if P["ion shells"] != "none" and P["counterpoise"] != "none":
            raise ValueError(
                "Counterpoise corrections are not available with ion shells yet: "
                "the ghost calculations are built per molecule, not per shell."
            )
        open_shell = sorted(
            t.name for t in system.types.values() if t.multiplicity != 1
        )
        if open_shell and rules.max_order >= 2:
            raise ValueError(
                f"Configuration {name} contains open-shell molecules "
                f"({', '.join(open_shell)}); fragments of several molecules with "
                "open shells are not supported yet. Use maximum order 1."
            )

    def _dfe0_offsets(self, level, frames):
        """The energy offsets onto the DfE0 scale, eV per molecule by type, from
        the thermochemistry database's atomic references at the high level (the
        level of the monomers' energies, which set the labels' scale)."""
        try:
            provider = self.flowchart.plugin_manager.get(level.mc["step"])
        except Exception:
            provider = None
        offsets = thermo.dfe0_offsets(
            level.mc, [frame["system"] for frame in frames], provider
        )
        printer.important(
            __(
                f"Energies on the DfE0 scale, with the atomic references at "
                f"{level.level}: offsets (eV per molecule) "
                + ", ".join(f"{name} {value:.5f}" for name, value in offsets.items())
                + ".",
                indent=4 * " ",
            )
        )
        printer.important("")
        return offsets

    def check_flowchart(self, preceding):
        """Check, before the job runs, that the task options name levels that
        this step uses with its settings (MBE.run refuses the others)."""
        P = self.parameters.values_to_dict()
        text = P["task options"]
        if not isinstance(text, str) or "$" in text:
            return []
        try:
            options = levels.parse_task_options(text)
        except ValueError as e:
            return [f"'{self.title}': the task options: {e}"]
        names = _level_names(P)
        unused = [key for key in options if key not in names]
        if not unused:
            return []
        these = "level" if len(unused) == 1 else "levels"
        return [
            f"'{self.title}': the task options are for the {these} "
            f"{', '.join(unused)}, which this step does not use with its settings "
            f"(the levels are {', '.join(sorted(names))}). Remove those entries, or "
            "give that level its own model chemistry."
        ]

    def _levels(self, P, context, frames):
        """Resolve each level's model chemistry, checking what is needed."""
        current = (
            self.get_variable("_model_chemistry")
            if self.variable_exists("_model_chemistry")
            else None
        )
        hours = P["bundle walltime"].m_as("s")
        archive = _truthy(P["archive"])
        tiers = levels.parse_resources(P["fragment resources"])
        task_options = levels.parse_task_options(P["task options"])

        def make(name, text, periodic, prefix):
            mc = levels.resolve(text, context, current=current, periodic=periodic)
            return levels.Level(
                name,
                mc,
                ranks=P[f"{prefix} ranks"],
                memory=P[f"{prefix} memory"],
                bundle=P[f"{prefix} bundle"] if prefix != "cell" else 1,
                walltime=hours,
                archive=archive,
                tiers=tiers if prefix == "molecular" else None,
            )

        def add_options(key):
            """The task options for level ``key``, onto its Level; a Level
            shared by two names takes both, and a conflict is refused."""
            options = task_options.get(key)
            if not options or key not in result:
                return
            level = result[key]
            for option, value in options.items():
                if level.options.get(option, value) != value:
                    raise ValueError(
                        f"The task options give {option} two values for the "
                        f"levels computed together with {key}: "
                        f"{level.options[option]} and {value}."
                    )
            level.options.update(options)

        def own(key, name, text, base_keys, setting, same):
            """A level of its own, or the same Level as one of ``base_keys`` that
            has its model chemistry (then computed once, for both ladders)."""
            candidate = make(name, text, False, "molecular")
            for base in base_keys:
                if base in result and result[base].level == candidate.level:
                    if base in ("high", "molecular"):
                        raise ValueError(
                            f"The {setting} ({text}) is the {base} level "
                            f"({result[base].level}): set it to '{same}', or "
                            "every monomer would be computed twice."
                        )
                    result[key] = result[base]
                    return
            result[key] = candidate

        result = {"high": make("high", P["high level"], False, "molecular")}
        triples = _triples_level(P)
        if triples:
            # Its own Evaluator (and directory): the triples' high level
            result["high:3"] = make("high_triples", triples, False, "molecular")
            if result["high:3"].level == result["high"].level:
                raise ValueError(
                    f"The triples' high level ({triples}) is the high level "
                    f"({result['high'].level}): set it to 'same as the high level', "
                    "or every monomer and pair would be computed twice."
                )
        if not P["molecular low level"]:
            raise ValueError("The MBE step needs a molecular low level.")
        result["molecular"] = make(
            "molecular", P["molecular low level"], False, "molecular"
        )
        triples_low = _triples_low_level(P)
        if triples_low:
            # Its own Evaluator (and directory): the triples' molecular low level
            result["molecular:3"] = make(
                "molecular_triples", triples_low, False, "molecular"
            )
            if result["molecular:3"].level == result["molecular"].level:
                raise ValueError(
                    f"The triples' low level ({triples_low}) is the molecular low "
                    f"level ({result['molecular'].level}): set it to 'same as the "
                    "molecular low level', or every monomer and pair would be "
                    "computed twice."
                )
        shell = _shell_level(P)
        if shell:
            own(
                "high:shell",
                "high_shells",
                shell,
                ("high", "high:3"),
                "shell fragments' high level",
                "same as the high level",
            )
        shell_low = _shell_low_level(P)
        if shell_low:
            own(
                "molecular:shell",
                "molecular_shells",
                shell_low,
                ("molecular", "molecular:3"),
                "shell fragments' low level",
                "same as the molecular low level",
            )
        if P["periodic low level"] != "none":
            result["periodic"] = make(
                "periodic", P["periodic low level"], False, "periodic"
            )
        periodic_frames = any(f["system"].periodic for f in frames)
        cluster_frames = any(not f["system"].periodic for f in frames)
        cell = P["cell low level"]
        if periodic_frames:
            if cell == "automatic":
                if "periodic" not in result:
                    raise ValueError(
                        "A periodic cell needs a low level for the whole cell: set "
                        "the periodic low level or the whole-system low level."
                    )
                text = P["periodic low level"]
            else:
                text = cell
            result["cell"] = make("cell", text, True, "cell")
            result["cell"].convention = levels.stress_convention(result["cell"].mc)
        if cluster_frames:
            if "periodic" in result:
                raise ValueError(
                    "A periodic low level applies to the fragments of periodic "
                    "cells, but the selection includes a cluster: set the periodic "
                    "low level to 'none' for clusters."
                )
            text = P["molecular low level"] if cell == "automatic" else cell
            result["cluster"] = make("cluster", text, False, "cell")

        # The periodic fragments and the cell must be the same calculation
        # (code, potentials, cutoff, grid) for their errors to cancel.
        if periodic_frames and "periodic" in result:
            if result["cell"].level != result["periodic"].level:
                raise ValueError(
                    f"The cell's low level ({result['cell'].level}) must be the "
                    f"periodic fragments' ({result['periodic'].level}): their "
                    "errors cancel only if they are the same calculation."
                )

        for key in task_options:
            if key not in result:
                raise ValueError(
                    f"The task options are for the {key} level, which this step "
                    f"does not use. With these settings the levels are "
                    f"{', '.join(sorted(result))}."
                )
            add_options(key)

        self._aux_levels(P, result)

        if P["counterpoise"] != "none":
            if periodic_frames:
                raise ValueError(
                    "Counterpoise is not available for periodic cells yet: their "
                    "stress would need the counterpoise correction too."
                )
            corrected = ["high"]
            if P["counterpoise"] == "pairwise, both levels":
                corrected.append("molecular")
            for name in corrected:
                options = result[name].mc.get("options") or {}
                if options.get("mdi_capable") and not options.get("prefers_batch"):
                    raise ValueError(
                        f"Counterpoise needs ghost atoms, which "
                        f"{result[name].level} cannot have here (it runs through "
                        "an MDI engine). Use a program that runs calculations as "
                        "tasks, such as ORCA."
                    )
        return result

    @staticmethod
    def _aux_levels(P, result):
        """Add to ``result`` the copies of the levels with explicit auxiliary
        basis sets, "<level>@aux", for the fragments listed to use them. Levels
        computed together (one Level) share one copy."""
        listed = _aux_listed(P)
        if not listed:
            return
        if P["counterpoise"] != "none":
            raise ValueError(
                "Explicit auxiliary basis sets for single fragments are not "
                "available with counterpoise yet: the counterpoise calculations of "
                "those pairs would need them too."
            )
        sets = levels.parse_aux_sets(P["auxiliary basis sets"])
        for level in sets:
            if level not in listed:
                raise ValueError(
                    f"Auxiliary basis sets are given for the {level} level, which has "
                    "no auxiliary basis fragments."
                )
        copies = {}  # id(Level): (copy, sets)
        for level in listed:
            if level not in result or level in ("periodic", "cell", "cluster"):
                raise ValueError(
                    f"The auxiliary basis fragments are for the {level} level, which "
                    "this step does not compute fragments at."
                )
            base = result[level]
            options = base.mc.get("options") or {}
            if base.mc.get("step") != "ORCA" or (
                options.get("mdi_capable") and not options.get("prefers_batch")
            ):
                raise ValueError(
                    f"The {level} level ({base.level}) cannot take explicit auxiliary "
                    "basis sets: they are an option of ORCA's calculations run as "
                    "tasks."
                )
            aux = sets.get(level, "explicit")
            if id(base) in copies:
                copy, other = copies[id(base)]
                if other != aux:
                    raise ValueError(
                        f"The levels computed together with {level} are given two "
                        f"different auxiliary basis sets: {other} and {aux}."
                    )
            else:
                copy = base.variant(f"{base.name}_aux", aux_basis=aux)
                copies[id(base)] = (copy, aux)
            result[f"{level}@aux"] = copy

    # ------------------------------------------------------------ evaluation
    def _evaluate(self, levels_, frames):
        """Run every level on what each frame needs.

        The levels that run as queued tasks run at the same time, so the slow
        periodic calculations (the cell and the periodic fragments, each hours
        on many cores) overlap the many molecular ones instead of following
        them; they are started first. Levels driven through an MDI engine run
        one at a time in this process.

        Returns {level name: {key: EvaluatorResult}}.
        """
        work = []
        aliases = {}  # name: the name it shares one Level (and its results) with
        first_name = {}
        for name, level in self._evaluation_order(levels_):
            if id(level) in first_name:
                aliases[name] = first_name[id(level)]
            else:
                first_name[id(level)] = name
        for name, level in self._evaluation_order(levels_):
            structures = {}
            for frame in frames:
                system = frame["system"]
                prefix = frame["prefix"]
                if name in ("cell", "cluster"):
                    if (name == "cell") == system.periodic:
                        whole = levels.whole(system)
                        if name == "cell" and "periodic" in levels_:
                            grid = {"max_spacing": self._grid["max_spacing"]}
                            whole = (whole, {"grid": grid})
                        structures[prefix + "cell"] = whole
                    continue
                fragments = frame["fragments"]
                options = None
                if name == "periodic" and system.periodic:
                    options = {"grid": {**self._grid, "reference_cell": system.cell}}
                for fname in frame["calculations"].get(name, []):
                    geometry = levels.geometry(system, fragments[fname])
                    if options is not None:
                        geometry = (geometry, options)
                    structures[prefix + fname] = geometry
                if self._counterpoise and (
                    name == "high" or (name == "molecular" and self._counterpoise_low)
                ):
                    for pair in fragments.by_order(2, in_sum=True):
                        if name == "molecular" and pair.level != "molecular":
                            continue
                        for label, job in counterpoise.ghost_jobs(system, pair).items():
                            structures[counterpoise.key(prefix, pair.name, label)] = job
            if not structures:
                continue
            if name in aliases:
                # The same calculations as another level's: run them once, there
                target = aliases[name]
                for entry in work:
                    if entry[0] == target:
                        entry[2].update(structures)
                        break
                else:
                    work.append((target, level, structures))
                continue
            work.append((name, level, structures))

        for name, level, structures in work:
            printer.important(
                __(
                    f"Running {len(structures)} calculations at the {name} level, "
                    f"{level.level}.",
                    indent=4 * " ",
                )
            )
        results = self._run_levels(work)
        for name, target in aliases.items():
            if target in results:
                results[name] = results[target]

        for name, level, structures in work:
            if name == "cell":
                for result in results[name].values():
                    if result.ok and result.stress is None:
                        result.ok = False
                        result.reason = f"{level.level} returned no stress for the cell"
            done = results[name].values()
            restored = sum(1 for r in done if getattr(r, "restored", False))
            failed = sum(1 for r in done if not r.ok)
            notes = []
            if restored:
                notes.append(f"{restored} finished in an earlier run")
            if failed:
                notes.append(f"{failed} failed")
            if notes:
                printer.important(
                    __(f"    ({name} level: " + "; ".join(notes) + ")", indent=4 * " ")
                )
        return results

    @staticmethod
    def _evaluation_order(levels_):
        """The levels with the long-running ones first: the whole cell (or
        cluster), then the periodic fragments, then the rest as given."""
        first = [n for n in ("cell", "cluster", "periodic") if n in levels_]
        rest = [n for n in levels_ if n not in first]
        return [(n, levels_[n]) for n in first + rest]

    def _run_levels(self, work):
        """Evaluate each level: those that run as queued tasks concurrently, in
        threads, and those that use an MDI engine one after the other here.

        A level that raises, or whose every calculation failed, means no
        configuration can be labelled, so the other levels are cancelled at
        once rather than left to run for hours: their calculations come back as
        failed ("cancelled"), and a rerun of the job submits them afresh. An
        exception is raised again once the levels have stopped.

        Returns {level name: {key: EvaluatorResult}}.
        """
        from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait

        results = {}
        batch = [w for w in work if w[1].runs_as_tasks(self)]
        serial = [w for w in work if w not in batch]
        if len(batch) < 2:
            serial = work
            batch = []
        failure = []  # (level name, exception or None)

        def failed(name, error=None):
            if not failure:
                failure.append((name, error))
                for other, level, _ in work:
                    if other != name:
                        level.cancel()

        def hopeless(level_results):
            return bool(level_results) and not any(r.ok for r in level_results.values())

        with ThreadPoolExecutor(max_workers=max(1, len(batch))) as pool:
            futures = {
                pool.submit(
                    level.evaluate, self, structures, stress=(name == "cell")
                ): name
                for name, level, structures in batch
            }
            collected = set()

            def collect(future):
                """Take a finished batch level's results, noting a failure."""
                collected.add(future)
                name = futures[future]
                error = future.exception()
                if error is not None:
                    failed(name, error)
                    return
                results[name] = future.result()
                if hopeless(results[name]):
                    failed(name)

            for name, level, structures in serial:
                # A batch level that has already failed stops the rest here too
                for future in [f for f in futures if f.done() and f not in collected]:
                    collect(future)
                if failure:
                    break
                try:
                    results[name] = level.evaluate(
                        self, structures, stress=(name == "cell")
                    )
                except Exception as error:
                    failed(name, error)
                    break
                if hopeless(results[name]):
                    failed(name)
            pending = set(futures) - collected
            while pending:
                done, pending = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    collect(future)
        if failure:
            name, error = failure[0]
            if error is not None:
                raise error
            if len(work) > 1:
                printer.important(
                    __(
                        f"Every calculation at the {name} level failed, so the other "
                        "levels were stopped early.",
                        indent=4 * " ",
                    )
                )
        return results

    def _assemble(self, frame, levels_, results, offsets):
        """The labels of one frame, or MissingFragmentsError."""
        system = frame["system"]
        prefix = frame["prefix"]
        n = len(prefix)

        def library(name):
            out = {}
            for key, result in results.get(name, {}).items():
                if key.startswith(prefix) and result.ok and "-cp-" not in key:
                    out[key[n:]] = levels.to_library(result)
            return out

        cell_name = "cell" if system.periodic else "cluster"
        cell_result = results.get(cell_name, {}).get(prefix + "cell")
        if cell_result is None or not cell_result.ok:
            raise seamm_mbe.MissingFragmentsError([("whole system", cell_name)])
        if system.periodic:
            cell = levels.to_library(
                cell_result,
                volume=system.volume,
                convention=levels_["cell"].convention,
            )
        else:
            cell = levels.to_library(cell_result)
        corrections = self._corrections(frame, results) if self._counterpoise else {}
        high_by_order, low_by_order = {}, {}
        for key, name in ((3, "high:3"), (seamm_mbe.SHELL, "high:shell")):
            if name in levels_:
                high_by_order[key] = library(name)
        for key, name in ((3, "molecular:3"), (seamm_mbe.SHELL, "molecular:shell")):
            if name in levels_:
                low_by_order[key] = library(name)
        high_by_order = high_by_order or None
        low_by_order = low_by_order or None
        variants = frame.get("variants") or {}
        correction = seamm_mbe.mbe_correction(
            frame["fragments"],
            library("high"),
            periodic=library("periodic"),
            molecular=library("molecular"),
            corrections=corrections,
            high_by_order=high_by_order,
            low_by_order=low_by_order,
            variants=variants or None,
            variant_results={v: library(v) for v, _ in variants.values()} or None,
        )
        frame["correction"] = correction
        term = seamm_mbe.CellTerm(
            levels_[cell_name].level, cell["energy"], cell["forces"], cell.get("virial")
        )
        return seamm_mbe.assemble(system, correction, [term], offsets=offsets)

    def _corrections(self, frame, results):
        """The pairwise counterpoise corrections of a frame's selected pairs:
        {pair name: (eV, eV/Å)}, at the high level and, with "pairwise, both
        levels" for pairs referenced to it, the molecular low level. A missing
        calculation makes the frame incomplete."""
        system = frame["system"]
        prefix = frame["prefix"]
        corrections = {}
        missing = []
        frame["counterpoise fallbacks"] = 0
        for pair in frame["fragments"].by_order(2, in_sum=True):
            total_e = 0.0
            total_f = np.zeros((pair.n_atoms, 3))
            for name, sign in (("high", 1.0), ("molecular", -1.0)):
                if name == "molecular" and (
                    not self._counterpoise_low or pair.level != "molecular"
                ):
                    continue
                try:
                    d_e, d_f, fallback = counterpoise.correction(
                        system, frame["fragments"], pair, results[name], prefix
                    )
                except (KeyError, ValueError) as e:
                    reason = e.args[0] if e.args else str(e)
                    missing.append((f"counterpoise of {pair.name}: {reason}", name))
                    continue
                frame["counterpoise fallbacks"] += int(fallback)
                total_e += sign * d_e
                total_f += sign * d_f
            corrections[pair.name] = (total_e, total_f)
        frame["counterpoise correction"] = sum(e for e, _ in corrections.values())
        if missing:
            raise seamm_mbe.MissingFragmentsError(missing)
        return corrections

    @staticmethod
    def _failure_text(frame, error, results):
        """The missing calculations of a frame, with reasons where known."""
        prefix = frame["prefix"]
        parts = []
        for name, level in error.missing[:5]:
            if name.startswith("counterpoise of "):
                parts.append(f"{name} ({level})")
                continue
            key = prefix + ("cell" if name == "whole system" else name)
            result = results.get(level, {}).get(key)
            reason = result.reason if result is not None else "not run"
            parts.append(f"{name} ({level}): {reason}")
        more = len(error.missing) - 5
        return "; ".join(parts) + (f"; and {more} more" if more > 0 else "")

    # ----------------------------------------------------------------- output
    def _store(self, frame, labels, levels_, P):
        """Store the labels as properties, on the atoms, via store_results and
        in the extxyz file."""
        configuration = frame["configuration"]
        system = frame["system"]
        data = labels_.to_seamm(labels, system.volume)
        data["configuration name"] = configuration.name
        data["model chemistry"] = self.model
        data["fragment counts"] = {
            str(k): v for k, v in frame["fragments"].counts().items()
        }
        if self._counterpoise:
            data["counterpoise correction"] = float(
                seamm_mbe.units.ev_to_kj_per_mol(frame["counterpoise correction"])
            )
            data["counterpoise fallbacks"] = frame.get("counterpoise fallbacks", 0)

        configuration.atoms.set_gradients(
            np.asarray(data["gradients"]), fractionals=False
        )
        self._put_property(configuration, "energy", data["energy"], "kJ/mol")
        self._put_property(
            configuration, "gradients", data["gradients"], "kJ/mol/Å", "json"
        )
        if "stress" in data:
            self._put_property(configuration, "stress", data["stress"], "GPa", "json")
        self.store_results(configuration=configuration, data=data)

        filename = P["extxyz file"]
        if filename and filename != "none":
            path = self.file_path(filename, relative_to=Path(self.directory))
            cell = levels_["cell" if system.periodic else "cluster"].level
            labels_.write_extxyz(
                path,
                system,
                labels,
                name=configuration.name,
                identifier=configuration.id,
                counterpoise=self._counterpoise,
                model=f"{self.model} MBE on {cell}",
            )
        self._write_increments(frame, P)
        frame["labels data"] = data

    def _write_increments(self, frame, P):
        """Each selected fragment's increment, for later analysis (e.g. a
        correction by triple topology): ``increments_c<configuration id>.json`` in
        the step's directory, replaced on a rerun. Energies in eV, forces in
        eV/Å in the fragment's atom order."""
        correction = frame.get("correction")
        if correction is None:
            return
        fragments = frame["fragments"]
        molecules = fragments.system.molecules
        triple_cutoff = None
        if P["cutoffs"] == "single value" and P["triple rule"] != "none":
            triple_cutoff = P["triple cutoff"].m_as("Å")
        records = []
        for name, inc in correction.increments.items():
            f = fragments[name]
            record = {
                "name": name,
                "order": inc.order,
                "molecules": [int(m) for m in f.molecules],
                "images": [[int(x) for x in image] for image in f.images],
                "distances": {
                    f"{i},{j}": float(d) for (i, j), d in f.distances.items()
                },
                "low level": inc.low_level or inc.level,
                "high level": inc.high_level,
                "energy": float(inc.energy),
                "forces": np.asarray(inc.forces).tolist(),
            }
            members = getattr(fragments.system, "members", None)
            if members is not None:
                # With ion shells the "molecules" are units: the real molecules
                # each holds (a shell's ion first)
                record["members"] = [[int(m) for m in members[u]] for u in f.molecules]
            if inc.order == 3:
                # The cutoffs the enumerator used, per slot pair (tables too)
                n = 0
                for (i, j), d in f.distances.items():
                    a = molecules[f.molecules[i]].type
                    b = molecules[f.molecules[j]].type
                    if d < fragments.rules.cutoff(3, a, b):
                        n += 1
                record["topology"] = "closed" if n == 3 else "chain"
            records.append(record)
        # The pairs of two shells at or beyond the cap: not computed, so no
        # energy or forces, and not in the sum
        cap = fragments.rules.shell_pair_cap
        for c in fragments.capped:
            record = {
                "name": c["name"],
                "order": 2,
                "molecules": [int(m) for m in c["molecules"]],
                "images": [[int(x) for x in image] for image in c["images"]],
                "distances": {"0,1": float(c["distance"])},
                "capped": f"{cap:g}A",
                "computed": False,
            }
            members = getattr(fragments.system, "members", None)
            if members is not None:
                record["members"] = [
                    [int(m) for m in members[u]] for u in c["molecules"]
                ]
            records.append(record)
        path = Path(self.directory) / f"increments_c{frame['configuration'].id}.json"
        path.write_text(
            json.dumps(
                {
                    "configuration": frame["configuration"].name,
                    "distance criterion": P["distance criterion"],
                    "triple cutoff (Å)": triple_cutoff,
                    "shell pair cap (Å)": cap,
                    "units": {"energy": "eV", "forces": "eV/Å", "distances": "Å"},
                    "increments": records,
                },
                indent=1,
            )
        )

    def _put_property(self, configuration, key, value, units, _type="float"):
        """Store a ``<key>#MBE#<model>`` property, defining it if needed."""
        name = self.metadata["results"][key]["property"].format(model=self.model)
        properties = configuration.properties
        if not properties.exists(name):
            description = self.metadata["results"][key]["description"]
            properties.add(
                name,
                _type,
                units=units,
                description=f"{description} ({self.model})",
                noerror=True,
            )
        properties.put(name, value)

    def _cite(self):
        """Cite the many-body expansion and the seamm_mbe library. The codes of
        the levels will be cited through the providers once the batch contract
        has a citation hook."""
        for key, note in (
            ("Herbert_2019", "A review of fragment-based quantum chemistry."),
            ("Richard_2014", "The many-body expansion for benchmark accuracy."),
            ("Beran_2016", "Fragment methods for molecular crystals."),
            ("seamm_mbe", "The MBE bookkeeping library."),
        ):
            if key in self._bibliography:
                self.references.cite(
                    raw=self._bibliography[key],
                    alias=key,
                    module="mbe_step",
                    level=1 if key != "seamm_mbe" else 2,
                    note=note,
                )

    def analyze(self, indent="", rows=(), levels=None, elapsed=0.0, **kwargs):
        """Print the results of each configuration."""
        if levels:
            lines = ["The levels:"]
            for name, level in levels.items():
                lines.append(f"    {name:10s} {level.level} ({level.mc['step']})")
            printer.important(__("\n".join(lines), indent=4 * " ", wrap=False))
            printer.important("")
        for frame, labels in rows:
            configuration = frame["configuration"]
            data = frame["labels data"]
            fragments = frame["fragments"]
            calcs = frame["calculations"]
            ghosts = {"high": 0, "molecular": 0}
            if self._counterpoise:
                for pair in fragments.by_order(2, in_sum=True):
                    ghosts["high"] += 2
                    if self._counterpoise_low and pair.level == "molecular":
                        ghosts["molecular"] += 2
            n = len(frame["system"].molecules)
            types = ", ".join(
                f"{c} {t}" for t, c in frame["system"].type_counts().items()
            )
            text = [f"{configuration.name}: {n} molecules ({types})"]
            shells = getattr(fragments.system, "shell_sizes", dict)()
            if shells:
                sizes = sorted(shells.values())
                span = (
                    f"{sizes[0]}"
                    if sizes[0] == sizes[-1]
                    else f"{sizes[0]}-{sizes[-1]}"
                )
                text.append(
                    f"    ion shells: {len(sizes)}, of {span} molecules besides the ion"
                )
            text += [
                f"    fragments: {labels_.counts_text(fragments)}",
                f"    calculations: {len(calcs['high']) + ghosts['high']} high, "
                + (
                    f"{len(calcs['high:3'])} high (triples), "
                    if "high:3" in calcs
                    else ""
                )
                + (
                    f"{len(calcs['high:shell'])} high (shells), "
                    if "high:shell" in calcs
                    else ""
                )
                + f"{len(calcs['molecular']) + ghosts['molecular']} molecular, "
                + (
                    f"{len(calcs['molecular:3'])} molecular (triples), "
                    if "molecular:3" in calcs
                    else ""
                )
                + (
                    f"{len(calcs['molecular:shell'])} molecular (shells), "
                    if "molecular:shell" in calcs
                    else ""
                )
                + f"{len(calcs['periodic'])} periodic, plus the whole system",
                f"    E = {data['energy']:.3f} kJ/mol, correction "
                f"{data['MBE energy']:.3f} kJ/mol",
            ]
            correction = frame.get("correction")
            if correction is not None and fragments.shells:
                by_order = {}
                for inc in correction.increments.values():
                    if inc.high_level == "high:shell" or inc.low_level == (
                        "molecular:shell"
                    ):
                        by_order[inc.order] = by_order.get(inc.order, 0) + 1
                if by_order:
                    text.append(
                        "    increments at the shell levels: "
                        + ", ".join(
                            f"{c} of order {k}" for k, c in sorted(by_order.items())
                        )
                    )
            per = labels.per_body
            text.append(
                "    per body (kJ/mol per molecule): "
                + ", ".join(
                    f"{k}-body {data['per-body energies'][str(k)] / n:+.3f}"
                    for k in sorted(per)
                )
            )
            if labels.pressure is not None:
                text.append(
                    f"    P = {labels.pressure:+.1f} atm atomic, "
                    f"{labels.molecular_pressure:+.1f} atm molecular "
                    f"(MBE {labels.breakdown['MBE']['pressure']:+.1f} atm atomic, "
                    f"{labels.breakdown['MBE']['molecular pressure']:+.1f} atm "
                    "molecular)"
                )
                text.append(
                    "    per body (atm, atomic): "
                    + ", ".join(
                        f"{k}-body {per[k]['pressure']:+.1f}" for k in sorted(per)
                    )
                )
                text.append(
                    "    per body (atm, molecular): "
                    + ", ".join(
                        f"{k}-body {per[k]['molecular pressure']:+.1f}"
                        for k in sorted(per)
                    )
                )
            text.append(
                "    largest increment net force "
                f"{labels.max_increment_net_force * 1000:.2f} meV/Å"
            )
            if self._counterpoise:
                n_pairs = len(frame["fragments"].by_order(2, in_sum=True))
                text.append(
                    f"    counterpoise: {n_pairs} pairs corrected "
                    f"({'both levels' if self._counterpoise_low else 'high level'}), "
                    f"{data['counterpoise correction']:+.3f} kJ/mol in all"
                    + (
                        f"; {frame['counterpoise fallbacks']} kept the uncorrected "
                        "gradient (unphysical ghost gradients)"
                        if frame.get("counterpoise fallbacks")
                        else ""
                    )
                )
            printer.important(__("\n".join(text), indent=4 * " ", wrap=False))
            printer.important("")
        if rows:
            printer.important(
                __(f"The calculations took {elapsed:.1f} s.", indent=4 * " ")
            )
            printer.important("")
