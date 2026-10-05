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


def _truthy(value):
    return value is True or (isinstance(value, str) and value.lower() == "yes")


def _triples_level(P):
    """The triples' own high level, or None when it is the high level's."""
    text = (P["triples high level"] or "").strip()
    if text in ("", "same as the high level") or text == P["high level"]:
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
        triples = (P["triples high level"] or "").strip()
        if order >= 3 and triples not in ("", "same as the high level"):
            text += (
                f" The triples' increments use {triples} as their high level, with "
                "their pairs and monomers at that level too."
            )
        if order >= 2:
            if P["cutoffs"] == "single value":
                text += f" Pairs closer than {P['pair cutoff']} are selected"
            else:
                text += " Pairs are selected with cutoffs by molecule type"
            if order >= 3 and P["triple rule"] != "none":
                text += f", and {P['triple rule']} triples"
            text += f" ({P['distance criterion']})."
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
        offsets = labels_.parse_offsets(P["energy offsets"])
        frames = [self._frame(c, rules, P) for c in configurations]
        levels_ = self._levels(P, context, frames)
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
        return seamm_mbe.SelectionRules(
            max_order=order,
            criterion=CRITERIA[P["distance criterion"]],
            cutoffs=cutoffs,
            rules=rules,
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
        fragments = seamm_mbe.enumerate_fragments(system, rules)
        periodic = {}
        if P["periodic low level"] != "none":
            periodic[1] = _truthy(P["periodic monomers"])
            r1 = P["periodic pair cutoff"].m_as("Å")
            if r1 > 0:
                periodic[2] = r1
        seamm_mbe.assign_levels(fragments, periodic)
        self._check_frame(configuration, system, rules, P)
        high_levels = {3: "high:3"} if _triples_level(P) else None
        return {
            "configuration": configuration,
            "system": system,
            "fragments": fragments,
            "calculations": fragments.calculations(high_levels=high_levels),
            "prefix": f"c{configuration.id}-",
        }

    def _check_frame(self, configuration, system, rules, P):
        """Refuse, before any calculation, what would fail or be wrong later:
        missing energy offsets, a configuration charge that disagrees with the
        molecules', and open-shell molecules inside larger fragments."""
        name = configuration.name
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
        open_shell = sorted(
            t.name for t in system.types.values() if t.multiplicity != 1
        )
        if open_shell and rules.max_order >= 2:
            raise ValueError(
                f"Configuration {name} contains open-shell molecules "
                f"({', '.join(open_shell)}); fragments of several molecules with "
                "open shells are not supported yet. Use maximum order 1."
            )

    def _levels(self, P, context, frames):
        """Resolve each level's model chemistry, checking what is needed."""
        current = (
            self.get_variable("_model_chemistry")
            if self.variable_exists("_model_chemistry")
            else None
        )
        hours = P["bundle walltime"].m_as("s")
        archive = _truthy(P["archive"])

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
            )

        result = {"high": make("high", P["high level"], False, "molecular")}
        triples = _triples_level(P)
        if triples:
            # Its own Evaluator (and directory): the triples' high level
            result["high:3"] = make("high_triples", triples, False, "molecular")
        if not P["molecular low level"]:
            raise ValueError("The MBE step needs a molecular low level.")
        result["molecular"] = make(
            "molecular", P["molecular low level"], False, "molecular"
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

    # ------------------------------------------------------------ evaluation
    def _evaluate(self, levels_, frames):
        """Run every level, one after the other, on what each frame needs.

        Returns {level name: {key: EvaluatorResult}}.
        """
        results = {}
        for name, level in levels_.items():
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
                for fname in frame["calculations"][name]:
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
            printer.important(
                __(
                    f"Running {len(structures)} calculations at the {name} level, "
                    f"{level.level}.",
                    indent=4 * " ",
                )
            )
            results[name] = level.evaluate(self, structures, stress=(name == "cell"))
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
                printer.important(__("    (" + "; ".join(notes) + ")", indent=4 * " "))
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
        high_by_order = {3: library("high:3")} if "high:3" in levels_ else None
        correction = seamm_mbe.mbe_correction(
            frame["fragments"],
            library("high"),
            periodic=library("periodic"),
            molecular=library("molecular"),
            corrections=corrections,
            high_by_order=high_by_order,
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
                "low level": inc.level,
                "high level": inc.high_level,
                "energy": float(inc.energy),
                "forces": np.asarray(inc.forces).tolist(),
            }
            if inc.order == 3 and triple_cutoff is not None:
                n = sum(1 for d in f.distances.values() if d < triple_cutoff)
                record["topology"] = "closed" if n == 3 else "chain"
            records.append(record)
        path = Path(self.directory) / f"increments_c{frame['configuration'].id}.json"
        path.write_text(
            json.dumps(
                {
                    "configuration": frame["configuration"].name,
                    "distance criterion": P["distance criterion"],
                    "triple cutoff (Å)": triple_cutoff,
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
            text = [
                f"{configuration.name}: {n} molecules ({types})",
                f"    fragments: {labels_.counts_text(fragments)}",
                f"    calculations: {len(calcs['high']) + ghosts['high']} high, "
                + (
                    f"{len(calcs['high:3'])} high (triples), "
                    if "high:3" in calcs
                    else ""
                )
                + f"{len(calcs['molecular']) + ghosts['molecular']} molecular, "
                f"{len(calcs['periodic'])} periodic, plus the whole system",
                f"    E = {data['energy']:.3f} kJ/mol, correction "
                f"{data['MBE energy']:.3f} kJ/mol",
            ]
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
                    f"(MBE {labels.breakdown['MBE']['pressure']:+.1f} atm)"
                )
                text.append(
                    "    per body (atm): "
                    + ", ".join(
                        f"{k}-body {per[k]['pressure']:+.1f}" for k in sorted(per)
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
