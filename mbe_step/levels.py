# -*- coding: utf-8 -*-

"""The MBE step's levels of theory: resolving their model chemistries and
evaluating structures with each through seamm_exec's Evaluator.

The step uses up to four levels:

* high: every fragment (by default the flowchart's current Model Chemistry);
* molecular: the molecular low level, on the fragments referenced to it;
* periodic: the periodic low level, on the compact fragments referenced to it;
* cell: the low level of the whole cell (or cluster).

Each level is one Evaluator with its own directory, so each keeps its own task
manifest (restart) and bundles; the Evaluator chooses how the level runs
(batch tasks or an MDI engine). In this phase the levels run one after the
other; running them concurrently is a seamm_exec follow-up.
"""

import logging
from pathlib import Path

import numpy as np

import seamm_mbe

logger = logging.getLogger(__name__)

#: The sentinel for the flowchart's current Model Chemistry
CURRENT = "current model chemistry"


def resolve(text, context, current=None, periodic=False):
    """The ``_model_chemistry`` wrapper for a level given as text.

    Parameters
    ----------
    text : str
        A level spec such as ``"ORCA:DFT@r2SCAN-D4/def2-TZVPPD"``, possibly with
        ``$variable`` components, or :data:`CURRENT`.
    context : dict
        The flowchart variables, to dereference any ``$variable``.
    current : dict or None
        The flowchart's ``_model_chemistry``, for :data:`CURRENT`.
    periodic : bool
        Whether the level must handle periodic systems.

    Returns
    -------
    dict
        The ``_model_chemistry`` wrapper: level, step, options, ...
    """
    from model_chemistry_step.model_chemistry import (
        availability_problem,
        discover_model_chemistries,
        match_model_chemistry,
        resolve_level,
    )

    if text == CURRENT:
        if current is None:
            raise ValueError(
                "The MBE step's high level is the current model chemistry, but "
                "there is none: add a Model Chemistry step before the MBE step, "
                "or give the high level explicitly."
            )
        return current
    text = resolve_level(text, context)
    available = discover_model_chemistries(periodic_only=periodic)
    if text in available:
        return available[text]
    mc = match_model_chemistry(text, available)
    if mc is None:
        raise ValueError(availability_problem(text, available, periodic=periodic))
    return mc


def stress_convention(mc):
    """The sign convention of a level's stress, as its provider declares it
    ("pressure" or "stress"); refuse a periodic level that does not."""
    convention = (mc.get("options") or {}).get("stress_convention")
    if convention not in seamm_mbe.units.CONVENTIONS:
        raise ValueError(
            f"The model chemistry '{mc['level']}' does not declare the sign "
            "convention of its stress (options['stress_convention'] = 'pressure' "
            "or 'stress'), so its stress cannot be used for a periodic cell."
        )
    return convention


def parse_resources(text):
    """Larger resources for larger fragments, from entries "atoms ranks memory"
    separated by ';' or new lines: fragments of at least that many atoms get
    that many MPI ranks and that much memory per rank (MB). Blank or "none" for
    none.

    Returns
    -------
    [(int, int, float)]
        (atoms, ranks, memory) ascending in atoms.
    """
    text = (text or "").strip()
    if text.lower() in ("", "none"):
        return []
    tiers = []
    for entry in text.replace("\n", ";").split(";"):
        entry = entry.strip()
        if not entry:
            continue
        fields = entry.split()
        if len(fields) != 3:
            raise ValueError(
                f"The fragment resources '{entry}' should be 'atoms ranks memory', "
                "e.g. '60 16 8000' for 16 ranks of 8,000 MB from 60 atoms."
            )
        try:
            atoms, ranks, memory = int(fields[0]), int(fields[1]), float(fields[2])
        except ValueError:
            raise ValueError(
                f"The fragment resources '{entry}' should be whole numbers of atoms "
                "and ranks and a memory in MB."
            ) from None
        if atoms < 1 or ranks < 1 or memory <= 0:
            raise ValueError(f"The fragment resources '{entry}' must be positive.")
        tiers.append((atoms, ranks, memory))
    atoms = [t[0] for t in tiers]
    if len(set(atoms)) != len(atoms):
        raise ValueError(f"The fragment resources give a size twice: '{text}'.")
    return sorted(tiers)


#: The level names "task options" entries may use, as the step names its levels
LEVEL_NAMES = (
    "high",
    "high:3",
    "high:shell",
    "molecular",
    "molecular:3",
    "molecular:shell",
    "periodic",
    "cell",
    "cluster",
)


def parse_task_options(text):
    """Options for every task of a level, from entries "level key=value ..."
    separated by ';' or new lines, e.g. "high:shell exact_exchange_elements=none".
    The values are text, as the program reads them. Blank or "none" for none.

    Returns
    -------
    {str: {str: str}}
        By level name (:data:`LEVEL_NAMES`).
    """
    text = (text or "").strip()
    if text.lower() in ("", "none"):
        return {}
    result = {}
    for entry in text.replace("\n", ";").split(";"):
        entry = entry.strip()
        if not entry:
            continue
        fields = entry.split()
        name, pairs = fields[0], fields[1:]
        if name not in LEVEL_NAMES:
            raise ValueError(
                f"The task options '{entry}' name the level '{name}', which is not "
                f"one of {', '.join(LEVEL_NAMES)}."
            )
        if not pairs or any("=" not in pair for pair in pairs):
            raise ValueError(
                f"The task options '{entry}' should be 'level key=value ...', e.g. "
                "'high:shell exact_exchange_elements=none'."
            )
        options = result.setdefault(name, {})
        for pair in pairs:
            key, value = pair.split("=", 1)
            if not key:
                raise ValueError(f"The task options '{entry}' have an empty key.")
            options[key] = value
    return result


def n_atoms(structure):
    """The number of atoms of a structure an Evaluator takes (a Geometry or a
    configuration, possibly with options as a tuple)."""
    if isinstance(structure, tuple):
        structure = structure[0]
    n = getattr(structure, "n_atoms", 0)
    return int(n() if callable(n) else n)


class Level:
    """One level of theory and how it runs.

    Parameters
    ----------
    name : str
        "high", "molecular", "periodic" or "cell"; also the directory name.
    mc : dict
        The ``_model_chemistry`` wrapper.
    ranks : int
        MPI ranks per calculation.
    memory : float
        Memory per rank, in MB.
    bundle : int
        Calculations per batch job on a queue.
    walltime : float
        Seconds per batch job.
    archive : bool
        Tar finished bundles.
    tiers : [(int, int, float)] or None
        Larger resources for larger fragments (:func:`parse_resources`): each
        size runs as its own Evaluator, in the directory ``<name>_<atoms>``, at
        the same time as the others.
    options : {str: str} or None
        Options added to every task of this level (:func:`parse_task_options`);
        a structure's own options take precedence.
    """

    def __init__(
        self,
        name,
        mc,
        *,
        ranks,
        memory,
        bundle,
        walltime,
        archive,
        tiers=None,
        options=None,
    ):
        self.name = name
        self.mc = mc
        self.ranks = int(ranks)
        self.memory = float(memory)
        self.bundle = int(bundle)
        self.walltime = float(walltime)
        self.archive = bool(archive)
        self.tiers = list(tiers or [])
        self.options = dict(options or {})
        self._evaluators = []  # while evaluate() runs, so cancel() can reach them
        self._cancel_requested = False

    @property
    def level(self):
        return self.mc["level"]

    def resources(self, ranks=None, memory=None):
        from seamm_exec import Resources

        ranks = self.ranks if ranks is None else ranks
        memory = self.memory if memory is None else memory
        return Resources(ntasks=int(ranks), mem_per_cpu=int(memory * 1e6))

    def split(self, structures):
        """The structures by size tier: [(directory name, ranks, memory,
        {key: structure})], the base tier first, empty tiers left out."""
        parts = [(self.name, self.ranks, self.memory, {})]
        for atoms, ranks, memory in self.tiers:
            parts.append((f"{self.name}_{atoms}", ranks, memory, {}))
        for key, structure in structures.items():
            n = n_atoms(structure)
            index = 0
            for k, (atoms, _, _) in enumerate(self.tiers):
                if n >= atoms:
                    index = k + 1
            parts[index][3][key] = structure
        return [part for part in parts if part[3]]

    def cancel(self):
        """Stop this level's calculations, from another thread: those in flight
        are killed and the rest dropped, all coming back as failed with the
        reason "cancelled" (a rerun submits them afresh). Results already
        produced stand."""
        self._cancel_requested = True
        for evaluator in list(self._evaluators):
            evaluator.cancel()

    def runs_as_tasks(self, node):
        """Whether this level runs its calculations as queued tasks (the batch
        path) rather than through an MDI engine. Batch levels can run at the same
        time as each other; an MDI engine is driven from this process, so those
        levels run one at a time."""
        import seamm_exec

        try:
            evaluator = seamm_exec.Evaluator(
                node,
                self.mc,
                directory=Path(node.directory) / self.name,
                name=f"MBE_{self.name}",
            )
            return evaluator.path == "batch"
        except Exception:
            return False

    def evaluate(self, node, structures, *, stress=False):
        """Evaluate structures at this level.

        Parameters
        ----------
        node : seamm.Node
            The MBE step.
        structures : {str: object}
            Configurations or ``seamm_exec.Geometry`` objects by key.
        stress : bool
            Request the stress (for a periodic cell).

        Returns
        -------
        {str: seamm_exec.EvaluatorResult}

        With size tiers (``tiers``) each tier is its own Evaluator, with its
        resources and directory; they run at the same time.
        """
        parts = self.split(structures)
        if len(parts) <= 1:
            name, ranks, memory, part = (
                parts[0] if parts else (self.name, None, None, {})
            )
            return self._evaluate(node, part, name, ranks, memory, stress)
        from concurrent.futures import ThreadPoolExecutor

        results = {}
        with ThreadPoolExecutor(max_workers=len(parts)) as pool:
            futures = [
                pool.submit(self._evaluate, node, part, name, ranks, memory, stress)
                for name, ranks, memory, part in parts
            ]
            for future in futures:
                results.update(future.result())
        return results

    def _evaluate(self, node, structures, directory, ranks, memory, stress):
        """One Evaluator on ``structures``, in ``directory`` under the step's."""
        import seamm_exec

        properties = ("energy", "gradients") + (("stress",) if stress else ())
        options = {"archive": self.archive, "bundle_tasks": self.bundle}
        if self.walltime > 0:
            options["bundle_walltime"] = self.walltime
        results = {}
        if not structures:
            return results
        with seamm_exec.Evaluator(
            node,
            self.mc,
            properties=properties,
            directory=Path(node.directory) / directory,
            task_set_options=options,
            resources=self.resources(ranks, memory),
            name=f"MBE_{directory}",  # no spaces: it goes into MDI_Init
        ) as evaluator:
            self._evaluators.append(evaluator)
            if self._cancel_requested:
                evaluator.cancel()
            try:
                for key, structure in structures.items():
                    options = None
                    if isinstance(structure, tuple):
                        structure, options = structure
                    if self.options:
                        options = {**self.options, **(options or {})}
                    evaluator.submit(structure, key=key, options=options)
                for result in evaluator.results():
                    results[result.key] = result
            finally:
                self._evaluators.remove(evaluator)
        return results


def geometry(system, fragment):
    """The isolated fragment as a ``seamm_exec.Geometry``: its atoms in
    fragment order at the fragment's coordinates, with its charge and
    multiplicity."""
    from seamm_exec import Geometry

    return Geometry(
        system.atomic_numbers[fragment.atoms],
        fragment.coordinates,
        charge=fragment.charge,
        multiplicity=fragment.multiplicity,
        name=fragment.name,
    )


def whole(system):
    """The whole cell or cluster as a ``seamm_exec.Geometry`` with the charge
    the molecules add up to (their types' charges, from the formal charges or
    the catalog) -- the same charges the fragments use -- and multiplicity 1
    (open-shell molecules are refused before this)."""
    from seamm_exec import Geometry

    charge = sum(system.types[m.type].charge for m in system.molecules)
    return Geometry(
        system.atomic_numbers,
        system.coordinates,
        charge=charge,
        multiplicity=1,
        cell=system.cell,
        name="whole system",
    )


def to_library(result, *, volume=None, convention=None):
    """An EvaluatorResult in seamm_mbe's units: (energy eV, forces eV/Å), plus
    the virial (eV) when the result has a stress."""
    data = {"energy": result.energy, "gradients": np.asarray(result.gradients)}
    if result.stress is not None and volume is not None:
        data["stress"] = result.stress
    out = seamm_mbe.units.from_analyze_task(
        data, volume=volume, stress_convention=convention
    )
    return out
