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
    """

    def __init__(self, name, mc, *, ranks, memory, bundle, walltime, archive):
        self.name = name
        self.mc = mc
        self.ranks = int(ranks)
        self.memory = float(memory)
        self.bundle = int(bundle)
        self.walltime = float(walltime)
        self.archive = bool(archive)
        self._evaluator = None  # while evaluate() runs, so cancel() can reach it
        self._cancel_requested = False

    @property
    def level(self):
        return self.mc["level"]

    def resources(self):
        from seamm_exec import Resources

        return Resources(ntasks=self.ranks, mem_per_cpu=int(self.memory * 1e6))

    def cancel(self):
        """Stop this level's calculations, from another thread: those in flight
        are killed and the rest dropped, all coming back as failed with the
        reason "cancelled" (a rerun submits them afresh). Results already
        produced stand."""
        self._cancel_requested = True
        evaluator = self._evaluator
        if evaluator is not None:
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
        """
        import seamm_exec

        properties = ("energy", "gradients") + (("stress",) if stress else ())
        options = {"archive": self.archive, "bundle_tasks": self.bundle}
        if self.walltime > 0:
            options["bundle_walltime"] = self.walltime
        results = {}
        with seamm_exec.Evaluator(
            node,
            self.mc,
            properties=properties,
            directory=Path(node.directory) / self.name,
            task_set_options=options,
            resources=self.resources(),
            name=f"MBE_{self.name}",  # no spaces: it goes into MDI_Init
        ) as evaluator:
            self._evaluator = evaluator
            if self._cancel_requested:
                evaluator.cancel()
            try:
                for key, structure in structures.items():
                    options = None
                    if isinstance(structure, tuple):
                        structure, options = structure
                    evaluator.submit(structure, key=key, options=options)
                for result in evaluator.results():
                    results[result.key] = result
            finally:
                self._evaluator = None
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
