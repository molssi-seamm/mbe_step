# -*- coding: utf-8 -*-

"""Pairwise counterpoise for the MBE step, through seamm_bsse.

For a selected pair, the Boys-Bernardi correction needs each monomer in the
pair's basis (the other monomer's atoms as ghosts): two extra calculations per
pair at each molecular level. With the pair itself and the two monomers (both
already computed as fragments), seamm_bsse.combine gives the counterpoise-
corrected pair energy and gradient, guarded against unphysical ghost gradients.

The correction enters the sum only (seamm_mbe's ``corrections``): with
[high - low] for the pair,

    correction = (E_high^CP - E_high) - (E_low^CP - E_low)

the second term only when the pair is referenced to the molecular low level (a
plane-wave periodic level has no basis-set superposition error). Triples still
subtract the uncorrected pairs.
"""

import numpy as np
import seamm_bsse
from seamm_util import Q_

from . import levels

HARTREE_KJ = Q_(1.0, "E_h").m_as("kJ/mol")
EV_KJ = Q_(1.0, "eV").m_as("kJ/mol")
BOHR = Q_(1.0, "bohr").m_as("Å")


def specs(system, pair):
    """The seamm_bsse job specs for a pair (atom indices within the pair)."""
    types = [system.types[system.molecules[m].type] for m in pair.molecules]
    fragments = [
        seamm_bsse.Fragment(
            label=label,
            atom_indices=[int(i) for i in pair.slot_atoms[slot]],
            charge=types[slot].charge,
            multiplicity=types[slot].multiplicity,
        )
        for slot, label in enumerate(("A", "B"))
    ]
    return seamm_bsse.generate_job_specs(fragments)


def ghost_jobs(system, pair):
    """The extra calculations for a pair: {label: (structure, options)} for
    each monomer in the pair's basis."""
    geometry = levels.geometry(system, pair)
    jobs = {}
    for spec in specs(system, pair):
        if spec.kind != seamm_bsse.FRAGMENT_IN_CLUSTER:
            continue
        jobs[spec.label] = (
            geometry,
            {
                "atom_indices": list(spec.atom_indices),
                "ghost_atoms": sorted(spec.ghost_indices),
                "charge": spec.charge,
                "multiplicity": spec.multiplicity,
            },
        )
    return jobs


def key(prefix, pair_name, label):
    """The task key of a ghost job."""
    return f"{prefix}{pair_name}-cp-{label}"


def _job_result(result):
    """An EvaluatorResult (kJ/mol, kJ/mol/Å) as a seamm_bsse.JobResult (E_h,
    E_h/bohr)."""
    gradient = np.asarray(result.gradients, dtype=float) * BOHR / HARTREE_KJ
    return seamm_bsse.JobResult(energy=result.energy / HARTREE_KJ, gradient=gradient)


def correction(system, fragments, pair, results, prefix):
    """The counterpoise correction of one pair at one level.

    Parameters
    ----------
    results : {str: EvaluatorResult}
        The level's results by key (the pair, its monomers and its ghost jobs).

    Returns
    -------
    energy : float
        E^CP - E for the pair, eV.
    forces : numpy.ndarray
        (n_atoms, 3) the change of the pair's forces, eV/Å.
    fallback : bool
        Whether seamm_bsse kept the uncorrected gradient (unphysical ghost
        gradients: the energy is corrected, the forces are not).
    """
    monomers = {slots[0]: name for name, slots in pair.subfragments if len(slots) == 1}
    by_label = {
        "cluster": results[prefix + pair.name],
        "A-alone": results[prefix + monomers[0]],
        "B-alone": results[prefix + monomers[1]],
        "A-in-cluster": results[key(prefix, pair.name, "A-in-cluster")],
        "B-in-cluster": results[key(prefix, pair.name, "B-in-cluster")],
    }
    for label, result in by_label.items():
        if not result.ok:
            raise KeyError(f"{pair.name} {label}: {result.reason}")
    job_results = {label: _job_result(r) for label, r in by_label.items()}
    cp = seamm_bsse.combine(specs(system, pair), job_results, n_atoms=pair.n_atoms)
    cluster = job_results["cluster"]
    d_energy = (cp.energy - cluster.energy) * HARTREE_KJ / EV_KJ
    d_gradient = np.asarray(cp.gradient) - np.asarray(cluster.gradient)
    d_forces = -d_gradient * HARTREE_KJ / BOHR / EV_KJ
    return d_energy, d_forces, bool(cp.gradient_fallback)
