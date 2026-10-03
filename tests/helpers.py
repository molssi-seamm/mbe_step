"""Helpers for the mbe_step tests: molsystem configurations, a stand-in for the
Evaluator results, and driving the step's frame/evaluate/assemble path."""

import itertools
from types import SimpleNamespace

import numpy as np

import mbe_step
from mbe_step import levels
from seamm_mbe.units import EV_TO_KJ_PER_MOL


def configuration(db, symbols, xyz, cell=None, bonds=None):
    """A molsystem configuration in ``db``."""
    system = db.create_system()
    periodicity = 0 if cell is None else 3
    configuration = system.create_configuration(periodicity=periodicity)
    if cell is not None:
        a, b, c = np.linalg.norm(np.asarray(cell), axis=1)
        configuration.cell.parameters = [a, b, c, 90.0, 90.0, 90.0]
    xyz = np.asarray(xyz, dtype=float)
    configuration.atoms.append(x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2], symbol=symbols)
    if cell is not None:
        configuration.coordinate_system = "Cartesian"
    if bonds:
        ids = configuration.atoms.ids
        configuration.bonds.append(
            i=[ids[i] for i, _ in bonds], j=[ids[j] for _, j in bonds]
        )
    return configuration


def parameters(**values):
    """The step's dereferenced parameters, as run() sees them."""
    node = mbe_step.Mbe()
    for key, value in values.items():
        node.parameters[key].value = value
    return node, node.parameters.current_values_to_dict(context={})


def level(name, convention=None):
    """A Level with a fake model chemistry."""
    mc = {"level": f"fake {name}", "step": "fake", "options": {}}
    result = levels.Level(
        name, mc, ranks=1, memory=100, bundle=10, walltime=0, archive=False
    )
    if convention is not None:
        result.convention = convention
    return result


def result(key, energy_ev, forces_ev, stress_gpa=None):
    """An EvaluatorResult stand-in in the contract's units (kJ/mol, gradients
    in kJ/mol/Å, stress in GPa)."""
    return SimpleNamespace(
        key=key,
        ok=True,
        energy=energy_ev * EV_TO_KJ_PER_MOL,
        gradients=-np.asarray(forces_ev) * EV_TO_KJ_PER_MOL,
        stress=stress_gpa,
        reason=None,
        restored=False,
        path="batch",
    )


def model_potential(z, xyz, charges, c3):
    """A potential with only 1-, 2- and 3-body terms between molecules: a
    harmonic intramolecular-like term (1-body), screened charges (2-body) and an
    Axilrod-Teller-like term between atoms (3-body). Energy in eV, forces in
    eV/Å, analytic for the 2-body part, central differences for the rest."""

    def energy(x):
        e = 0.0
        for i, j in itertools.combinations(range(len(x)), 2):
            r = np.linalg.norm(x[i] - x[j])
            e += charges[i] * charges[j] * np.exp(-0.3 * r) / r
        if c3:
            for i, j, k in itertools.combinations(range(len(x)), 3):
                rij = np.linalg.norm(x[i] - x[j])
                rjk = np.linalg.norm(x[j] - x[k])
                rik = np.linalg.norm(x[i] - x[k])
                e += c3 / (rij * rjk * rik) ** 2
        return e

    x = np.asarray(xyz, dtype=float)
    forces = np.zeros_like(x)
    h = 1e-5
    for a in range(len(x)):
        for d in range(3):
            xp = x.copy()
            xm = x.copy()
            xp[a, d] += h
            xm[a, d] -= h
            forces[a, d] = -(energy(xp) - energy(xm)) / (2 * h)
    return energy(x), forces
