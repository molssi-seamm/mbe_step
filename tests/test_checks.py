"""Checks made before any calculation, the whole system's charge, a cell
without a stress, and a smoke test of the real Evaluator."""

import logging
import types

import numpy as np
import pytest

import seamm_exec
import seamm_mbe
from mbe_step import levels

from .helpers import configuration, level, parameters, result

WATERS = [
    [0.000, 0.000, 0.000],
    [0.957, 0.000, 0.000],
    [-0.240, 0.927, 0.000],
    [2.900, 0.200, 0.100],
    [3.400, 1.000, 0.200],
    [3.500, -0.500, 0.300],
]


@pytest.fixture
def db():
    from molsystem import SystemDB

    db = SystemDB(filename="file:mbe_checks?mode=memory&cache=shared")
    yield db
    db.close()


def water_dimer(db):
    return configuration(
        db, ["O", "H", "H"] * 2, WATERS, bonds=[(0, 1), (0, 2), (3, 4), (3, 5)]
    )


def test_missing_offsets_are_refused_before_running(db):
    conf = configuration(
        db,
        ["F", "H", "F", "H"],
        [[0, 0, 0], [0.92, 0, 0], [3, 0, 0], [3.92, 0, 0]],
        bonds=[(0, 1), (2, 3)],
    )
    node, P = parameters(**{"energy offsets": "water 2074.69325"})
    with pytest.raises(ValueError, match=r"type\(s\) \['FH'\] with no energy offset"):
        node._frame(conf, node._rules(P), P)
    node, P = parameters()  # default: none
    node._frame(conf, node._rules(P), P)


def test_whole_system_uses_the_molecules_charges(db):
    """Li+ and two waters from plain coordinates: the configuration says
    nothing about the charge, but the whole system runs at +1, singlet, the
    same charges the fragments use."""
    xyz = [[-2.0, 0.0, 0.0]] + WATERS
    conf = configuration(
        db, ["Li"] + ["O", "H", "H"] * 2, xyz, bonds=[(1, 2), (1, 3), (4, 5), (4, 6)]
    )
    node, P = parameters()
    frame = node._frame(conf, node._rules(P), P)
    whole = levels.whole(frame["system"])
    data = seamm_exec.structure_data(whole)
    assert data["charge"] == 1 and data["multiplicity"] == 1
    assert data["periodicity"] == 0
    monomer = levels.geometry(frame["system"], frame["fragments"]["m00"])
    assert seamm_exec.structure_data(monomer)["charge"] == 1

    conf.charge = 2
    with pytest.raises(ValueError, match="charge \\+2, but its molecules'"):
        node._frame(conf, node._rules(P), P)


def test_open_shell_molecules_in_fragments_are_refused():
    system = seamm_mbe.System(
        ["O", "H", "O", "H"],
        [[0, 0, 0], [0.97, 0, 0], [3, 0, 0], [3.97, 0, 0]],
        bonds=[(0, 1), (2, 3)],
        multiplicities={"HO": 2},
    )
    node, P = parameters()
    conf = types.SimpleNamespace(name="radicals", charge=0)
    with pytest.raises(ValueError, match="open-shell molecules \\(HO\\)"):
        node._check_frame(conf, system, node._rules(P), P)
    node, P = parameters(**{"maximum order": "1"})
    node._check_frame(conf, system, node._rules(P), P)


def test_cell_without_stress_fails_early(db, monkeypatch):
    def evaluate(self, node, structures, stress=False):
        out = {}
        for key, structure in structures.items():
            n = len(seamm_exec.structure_data(structure)["atomic_numbers"])
            out[key] = result(key, -1.0, np.zeros((n, 3)))  # no stress
        return out

    monkeypatch.setattr(levels.Level, "evaluate", evaluate)
    box = 12.0
    conf = configuration(
        db,
        ["O", "H", "H"] * 2,
        np.array(WATERS) + 3.0,
        cell=np.eye(3) * box,
        bonds=[(0, 1), (0, 2), (3, 4), (3, 5)],
    )
    node, P = parameters(**{"maximum order": "2", "pair cutoff": 5.0})
    frame = node._frame(conf, node._rules(P), P)
    levels_ = {
        "high": level("high"),
        "molecular": level("molecular"),
        "cell": level("cell", convention="pressure"),
    }
    results = node._evaluate(levels_, [frame])
    (cell,) = results["cell"].values()
    assert not cell.ok and "returned no stress" in cell.reason
    with pytest.raises(seamm_mbe.MissingFragmentsError):
        node._assemble(frame, levels_, results, None)


# ---- the real Evaluator ------------------------------------------------------
class FakeProvider:
    """A program whose 'calculation' is a shell task writing a number."""

    seen = []

    @classmethod
    def get_task(
        cls, configuration, model_chemistry, *, key, properties, options, resources=None
    ):
        cls.seen.append(resources)
        return seamm_exec.Task(
            key=key,
            program="fake",
            cmd=["echo", "-1.5", ">", "energy.txt"],
            shell=True,
            return_files=["energy.txt"],
            config={"installation": "local"},
        )

    @classmethod
    def analyze_task(
        cls, result, model_chemistry, configuration, *, properties, options
    ):
        n = len(seamm_exec.structure_data(configuration)["atomic_numbers"])
        return {
            "energy": float(result.files["energy.txt"]),
            "gradients": np.zeros((n, 3)),
        }


def test_real_evaluator_with_resources(tmp_path):
    """Level.evaluate through the real seamm_exec Evaluator, batch path. This
    fails (TypeError) with a seamm_exec that predates Evaluator(resources=):
    mbe_step must pin a seamm-exec that has it."""
    node = types.SimpleNamespace(
        directory=str(tmp_path / "step"),
        global_options={"root": str(tmp_path)},
        logger=logging.getLogger("test"),
        variable_exists=lambda name: False,
        flowchart=types.SimpleNamespace(
            executor=seamm_exec.Local(),
            plugin_manager=types.SimpleNamespace(get=lambda name: FakeProvider),
            root_directory=str(tmp_path),
        ),
    )
    mc = {
        "level": "FAKE:X@Y",
        "method": "Y",
        "basis": None,
        "step": "fake",
        "options": {"mdi_capable": False},
    }
    lvl = levels.Level(
        "high", mc, ranks=4, memory=1500, bundle=10, walltime=0, archive=True
    )
    geometry = seamm_exec.Geometry([8, 1, 1], WATERS[:3])
    results = lvl.evaluate(node, {"c1-m00": geometry})
    assert results["c1-m00"].ok, results["c1-m00"].reason
    assert results["c1-m00"].energy == pytest.approx(-1.5)
    (resources,) = FakeProvider.seen
    assert resources.ntasks == 4 and resources.mem_per_cpu == 1_500_000_000
