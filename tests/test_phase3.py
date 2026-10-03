"""Phase 3: pairwise counterpoise, the grid options of the periodic levels, and
the checks on the levels."""

import numpy as np
import pytest
from seamm_exec import structure_data

from mbe_step import levels

from .helpers import configuration, level, model_potential, parameters, result
from .test_cluster_identity import WATERS, charges_for

BSSE = {"high": 0.003, "molecular": 0.001}  # eV per monomer in a pair's basis


def evaluate_with_bsse(self, node, structures, stress=False):
    """Analytic levels; a monomer in a pair's basis is lower by the level's
    BSSE, its gradient on the ghost atoms zero."""
    out = {}
    c3 = 0.5 if self.name == "high" else 0.0
    scale = 1.0 if self.name == "high" else 0.7
    for key, structure in structures.items():
        options = {}
        if isinstance(structure, tuple):
            structure, options = structure
        data = structure_data(structure)
        z = np.asarray(data["atomic_numbers"])
        xyz = np.asarray(data["coordinates"], dtype=float)
        ghosts = set(options.get("ghost_atoms") or ())
        real = [i for i in range(len(z)) if i not in ghosts]
        e, f = model_potential(z[real], xyz[real], charges_for(z[real], scale), c3)
        forces = np.zeros_like(xyz)
        forces[real] = f
        if ghosts:
            e -= BSSE[self.name]
        out[key] = result(key, e, forces)
    return out


@pytest.fixture
def db():
    from molsystem import SystemDB

    db = SystemDB(filename="file:mbe_phase3?mode=memory&cache=shared")
    yield db
    db.close()


def cluster(db):
    return configuration(
        db,
        ["O", "H", "H"] * 4,
        WATERS,
        bonds=[(3 * m, 3 * m + k) for m in range(4) for k in (1, 2)],
    )


def run(node, P, frame, levels_):
    results = node._evaluate(levels_, [frame])
    return results, node._assemble(frame, levels_, results, None)


def test_pairwise_counterpoise(db, monkeypatch):
    monkeypatch.setattr(levels.Level, "evaluate", evaluate_with_bsse)
    settings = {"pair cutoff": 100.0, "triple cutoff": 100.0, "triple rule": "compact"}
    conf = cluster(db)
    levels_ = {
        "high": level("high"),
        "molecular": level("molecular"),
        "cluster": level("cluster"),
    }

    node, P = parameters(**settings)
    frame = node._frame(conf, node._rules(P), P)
    _, plain = run(node, P, frame, levels_)

    node, P = parameters(**settings, counterpoise="pairwise")
    node._counterpoise = True
    frame = node._frame(conf, node._rules(P), P)
    results, labels = run(node, P, frame, levels_)

    pairs = frame["fragments"].by_order(2, in_sum=True)
    ghost_keys = [k for k in results["high"] if "-cp-" in k]
    assert len(ghost_keys) == 2 * len(pairs) == 12
    assert len([k for k in results["molecular"] if "-cp-" in k]) == 12
    # E^CP = E - sum(E_in - E_alone) = E + 2 b at each level; [high - low]
    expected = 2 * (BSSE["high"] - BSSE["molecular"]) * len(pairs)
    assert labels.energy - plain.energy == pytest.approx(expected, abs=1e-12)
    assert labels.per_body[2]["energy"] - plain.per_body[2]["energy"] == (
        pytest.approx(expected, abs=1e-12)
    )
    assert labels.per_body[3]["energy"] == pytest.approx(
        plain.per_body[3]["energy"], abs=1e-12
    )
    # The ghost gradients were zero: no change of the forces, no fallback
    assert np.allclose(labels.forces, plain.forces, atol=1e-10)
    assert frame["counterpoise fallbacks"] == 0


def test_periodic_levels_get_the_grid(db, monkeypatch):
    seen = {}

    def evaluate(self, node, structures, stress=False):
        out = {}
        for key, structure in structures.items():
            options = {}
            if isinstance(structure, tuple):
                structure, options = structure
            seen[(self.name, key)] = options
            n = len(structure_data(structure)["atomic_numbers"])
            stress_gpa = np.zeros((3, 3)) if self.name == "cell" else None
            out[key] = result(key, -1.0, np.zeros((n, 3)), stress_gpa)
        return out

    monkeypatch.setattr(levels.Level, "evaluate", evaluate)
    conf = configuration(
        db,
        ["O", "H", "H"] * 2,
        np.array(WATERS[:6]) + 3.0,
        cell=np.eye(3) * 12.0,
        bonds=[(0, 1), (0, 2), (3, 4), (3, 5)],
    )
    node, P = parameters(
        **{"maximum order": "2", "periodic low level": "fake", "grid spacing": 0.08}
    )
    node._grid = {"max_spacing": 0.08, "padding": 7.5}
    frame = node._frame(conf, node._rules(P), P)
    levels_ = {
        "high": level("high"),
        "molecular": level("molecular"),
        "periodic": level("periodic"),
        "cell": level("cell", convention="stress"),
    }
    node._evaluate(levels_, [frame])
    periodic = [o for (name, _), o in seen.items() if name == "periodic"]
    assert periodic and all(
        o["grid"]["max_spacing"] == 0.08
        and o["grid"]["padding"] == 7.5
        and np.allclose(o["grid"]["reference_cell"], np.eye(3) * 12.0)
        for o in periodic
    )
    (cell,) = [o for (name, _), o in seen.items() if name == "cell"]
    assert cell == {"grid": {"max_spacing": 0.08}}
    assert all(o == {} for (name, _), o in seen.items() if name == "high")


def fake_resolve(text, context, current=None, periodic=False):
    options = {}
    if text.startswith("MDI"):
        options["mdi_capable"] = True
    return {"level": text, "step": "fake", "options": options}


def test_level_checks(db, monkeypatch):
    import mbe_step

    monkeypatch.setattr(levels, "resolve", fake_resolve)
    monkeypatch.setattr(mbe_step.Mbe, "variable_exists", lambda self, name: False)
    monkeypatch.setattr(levels, "stress_convention", lambda mc: "stress")
    cell = configuration(
        db,
        ["O", "H", "H"] * 2,
        np.array(WATERS[:6]) + 3.0,
        cell=np.eye(3) * 12.0,
        bonds=[(0, 1), (0, 2), (3, 4), (3, 5)],
    )
    base = {
        "high level": "HIGH",
        "molecular low level": "LOW",
        "periodic low level": "VASP:A",
        "maximum order": "2",
    }
    node, P = parameters(**base, **{"cell low level": "VASP:B"})
    frame = node._frame(cell, node._rules(P), P)
    with pytest.raises(ValueError, match="must be the periodic fragments'"):
        node._levels(P, {}, [frame])
    node, P = parameters(**base)
    node._levels(P, {}, [frame])  # automatic: the same level

    node, P = parameters(**base, counterpoise="pairwise")
    with pytest.raises(ValueError, match="not available for periodic cells"):
        node._levels(P, {}, [frame])

    water = cluster(db)
    node, P = parameters(
        **{"high level": "MDI-HIGH", "molecular low level": "LOW"},
        counterpoise="pairwise",
    )
    frame = node._frame(water, node._rules(P), P)
    with pytest.raises(ValueError, match="needs ghost atoms"):
        node._levels(P, {}, [frame])
