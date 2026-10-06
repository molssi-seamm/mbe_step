"""A high level of their own for the triples (e.g. QZ pairs, TZ triples), and the
per-increment file."""

import json
import types

import numpy as np
import pytest

import mbe_step
from mbe_step import levels

from .helpers import level, parameters
from .test_phase3 import cluster, evaluate_with_bsse, run

SHIFT = 0.004  # eV added to every energy at the triples' level


def evaluate_with_triples_shift(self, node, structures, stress=False):
    """The analytic levels; the triples' level is the high level plus SHIFT."""
    if self.name == "high_triples":
        as_high = types.SimpleNamespace(name="high")
        out = evaluate_with_bsse(as_high, node, structures, stress)
        for result in out.values():
            result.energy += SHIFT * 96.48533212331002  # kJ/mol
        return out
    return evaluate_with_bsse(self, node, structures, stress)


SETTINGS = {"pair cutoff": 100.0, "triple cutoff": 100.0, "triple rule": "compact"}


@pytest.fixture
def db():
    from molsystem import SystemDB

    db = SystemDB(filename="file:mbe_triples?mode=memory&cache=shared")
    yield db
    db.close()


def _levels(triples):
    out = {
        "high": level("high"),
        "molecular": level("molecular"),
        "cluster": level("cluster"),
    }
    if triples:
        out["high:3"] = level("high_triples")
    return out


def test_same_triples_level_is_the_single_high_level(db, monkeypatch):
    monkeypatch.setattr(levels.Level, "evaluate", evaluate_with_triples_shift)
    conf = cluster(db)
    node, P = parameters(**SETTINGS)
    frame = node._frame(conf, node._rules(P), P)
    _, plain = run(node, P, frame, _levels(False))
    assert "high:3" not in frame["calculations"]

    # Naming the high level itself is the same as "same as the high level"
    node, P = parameters(**SETTINGS, **{"triples high level": P["high level"]})
    assert mbe_step.mbe._triples_level(P) is None
    frame = node._frame(conf, node._rules(P), P)
    assert "high:3" not in frame["calculations"]
    _, same = run(node, P, frame, _levels(False))
    assert same.energy == plain.energy
    assert np.array_equal(same.forces, plain.forces)


def test_a_triples_level_moves_only_the_triples(db, monkeypatch):
    monkeypatch.setattr(levels.Level, "evaluate", evaluate_with_triples_shift)
    conf = cluster(db)
    node, P = parameters(**SETTINGS)
    frame = node._frame(conf, node._rules(P), P)
    _, plain = run(node, P, frame, _levels(False))

    node, P = parameters(**SETTINGS, **{"triples high level": "TRIPLES"})
    frame = node._frame(conf, node._rules(P), P)
    calcs = frame["calculations"]
    assert "high:3" in calcs
    # The triples with all their sub-pairs and monomers at the triples' level
    triples = [f for f in frame["fragments"].selected() if f.order == 3]
    assert set(calcs["high:3"]) == set(frame["fragments"].closure(triples))
    results, labels = run(node, P, frame, _levels(True))
    assert results["high:3"]
    n_triples = len(triples)
    assert n_triples == 4
    # Each triple's increment moves by the shift (a constant cancels below)
    assert labels.energy - plain.energy == pytest.approx(n_triples * SHIFT, abs=1e-9)
    assert labels.per_body[2]["energy"] == pytest.approx(
        plain.per_body[2]["energy"], abs=1e-12
    )
    assert np.allclose(labels.forces, plain.forces, atol=1e-10)


def test_the_increments_file(db, monkeypatch, tmp_path):
    monkeypatch.setattr(levels.Level, "evaluate", evaluate_with_triples_shift)
    monkeypatch.setattr(mbe_step.Mbe, "directory", str(tmp_path), raising=False)
    conf = cluster(db)
    node, P = parameters(
        **{**SETTINGS, "triple cutoff": 3.2, "triple rule": "connected"},
        **{"triples high level": "TRIPLES"},
    )
    frame = node._frame(conf, node._rules(P), P)
    run(node, P, frame, _levels(True))
    node._write_increments(frame, P)
    data = json.loads((tmp_path / f"increments_c{conf.id}.json").read_text())
    records = data["increments"]
    assert {r["order"] for r in records} == {1, 2, 3}
    assert all(r["high level"] == "high:3" for r in records if r["order"] == 3)
    assert all(r["high level"] == "high" for r in records if r["order"] < 3)
    for r in records:
        if r["order"] == 3:
            n = sum(1 for d in r["distances"].values() if d < 3.2)
            assert r["topology"] == ("closed" if n == 3 else "chain")
            assert len(r["forces"]) == frame["fragments"][r["name"]].n_atoms
    assert data["units"] == {"energy": "eV", "forces": "eV/Å", "distances": "Å"}
    total = sum(r["energy"] for r in records)
    assert total == pytest.approx(frame["correction"].energy, abs=1e-10)


def test_no_triples_level_without_triples():
    """Order 3 with triple rule 'none' (e.g. a hand-edited flowchart) has no
    triples, so a triples' level set anyway is ignored."""
    node, P = parameters(
        **{"triple rule": "none", "maximum order": "3"},
        **{"triples high level": "TRIPLES"},
    )
    assert mbe_step.mbe._triples_level(P) is None
    node._id = (1,)  # the header needs a step id, as in a flowchart
    assert "TRIPLES" not in node.description_text(P)
    node, P = parameters(
        **{"triple rule": "connected"}, **{"triples high level": "TRIPLES"}
    )
    assert mbe_step.mbe._triples_level(P) == "TRIPLES"
    node._id = (1,)
    assert "TRIPLES" in node.description_text(P)


def test_a_triples_level_that_is_the_high_level_is_refused(db, monkeypatch):
    """Two spellings of one level would compute every monomer and pair twice."""

    def resolve(text, context, current=None, periodic=False):
        return {"level": text.upper(), "step": "fake", "options": {}}

    monkeypatch.setattr(levels, "resolve", resolve)
    monkeypatch.setattr(mbe_step.Mbe, "variable_exists", lambda self, name: False)
    conf = cluster(db)
    node, P = parameters(
        **SETTINGS,
        **{"high level": "orca:x", "triples high level": "ORCA:X"},
        **{"molecular low level": "LOW"},
    )
    frame = node._frame(conf, node._rules(P), P)
    with pytest.raises(ValueError, match="is the high level"):
        node._levels(P, {}, [frame])


def test_topology_with_a_cutoff_table(db, monkeypatch, tmp_path):
    """Chains and closed triples are classified with the cutoffs the enumerator
    used, by type pair, so table mode gets the topology too."""
    monkeypatch.setattr(levels.Level, "evaluate", evaluate_with_triples_shift)
    monkeypatch.setattr(mbe_step.Mbe, "directory", str(tmp_path), raising=False)
    conf = cluster(db)
    node, P = parameters(
        **{
            "cutoffs": "table by type pair",
            "pair cutoff table": "* * 100.0",
            "triple cutoff table": "* * 3.2",
            "triple rule": "connected",
        }
    )
    frame = node._frame(conf, node._rules(P), P)
    run(node, P, frame, _levels(False))
    node._write_increments(frame, P)
    data = json.loads((tmp_path / f"increments_c{conf.id}.json").read_text())
    triples = [r for r in data["increments"] if r["order"] == 3]
    assert triples
    for r in triples:
        n = sum(1 for d in r["distances"].values() if d < 3.2)
        assert r["topology"] == ("closed" if n == 3 else "chain")


# ------------------------------------------------- the triples' own low level
def evaluate_with_low_triples_shift(self, node, structures, stress=False):
    """The analytic levels; the triples' low level is the molecular low level
    plus SHIFT."""
    if self.name == "molecular_triples":
        as_low = types.SimpleNamespace(name="molecular")
        out = evaluate_with_bsse(as_low, node, structures, stress)
        for result in out.values():
            result.energy += SHIFT * 96.48533212331002  # kJ/mol
        return out
    return evaluate_with_bsse(self, node, structures, stress)


def _low_levels(triples):
    out = _levels(False)
    if triples:
        out["molecular:3"] = level("molecular_triples")
    return out


def test_a_triples_low_level_moves_only_the_triples(db, monkeypatch):
    monkeypatch.setattr(levels.Level, "evaluate", evaluate_with_low_triples_shift)
    conf = cluster(db)
    node, P = parameters(**SETTINGS)
    frame = node._frame(conf, node._rules(P), P)
    _, plain = run(node, P, frame, _low_levels(False))

    node, P = parameters(**SETTINGS, **{"triples low level": "TRIPLES"})
    assert mbe_step.mbe._triples_low_level(P) == "TRIPLES"
    frame = node._frame(conf, node._rules(P), P)
    calcs = frame["calculations"]
    assert "molecular:3" in calcs
    triples = [f for f in frame["fragments"].selected() if f.order == 3]
    assert set(calcs["molecular:3"]) == set(frame["fragments"].closure(triples))
    assert all(frame["fragments"][n].order <= 2 for n in calcs["molecular"])
    _, labels = run(node, P, frame, _low_levels(True))
    # The low level is subtracted: each triple's increment moves by -SHIFT
    assert labels.energy - plain.energy == pytest.approx(
        -len(triples) * SHIFT, abs=1e-9
    )
    assert labels.per_body[2]["energy"] == pytest.approx(
        plain.per_body[2]["energy"], abs=1e-12
    )
    assert np.allclose(labels.forces, plain.forces, atol=1e-10)


def test_same_triples_low_level_is_the_single_low_level(db, monkeypatch):
    monkeypatch.setattr(levels.Level, "evaluate", evaluate_with_low_triples_shift)
    conf = cluster(db)
    node, P = parameters(**SETTINGS)
    frame = node._frame(conf, node._rules(P), P)
    _, plain = run(node, P, frame, _low_levels(False))
    node, P = parameters(
        **SETTINGS, **{"triples low level": "same as the molecular low level"}
    )
    assert mbe_step.mbe._triples_low_level(P) is None
    frame = node._frame(conf, node._rules(P), P)
    assert "molecular:3" not in frame["calculations"]
    _, same = run(node, P, frame, _low_levels(False))
    assert same.energy == plain.energy
    assert np.array_equal(same.forces, plain.forces)


def test_the_increments_file_records_the_low_level(db, monkeypatch, tmp_path):
    monkeypatch.setattr(levels.Level, "evaluate", evaluate_with_low_triples_shift)
    monkeypatch.setattr(mbe_step.Mbe, "directory", str(tmp_path), raising=False)
    conf = cluster(db)
    node, P = parameters(**SETTINGS, **{"triples low level": "TRIPLES"})
    frame = node._frame(conf, node._rules(P), P)
    run(node, P, frame, _low_levels(True))
    node._write_increments(frame, P)
    data = json.loads((tmp_path / f"increments_c{conf.id}.json").read_text())
    records = data["increments"]
    assert all(r["low level"] == "molecular:3" for r in records if r["order"] == 3)
    assert all(r["low level"] == "molecular" for r in records if r["order"] < 3)
    node._id = (1,)
    assert "TRIPLES" in node.description_text(P)


def test_a_triples_low_level_that_is_the_low_level_is_refused(db, monkeypatch):
    def resolve(text, context, current=None, periodic=False):
        return {"level": text.upper(), "step": "fake", "options": {}}

    monkeypatch.setattr(levels, "resolve", resolve)
    monkeypatch.setattr(mbe_step.Mbe, "variable_exists", lambda self, name: False)
    conf = cluster(db)
    node, P = parameters(
        **SETTINGS,
        **{"molecular low level": "orca:low", "triples low level": "ORCA:LOW"},
    )
    frame = node._frame(conf, node._rules(P), P)
    with pytest.raises(ValueError, match="is the molecular low level"):
        node._levels(P, {}, [frame])
