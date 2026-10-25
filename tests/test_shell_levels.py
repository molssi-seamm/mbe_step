"""Levels of their own for the fragments holding an ion shell (e.g. Li+ shell pairs
at def2-TZVPPD with the solvent pairs at def2-QZVPPD), and more cores and memory for
the larger fragments."""

import json
import types

import numpy as np
import pytest

import mbe_step
from mbe_step import levels

from .helpers import configuration, level, parameters
from .test_phase3 import evaluate_with_bsse, run

SHIFT = 0.004  # eV added to every energy at the shifted level
SHIFTED = "high_shells"  # the Level (by name) whose energies are shifted

#: A water's H atoms relative to its O (Å), pointing away from +x
H_OFFSETS = np.array([[0.59, 0.76, 0.0], [0.59, -0.76, 0.0]])


def evaluate_with_shift(self, node, structures, stress=False):
    """The analytic levels; the SHIFTED level is the high level plus SHIFT."""
    if self.name == SHIFTED:
        out = evaluate_with_bsse(types.SimpleNamespace(name="high"), node, structures)
        for result in out.values():
            result.energy += SHIFT * 96.48533212331002  # kJ/mol
        return out
    return evaluate_with_bsse(self, node, structures, stress)


def li_cluster(db):
    """Li+ with four waters in its first shell (Li-O 1.9-2.05 Å) and three
    other waters."""
    symbols, xyz, bonds = ["Li"], [[0.0, 0.0, 0.0]], []
    waters = []
    for k, d in enumerate((1.9, 1.95, 2.0, 2.05)):
        theta = k * np.pi / 2
        u = np.array([np.cos(theta), np.sin(theta), 0.6 * (-1) ** k])
        waters.append(d * u / np.linalg.norm(u))
    waters += [np.array(p) for p in ((6.0, 0.0, 0.0), (0.0, 6.5, 0.0), (0.0, 0.0, 7.0))]
    for o in waters:
        n = len(symbols)
        away = o / np.linalg.norm(o)
        # H atoms pointing away from the ion: rotate +x onto the O's direction
        x = np.array([1.0, 0.0, 0.0])
        v = np.cross(x, away)
        s, c = np.linalg.norm(v), float(x @ away)
        if s < 1e-12:
            rot = np.eye(3)
        else:
            k = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]]) / s
            rot = np.eye(3) + s * k + (1 - c) * k @ k
        symbols += ["O", "H", "H"]
        xyz += [o, *(o + H_OFFSETS @ rot.T)]
        bonds += [(n, n + 1), (n, n + 2)]
    return configuration(db, symbols, np.array(xyz), bonds=bonds)


SETTINGS = {
    "pair cutoff": 100.0,
    "triple cutoff": 100.0,
    "triple rule": "compact",
    "distance criterion": "closest contact",
    "ion shells": "Li+ first shell",
    "shell truncation": "pairs",
}


@pytest.fixture
def db():
    from molsystem import SystemDB

    db = SystemDB(filename="file:mbe_shells?mode=memory&cache=shared")
    yield db
    db.close()


def _levels(**extra):
    out = {
        "high": level("high"),
        "molecular": level("molecular"),
        "cluster": level("cluster"),
    }
    for key, name in extra.items():
        out[key] = level(name)
    return out


def holds_shell(fragments, f):
    return bool(set(f.molecules) & fragments.shells)


def test_no_shell_levels_without_shells():
    node, P = parameters(**{"shell high level": "TZ", "shell low level": "TZLOW"})
    assert P["ion shells"] == "none"
    assert mbe_step.mbe._shell_level(P) is None
    assert mbe_step.mbe._shell_low_level(P) is None
    node, P = parameters(**SETTINGS, **{"shell high level": "same as the high level"})
    assert mbe_step.mbe._shell_level(P) is None
    node, P = parameters(**SETTINGS, **{"shell high level": "TZ"})
    assert mbe_step.mbe._shell_level(P) == "TZ"
    assert mbe_step.mbe._level_maps(P) == ({"shell": "high:shell"}, None)


def test_a_shell_level_moves_only_the_shell_increments(db, monkeypatch):
    monkeypatch.setattr(levels.Level, "evaluate", evaluate_with_shift)
    conf = li_cluster(db)
    node, P = parameters(**SETTINGS)
    frame = node._frame(conf, node._rules(P), P)
    assert frame["fragments"].shells
    assert "high:shell" not in frame["calculations"]
    _, plain = run(node, P, frame, _levels())

    node, P = parameters(**SETTINGS, **{"shell high level": "TZ"})
    frame = node._frame(conf, node._rules(P), P)
    fragments = frame["fragments"]
    calcs = frame["calculations"]
    # the fragments of two or more units holding a shell; the shell's own
    # monomer keeps the 1-body level (the labels' absolute energy scale)
    shell = [
        f for f in fragments.selected() if f.order >= 2 and holds_shell(fragments, f)
    ]
    assert set(calcs["high:shell"]) == set(fragments.closure(shell))
    results, labels = run(node, P, frame, _levels(**{"high:shell": "high_shells"}))
    assert results["high:shell"]
    n2 = len(shell)
    assert n2 == 3
    # -SHIFT for each shell pair's increment; the shell monomer's is unchanged
    assert labels.energy - plain.energy == pytest.approx(-n2 * SHIFT, abs=1e-9)
    assert labels.per_body[1]["energy"] == pytest.approx(
        plain.per_body[1]["energy"], abs=1e-12
    )
    assert labels.per_body[3]["energy"] == pytest.approx(
        plain.per_body[3]["energy"], abs=1e-12
    )
    names = {f.name for f in shell}
    for name, inc in frame["correction"].increments.items():
        expected = "high:shell" if name in names else "high"
        if inc.order < 3:
            assert inc.high_level == expected


def fake_resolve(text, context, current=None, periodic=False):
    return {"level": text, "step": "fake", "options": {}}


def test_a_shell_level_equal_to_the_triples_level_is_computed_once(db, monkeypatch):
    monkeypatch.setattr(levels, "resolve", fake_resolve)
    conf = li_cluster(db)
    settings = {**SETTINGS, "high level": "QZ", "molecular low level": "QZLOW"}
    node, P = parameters(
        **settings, **{"triples high level": "TZ", "shell high level": "TZ"}
    )
    frame = node._frame(conf, node._rules(P), P)
    monkeypatch.setattr(node, "variable_exists", lambda name: False)
    levels_ = node._levels(P, {}, [frame])
    assert levels_["high:shell"] is levels_["high:3"]
    assert levels_["high:3"].name == "high_triples"

    # Run once at the shared level: the union of the triples' and shells' ladders
    seen = {}

    def recording(self, node_, structures, stress=False):
        seen.setdefault(self.name, []).append(set(structures))
        return evaluate_with_bsse(self, node_, structures, stress)

    monkeypatch.setattr(levels.Level, "evaluate", recording)
    calcs = frame["calculations"]
    results = node._evaluate(levels_, [frame])
    prefix = frame["prefix"]
    union = {prefix + n for n in set(calcs["high:3"]) | set(calcs["high:shell"])}
    assert len(seen["high_triples"]) == 1
    assert seen["high_triples"][0] == union
    assert results["high:shell"] is results["high:3"]


def test_a_shell_level_equal_to_the_high_level_is_refused(db, monkeypatch):
    monkeypatch.setattr(levels, "resolve", fake_resolve)
    conf = li_cluster(db)
    settings = {**SETTINGS, "high level": "QZ", "molecular low level": "QZLOW"}
    node, P = parameters(**settings, **{"shell low level": "QZLOW2"})
    frame = node._frame(conf, node._rules(P), P)
    monkeypatch.setattr(mbe_step.Mbe, "variable_exists", lambda self, name: False)
    levels_ = node._levels(P, {}, [frame])
    assert levels_["molecular:shell"].name == "molecular_shells"
    # the same text as the high level is "same as the high level"
    node, P = parameters(**settings, **{"shell high level": "QZ"})
    assert mbe_step.mbe._shell_level(P) is None
    # a different spelling of the same level is refused
    monkeypatch.setattr(
        levels, "resolve", lambda text, *a, **k: fake_resolve(text.upper(), None)
    )
    node, P = parameters(**settings, **{"shell high level": "qz"})
    with pytest.raises(ValueError, match="computed twice"):
        node._levels(P, {}, [frame])


# ------------------------------------------------------- fragment resources
def test_parse_resources():
    assert levels.parse_resources("") == []
    assert levels.parse_resources("none") == []
    assert levels.parse_resources("60 16 8000 1; 40 4 6000") == [
        (40, 4, 6000.0, None),
        (60, 16, 8000.0, 1),
    ]
    for bad in (
        "60 16",
        "a 4 6000",
        "40 0 6000",
        "40 4 6000; 40 8 8000",
        "40 4 6000 0",
        "40 4 6000 1 2",
    ):
        with pytest.raises(ValueError):
            levels.parse_resources(bad)


def _structure(n):
    from seamm_exec import Geometry

    return Geometry([1] * n, np.zeros((n, 3)))


def test_split_by_size():
    lv = levels.Level(
        "high",
        {"level": "x"},
        ranks=2,
        memory=3800,
        bundle=8,
        walltime=0,
        archive=False,
        tiers=levels.parse_resources("40 4 6000; 60 16 8000 1"),
    )
    structures = {f"s{n}": _structure(n) for n in (10, 39, 40, 59, 60, 114)}
    structures["opt"] = (_structure(70), {"grid": {}})
    parts = {name: (r, m, b, sorted(s)) for name, r, m, b, s in lv.split(structures)}
    # a tier's own bundle size, else the level's (8)
    assert parts == {
        "high": (2, 3800, 8, ["s10", "s39"]),
        "high_40": (4, 6000.0, 8, ["s40", "s59"]),
        "high_60": (16, 8000.0, 1, ["opt", "s114", "s60"]),
    }
    # a tier with nothing in it is left out
    small = {"a": _structure(5)}
    assert [p[0] for p in lv.split(small)] == ["high"]


def test_tiers_run_as_their_own_evaluators(monkeypatch):
    calls = []

    def fake(self, node, structures, directory, ranks, memory, stress, bundle=None):
        calls.append((directory, ranks, memory, bundle, sorted(structures)))
        return {key: key for key in structures}

    monkeypatch.setattr(levels.Level, "_evaluate", fake)
    lv = levels.Level(
        "molecular",
        {"level": "x"},
        ranks=2,
        memory=3800,
        bundle=8,
        walltime=0,
        archive=False,
        tiers=[(60, 16, 8000.0, 1)],
    )
    structures = {"a": _structure(20), "b": _structure(90)}
    results = lv.evaluate(None, structures)
    assert results == {"a": "a", "b": "b"}
    assert sorted(calls) == [
        ("molecular", 2, 3800, 8, ["a"]),
        ("molecular_60", 16, 8000.0, 1, ["b"]),
    ]
    # without tiers, one Evaluator in the level's own directory
    calls.clear()
    lv.tiers = []
    lv.evaluate(None, structures)
    assert calls == [("molecular", 2, 3800.0, 8, ["a", "b"])]


def two_shells(db, separation=8.0):
    """Two Li+ with four waters each in their first shells, and a free water."""
    symbols, xyz, bonds = [], [], []
    centers = [np.zeros(3), np.array([separation, 0.0, 0.0])]
    for c in centers:
        symbols.append("Li")
        xyz.append(c)
    waters = []
    for c in centers:
        for k, d in enumerate((1.9, 1.95, 2.0, 2.05)):
            theta = k * np.pi / 2
            u = np.array([0.6 * (-1) ** k, np.cos(theta), np.sin(theta)])
            waters.append((c + d * u / np.linalg.norm(u), c))
    waters.append(
        (np.array([separation / 2, 6.0, 0.0]), np.array([separation / 2, 7.0, 0.0]))
    )
    for o, center in waters:
        n = len(symbols)
        away = (o - center) / np.linalg.norm(o - center)
        x = np.array([1.0, 0.0, 0.0])
        v = np.cross(x, away)
        s, c = np.linalg.norm(v), float(x @ away)
        if s < 1e-12:
            rot = np.eye(3) if c > 0 else np.diag([-1.0, -1.0, 1.0])
        else:
            k = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]]) / s
            rot = np.eye(3) + s * k + (1 - c) * k @ k
        symbols += ["O", "H", "H"]
        xyz += [o, *(o + H_OFFSETS @ rot.T)]
        bonds += [(n, n + 1), (n, n + 2)]
    return configuration(db, symbols, np.array(xyz), bonds=bonds)


def test_the_shell_pair_cap(db, monkeypatch, tmp_path):
    """A pair of two shells at or beyond the cap is not computed at any level,
    its increment leaves the sum, and increments_c<id>.json lists it as capped
    (Li-salt production, 2026-10-25)."""
    monkeypatch.setattr(levels.Level, "evaluate", evaluate_with_shift)
    monkeypatch.setattr(mbe_step.Mbe, "directory", str(tmp_path), raising=False)
    conf = two_shells(db)
    shell_levels = {"shell high level": "TZ", "shell low level": "TZLOW"}
    extra = {"high:shell": "high_shells", "molecular:shell": "molecular_shells"}

    node, P = parameters(**SETTINGS, **shell_levels)
    frame = node._frame(conf, node._rules(P), P)
    fragments = frame["fragments"]
    assert len(fragments.shells) == 2 and fragments.capped == []
    (pair,) = [f for f in fragments.by_order(2) if set(f.molecules) <= fragments.shells]
    _, plain = run(node, P, frame, _levels(**extra))
    increment = frame["correction"].increments[pair.name]

    node, P = parameters(
        **SETTINGS, **shell_levels, **{"shell pair cap": pair.max_distance - 0.01}
    )
    node._id = (1,)  # the header needs a step id, as in a flowchart
    text = " ".join(node.description_text(P).split())
    assert "are not computed (capped)" in text
    frame = node._frame(conf, node._rules(P), P)
    fragments = frame["fragments"]
    assert [c["name"] for c in fragments.capped] == [pair.name]
    assert pair.name not in fragments
    for level_, names in frame["calculations"].items():
        assert pair.name not in names, level_
    _, capped = run(node, P, frame, _levels(**extra))
    assert capped.energy == pytest.approx(plain.energy - increment.energy, abs=1e-9)
    assert np.allclose(
        capped.forces, plain.forces - _on_atoms(increment, pair, len(plain.forces))
    )

    node._write_increments(frame, P)
    data = json.loads((tmp_path / f"increments_c{conf.id}.json").read_text())
    assert data["shell pair cap (Å)"] == pytest.approx(pair.max_distance - 0.01)
    records = {r["name"]: r for r in data["increments"]}
    record = records[pair.name]
    assert record["capped"] == f"{pair.max_distance - 0.01:g}A"
    assert record["computed"] is False
    assert "energy" not in record and "forces" not in record
    assert record["distances"]["0,1"] == pytest.approx(pair.max_distance)
    assert len(record["members"]) == 2
    assert all("capped" not in r for n, r in records.items() if n != pair.name)

    # A cap beyond the pair keeps it
    node, P = parameters(
        **SETTINGS, **shell_levels, **{"shell pair cap": pair.max_distance + 0.01}
    )
    frame = node._frame(conf, node._rules(P), P)
    assert frame["fragments"].capped == [] and pair.name in frame["fragments"]


def _on_atoms(increment, fragment, n_atoms):
    forces = np.zeros((n_atoms, 3))
    np.add.at(forces, fragment.atoms, increment.forces)
    return forces


def test_the_shell_pair_cap_needs_shells():
    node, P = parameters(**{"shell pair cap": 2.5})
    assert node._rules(P).shell_pair_cap is None
    node, P = parameters(**SETTINGS, **{"shell pair cap": 2.5})
    assert node._rules(P).shell_pair_cap == pytest.approx(2.5)
    node, P = parameters(**SETTINGS)
    assert node._rules(P).shell_pair_cap is None
