"""Explicit auxiliary basis sets for single fragments: a shell pair whose AutoAux
calculation fails (MBE cell f04, 2026-10-23) is computed, with its monomers, at a
copy of its level with explicit sets, and its increment built from those alone."""

import types

import pytest

import mbe_step
from mbe_step import levels

from .helpers import parameters
from .test_phase3 import evaluate_with_bsse, run
from .test_shell_levels import SETTINGS, _levels, holds_shell, li_cluster


@pytest.fixture
def db():
    from molsystem import SystemDB

    db = SystemDB(filename="file:mbe_aux?mode=memory&cache=shared")
    yield db
    db.close()


SHIFT = 0.004  # eV added to every energy at the explicit-aux level
KJ_PER_EV = 96.48533212331002


def evaluate(self, node, structures, stress=False):
    """The shell levels are the high level; the explicit-aux copy is shifted."""
    if self.name.startswith("high_shells"):
        out = evaluate_with_bsse(types.SimpleNamespace(name="high"), node, structures)
        if self.name.endswith("_aux"):
            for result in out.values():
                result.energy += SHIFT * KJ_PER_EV
        return out
    return evaluate_with_bsse(self, node, structures, stress)


SHELL_TZ = {**SETTINGS, "shell high level": "TZ"}


def shell_pair(node, P, conf):
    frame = node._frame(conf, node._rules(P), P)
    fragments = frame["fragments"]
    pairs = [
        f for f in fragments.selected() if f.order == 2 and holds_shell(fragments, f)
    ]
    return frame, pairs[0]


def test_parse():
    assert levels.parse_aux_fragments("") == {}
    assert levels.parse_aux_fragments(
        "high:shell c1-d01_09 c1-d00_02; molecular:shell c1-d01_09"
    ) == {"high:shell": {"c1-d01_09", "c1-d00_02"}, "molecular:shell": {"c1-d01_09"}}
    with pytest.raises(ValueError, match="not one of"):
        levels.parse_aux_fragments("low:shell c1-d01_09")
    assert levels.parse_aux_sets("explicit") == {}
    assert levels.parse_aux_sets("high:shell def2/J,def2-TZVPP/C") == {
        "high:shell": "def2/J,def2-TZVPP/C"
    }
    with pytest.raises(ValueError, match="separated by commas"):
        levels.parse_aux_sets("high:shell def2/J def2-TZVPP/C")


def test_only_the_listed_increment_moves(db, monkeypatch):
    monkeypatch.setattr(levels.Level, "evaluate", evaluate)
    conf = li_cluster(db)
    node, P = parameters(**SHELL_TZ)
    frame, pair = shell_pair(node, P, conf)
    _, plain = run(node, P, frame, _levels(**{"high:shell": "high_shells"}))
    increments = dict(frame["correction"].increments)

    key = f"{frame['prefix']}{pair.name}"
    node, P = parameters(
        **SHELL_TZ, **{"auxiliary basis fragments": f"high:shell {key}"}
    )
    frame = node._frame(conf, node._rules(P), P)
    fragments = frame["fragments"]
    calcs = frame["calculations"]
    assert frame["variants"] == {"high:shell": ("high:shell@aux", {pair.name})}
    assert calcs["high:shell@aux"] == fragments.closure([pair])
    assert pair.name not in calcs["high:shell"]
    levels_ = _levels(
        **{"high:shell": "high_shells", "high:shell@aux": "high_shells_aux"}
    )
    _, labels = run(node, P, frame, levels_)
    # the pair's increment is (E + c) - 2 (E_m + c) = -c off; nothing else moves
    assert labels.energy - plain.energy == pytest.approx(-SHIFT, abs=1e-9)
    for name, inc in frame["correction"].increments.items():
        delta = inc.energy - increments[name].energy
        if name == pair.name:
            assert inc.high_level == "high:shell@aux"
            assert delta == pytest.approx(-SHIFT, abs=1e-9)
        else:
            assert inc.high_level != "high:shell@aux"
            assert delta == pytest.approx(0.0, abs=1e-12)


def orca_resolve(text, context, current=None, periodic=False):
    return {"level": text, "step": "ORCA", "options": {"prefers_batch": True}}


def test_the_aux_level_is_a_copy_with_the_sets(db, monkeypatch):
    monkeypatch.setattr(levels, "resolve", orca_resolve)
    monkeypatch.setattr(mbe_step.Mbe, "variable_exists", lambda self, name: False)
    conf = li_cluster(db)
    settings = {**SHELL_TZ, "high level": "QZ", "molecular low level": "QZLOW"}
    node, P = parameters(**settings)
    frame, pair = shell_pair(node, P, conf)
    key = f"{frame['prefix']}{pair.name}"
    node, P = parameters(
        **settings, **{"auxiliary basis fragments": f"high:shell {key}"}
    )
    frame = node._frame(conf, node._rules(P), P)
    levels_ = node._levels(P, {}, [frame])
    aux = levels_["high:shell@aux"]
    assert aux.name == "high_shells_aux"
    assert aux.options["aux_basis"] == "explicit"
    assert "aux_basis" not in levels_["high:shell"].options
    node, P = parameters(
        **settings,
        **{
            "auxiliary basis fragments": f"high:shell {key}",
            "auxiliary basis sets": "high:shell def2/J,def2-TZVPP/C",
        },
    )
    levels_ = node._levels(P, {}, [frame])
    assert levels_["high:shell@aux"].options["aux_basis"] == "def2/J,def2-TZVPP/C"


def test_refusals(db, monkeypatch):
    monkeypatch.setattr(mbe_step.Mbe, "variable_exists", lambda self, name: False)
    conf = li_cluster(db)
    node, P = parameters(**SHELL_TZ)
    frame, pair = shell_pair(node, P, conf)
    prefix = frame["prefix"]

    def frame_with(text):
        node, P = parameters(
            **SHELL_TZ,
            **{
                "high level": "QZ",
                "molecular low level": "QZLOW",
                "auxiliary basis fragments": text,
            },
        )
        return node, P, node._frame(conf, node._rules(P), P)

    with pytest.raises(ValueError, match="is a monomer"):
        frame_with(f"high:shell {prefix}m00")
    with pytest.raises(ValueError, match="is not a fragment"):
        frame_with(f"high:shell {prefix}d77_99")
    with pytest.raises(ValueError, match="its increment uses"):
        frame_with(f"high {prefix}{pair.name}")
    node, P, frame = frame_with(f"high:shell c999-{pair.name}")
    with pytest.raises(ValueError, match="not in any of the configurations"):
        node._check_aux_fragments(P, [frame])
    # not ORCA as tasks
    monkeypatch.setattr(
        levels,
        "resolve",
        lambda text, *a, **k: {"level": text, "step": "fake", "options": {}},
    )
    node, P, frame = frame_with(f"high:shell {prefix}{pair.name}")
    with pytest.raises(ValueError, match="cannot take explicit auxiliary"):
        node._levels(P, {}, [frame])


def evaluate_monomer_aux(self, node, structures, stress=False):
    """As evaluate, and the high level's explicit-aux copy is the high level shifted."""
    if self.name == "high_aux":
        out = evaluate_with_bsse(types.SimpleNamespace(name="high"), node, structures)
        for result in out.values():
            result.energy += SHIFT * KJ_PER_EV
        return out
    return evaluate(self, node, structures, stress)


def test_a_listed_1_body_term_moves_alone(db, monkeypatch):
    """'Explicit auxiliary basis for 1-body terms' runs a unit's own 1-body term
    at the explicit-aux copy of its level (Li-salt shells where AutoAux fails at
    revDSD/def2-QZVPPD, 2026-10-29): only that term moves, by the copy's offset."""
    monkeypatch.setattr(levels.Level, "evaluate", evaluate_monomer_aux)
    conf = li_cluster(db)
    node, P = parameters(**SHELL_TZ)
    frame = node._frame(conf, node._rules(P), P)
    fragments = frame["fragments"]
    (shell,) = [f for f in fragments.by_order(1) if holds_shell(fragments, f)]
    _, plain = run(node, P, frame, _levels(**{"high:shell": "high_shells"}))
    increments = dict(frame["correction"].increments)

    key = f"{frame['prefix']}{shell.name}"
    node, P = parameters(**SHELL_TZ, **{"auxiliary basis monomers": f"high {key}"})
    node._id = (1,)
    text = " ".join(node.description_text(node.parameters.values_to_dict()).split())
    assert "not on the AutoAux scale" in text
    frame = node._frame(conf, node._rules(P), P)
    calcs = frame["calculations"]
    assert frame["variants"] == {"high": ("high@aux", {shell.name})}
    assert calcs["high@aux"] == [shell.name]
    assert shell.name not in calcs["high"]
    levels_ = _levels(**{"high:shell": "high_shells", "high@aux": "high_aux"})
    _, labels = run(node, P, frame, levels_)
    assert labels.energy - plain.energy == pytest.approx(SHIFT, abs=1e-9)
    for name, inc in frame["correction"].increments.items():
        delta = inc.energy - increments[name].energy
        if name == shell.name:
            assert inc.high_level == "high@aux"
            assert delta == pytest.approx(SHIFT, abs=1e-9)
        else:
            assert delta == pytest.approx(0.0, abs=1e-12)


def test_1_body_entries_hold_only_1_body_terms(db, monkeypatch):
    monkeypatch.setattr(mbe_step.Mbe, "variable_exists", lambda self, name: False)
    conf = li_cluster(db)
    node, P = parameters(**SHELL_TZ)
    frame, pair = shell_pair(node, P, conf)
    prefix = frame["prefix"]
    node, P = parameters(
        **SHELL_TZ, **{"auxiliary basis monomers": f"high:shell {prefix}{pair.name}"}
    )
    with pytest.raises(ValueError, match="is a fragment of 2 units"):
        node._frame(conf, node._rules(P), P)
    # a monomer under the increments' entry is still refused, with the pointer
    node, P = parameters(
        **SHELL_TZ, **{"auxiliary basis fragments": f"high {prefix}m00"}
    )
    with pytest.raises(ValueError, match="Explicit auxiliary basis for 1-body terms"):
        node._frame(conf, node._rules(P), P)
