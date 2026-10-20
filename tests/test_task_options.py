"""Options for every task of a level, e.g. RIJCOSX kept for the shell pairs holding
B or P (``high:shell exact_exchange_elements=none``)."""

import numpy as np
import pytest

import mbe_step
from mbe_step import levels

from .helpers import parameters
from .test_shell_levels import SETTINGS, fake_resolve, li_cluster


@pytest.fixture
def db():
    from molsystem import SystemDB

    db = SystemDB(filename="file:mbe_options?mode=memory&cache=shared")
    yield db
    db.close()


def test_parse_task_options():
    assert levels.parse_task_options("") == {}
    assert levels.parse_task_options("none") == {}
    assert levels.parse_task_options(
        "high:shell exact_exchange_elements=none; high a=1 b=x=y"
    ) == {
        "high:shell": {"exact_exchange_elements": "none"},
        "high": {"a": "1", "b": "x=y"},
    }
    for bad in ("nowhere a=1", "high", "high a", "high =1"):
        with pytest.raises(ValueError):
            levels.parse_task_options(bad)


def _levels_for(db, monkeypatch, **values):
    monkeypatch.setattr(levels, "resolve", fake_resolve)
    monkeypatch.setattr(mbe_step.Mbe, "variable_exists", lambda self, name: False)
    conf = li_cluster(db)
    settings = {**SETTINGS, "high level": "QZ", "molecular low level": "QZLOW"}
    node, P = parameters(**settings, **values)
    frame = node._frame(conf, node._rules(P), P)
    return node._levels(P, {}, [frame])


def test_options_go_to_their_level(db, monkeypatch):
    levels_ = _levels_for(
        db,
        monkeypatch,
        **{
            "triples high level": "TZ",
            "shell high level": "TZ",
            "task options": "high:shell exact_exchange_elements=none",
        },
    )
    # the shell level is the triples' Level, which takes the options
    assert levels_["high:shell"] is levels_["high:3"]
    assert levels_["high:3"].options == {"exact_exchange_elements": "none"}
    assert levels_["high"].options == {}
    assert levels_["molecular"].options == {}


def test_options_for_an_unused_level_are_refused(db, monkeypatch):
    with pytest.raises(ValueError, match="does not use"):
        _levels_for(db, monkeypatch, **{"task options": "high:shell a=1"})


def test_conflicting_options_on_a_shared_level_are_refused(db, monkeypatch):
    with pytest.raises(ValueError, match="two values"):
        _levels_for(
            db,
            monkeypatch,
            **{
                "triples high level": "TZ",
                "shell high level": "TZ",
                "task options": "high:3 a=1; high:shell a=2",
            },
        )


def test_each_task_gets_the_options(monkeypatch):
    """The level's options join each structure's own, which take precedence."""
    submitted = []

    class FakeEvaluator:
        def __init__(self, *args, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def cancel(self):
            pass

        def submit(self, structure, key=None, options=None):
            submitted.append((key, options))

        def results(self):
            return iter(())

    import seamm_exec

    monkeypatch.setattr(seamm_exec, "Evaluator", FakeEvaluator)
    from seamm_exec import Geometry

    lv = levels.Level(
        "high_shells",
        {"level": "x"},
        ranks=2,
        memory=3800,
        bundle=8,
        walltime=0,
        archive=False,
        options={"exact_exchange_elements": "none", "a": "level"},
    )
    node = type("Node", (), {"directory": "."})()
    g = Geometry([1, 1], np.zeros((2, 3)))
    lv.evaluate(node, {"plain": g, "own": (g, {"a": "structure", "ghost_atoms": [1]})})
    assert dict(submitted) == {
        "plain": {"exact_exchange_elements": "none", "a": "level"},
        "own": {
            "exact_exchange_elements": "none",
            "a": "structure",
            "ghost_atoms": [1],
        },
    }
