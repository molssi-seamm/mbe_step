"""The levels that run as queued tasks run at the same time, the long periodic
ones started first; MDI levels run one at a time."""

import threading

import numpy as np
import pytest

from mbe_step import levels

from .helpers import level, parameters
from .test_phase3 import cluster, evaluate_with_bsse, run

SETTINGS = {"pair cutoff": 100.0, "triple cutoff": 100.0, "triple rule": "compact"}


@pytest.fixture
def db():
    from molsystem import SystemDB

    db = SystemDB(filename="file:mbe_concurrent?mode=memory&cache=shared")
    yield db
    db.close()


def _levels():
    return {
        "high": level("high"),
        "molecular": level("molecular"),
        "cluster": level("cluster"),
    }


def test_batch_levels_run_at_the_same_time(db, monkeypatch):
    """Each batch level waits at a barrier that only opens when all three are
    inside evaluate() at once: a serial run would time out."""
    barrier = threading.Barrier(3, timeout=10)
    started = []

    def evaluate(self, node, structures, stress=False):
        started.append(self.name)
        barrier.wait()
        return evaluate_with_bsse(self, node, structures, stress)

    monkeypatch.setattr(levels.Level, "evaluate", evaluate)
    monkeypatch.setattr(levels.Level, "runs_as_tasks", lambda self, node: True)
    conf = cluster(db)
    node, P = parameters(**SETTINGS)
    frame = node._frame(conf, node._rules(P), P)
    results, labels = run(node, P, frame, _levels())
    assert set(started) == {"high", "molecular", "cluster"}
    assert set(results) == {"high", "molecular", "cluster"}

    # The same labels as a one-at-a-time run
    monkeypatch.setattr(levels.Level, "evaluate", evaluate_with_bsse)
    monkeypatch.setattr(levels.Level, "runs_as_tasks", lambda self, node: False)
    _, serial = run(node, P, frame, _levels())
    assert labels.energy == serial.energy
    assert np.array_equal(labels.forces, serial.forces)


def test_mdi_levels_run_one_at_a_time(db, monkeypatch):
    active, peak = [0], [0]
    lock = threading.Lock()

    def evaluate(self, node, structures, stress=False):
        with lock:
            active[0] += 1
            peak[0] = max(peak[0], active[0])
        try:
            return evaluate_with_bsse(self, node, structures, stress)
        finally:
            with lock:
                active[0] -= 1

    monkeypatch.setattr(levels.Level, "evaluate", evaluate)
    monkeypatch.setattr(levels.Level, "runs_as_tasks", lambda self, node: False)
    conf = cluster(db)
    node, P = parameters(**SETTINGS)
    frame = node._frame(conf, node._rules(P), P)
    run(node, P, frame, _levels())
    assert peak[0] == 1


def test_the_long_levels_start_first():
    order = [
        name
        for name, _ in levels_order(
            {"high": 1, "high:3": 2, "molecular": 3, "periodic": 4, "cell": 5}
        )
    ]
    assert order[:2] == ["cell", "periodic"]
    assert order[2:] == ["high", "high:3", "molecular"]


def levels_order(levels_):
    import mbe_step

    return mbe_step.Mbe._evaluation_order(levels_)


# ------------------------------------------------------------- fail fast
def _waiting_evaluate(failing):
    """An evaluate that returns at once for the ``failing`` level (raising or
    with every result failed) and otherwise waits until it is cancelled, as a
    long VASP level would."""
    import time
    from types import SimpleNamespace

    def evaluate(self, node, structures, stress=False):
        if self.name == failing[0]:
            if failing[1] == "raise":
                raise RuntimeError("bad input for this level")
            return {
                key: SimpleNamespace(key=key, ok=False, reason="no task")
                for key in structures
            }
        deadline = time.time() + 20
        while not self._cancel_requested:
            if time.time() > deadline:
                raise AssertionError(f"{self.name} was never cancelled")
            time.sleep(0.01)
        return {
            key: SimpleNamespace(key=key, ok=False, reason="cancelled")
            for key in structures
        }

    return evaluate


@pytest.mark.parametrize("how", ["raise", "all failed"])
def test_a_failed_level_stops_the_others(db, monkeypatch, how):
    monkeypatch.setattr(levels.Level, "evaluate", _waiting_evaluate(("molecular", how)))
    monkeypatch.setattr(levels.Level, "runs_as_tasks", lambda self, node: True)
    conf = cluster(db)
    node, P = parameters(**SETTINGS)
    frame = node._frame(conf, node._rules(P), P)
    levels_ = _levels()
    if how == "raise":
        with pytest.raises(RuntimeError, match="bad input"):
            node._evaluate(levels_, [frame])
    else:
        results = node._evaluate(levels_, [frame])
        assert {r.reason for r in results["high"].values()} == {"cancelled"}
        assert {r.reason for r in results["cluster"].values()} == {"cancelled"}
    assert levels_["high"]._cancel_requested
    assert levels_["cluster"]._cancel_requested
    assert not levels_["molecular"]._cancel_requested


def test_a_serial_failure_skips_the_remaining_levels(db, monkeypatch):
    calls = []

    def evaluate(self, node, structures, stress=False):
        calls.append(self.name)
        if self.name == "cluster":
            raise RuntimeError("bad input for this level")
        return evaluate_with_bsse(self, node, structures, stress)

    monkeypatch.setattr(levels.Level, "evaluate", evaluate)
    monkeypatch.setattr(levels.Level, "runs_as_tasks", lambda self, node: False)
    conf = cluster(db)
    node, P = parameters(**SETTINGS)
    frame = node._frame(conf, node._rules(P), P)
    with pytest.raises(RuntimeError, match="bad input"):
        node._evaluate(_levels(), [frame])
    # The cluster runs first; the molecular levels never start
    assert calls == ["cluster"]


def test_cancel_reaches_a_running_evaluator():
    from types import SimpleNamespace

    cancelled = []
    lv = level("high")
    lv._evaluator = SimpleNamespace(cancel=lambda: cancelled.append(True))
    lv.cancel()
    assert cancelled == [True] and lv._cancel_requested
