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
