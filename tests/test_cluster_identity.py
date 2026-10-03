"""The end-to-end identity check through the step's own code.

For a cluster, the step's labels are E_low(cluster) + sum of increments. With
every pair and triple selected and a high level whose interactions are at most
3-body, the expansion truncated at triples is exact: the labels must equal the
high level on the whole cluster, energy and forces, to rounding.

The levels are analytic potentials behind a stand-in for the Evaluator, so the
test exercises the step's frame building, keys, structure submission, unit
conversion and assembly without running a quantum-chemistry code.
"""

import numpy as np
import pytest
from seamm_exec import structure_data

from mbe_step import levels

from .helpers import configuration, level, model_potential, parameters, result

WATERS = np.array(
    [
        [0.000, 0.000, 0.000],
        [0.957, 0.000, 0.000],
        [-0.240, 0.927, 0.000],
        [2.900, 0.200, 0.100],
        [3.400, 1.000, 0.200],
        [3.500, -0.500, 0.300],
        [0.300, 2.800, 0.500],
        [1.200, 3.100, 0.600],
        [-0.200, 3.600, 0.700],
        [1.500, 1.500, 2.700],
        [2.300, 1.400, 3.200],
        [0.900, 1.900, 3.400],
    ]
)


def charges_for(z, scale):
    return np.where(np.asarray(z) == 8, -0.8, 0.4) * scale


def fake_evaluate(self, node, structures, stress=False):
    """Each level's potential: high = charges + a 3-body term; the low levels
    = scaled charges, no 3-body term."""
    out = {}
    for key, structure in structures.items():
        data = structure_data(structure)
        z = data["atomic_numbers"]
        xyz = np.asarray(data["coordinates"], dtype=float)
        if self.name == "high":
            e, f = model_potential(z, xyz, charges_for(z, 1.0), c3=0.5)
        else:
            e, f = model_potential(z, xyz, charges_for(z, 0.7), c3=0.0)
        out[key] = result(key, e, f)
    return out


@pytest.fixture
def db():
    from molsystem import SystemDB

    db = SystemDB(filename="file:mbe_identity?mode=memory&cache=shared")
    yield db
    db.close()


def test_cluster_identity(db, monkeypatch, tmp_path):
    monkeypatch.setattr(levels.Level, "evaluate", fake_evaluate)
    conf = configuration(
        db,
        ["O", "H", "H"] * 4,
        WATERS,
        bonds=[(3 * m, 3 * m + k) for m in range(4) for k in (1, 2)],
    )
    node, P = parameters(
        **{
            "pair cutoff": 100.0,
            "triple cutoff": 100.0,
            "triple rule": "compact",
            "energy offsets": "none",
        }
    )
    rules = node._rules(P)
    frame = node._frame(conf, rules, P)
    counts = frame["fragments"].counts()
    assert counts[2]["selected"] == 6 and counts[3]["selected"] == 4
    levels_ = {
        "high": level("high"),
        "molecular": level("molecular"),
        "cluster": level("cluster"),
    }
    results = node._evaluate(levels_, [frame])
    labels = node._assemble(frame, levels_, results, None)

    z = [8, 1, 1] * 4
    e_high, f_high = model_potential(z, WATERS, charges_for(z, 1.0), c3=0.5)
    assert labels.energy == pytest.approx(e_high, abs=1e-9)
    assert np.abs(labels.forces - f_high).max() < 1e-8
    assert labels.virial is None
    # The 3-body term shows up only in the 3-body increments
    assert abs(labels.per_body[3]["energy"]) > 1e-6


def test_missing_calculation_gives_no_labels(db, monkeypatch, tmp_path):
    def failing(self, node, structures, stress=False):
        out = fake_evaluate(self, node, structures, stress)
        if self.name == "molecular":
            key = next(k for k in out if k.endswith("d00_01"))
            out[key].ok = False
            out[key].reason = "SCF did not converge"
        return out

    monkeypatch.setattr(levels.Level, "evaluate", failing)
    conf = configuration(
        db,
        ["O", "H", "H"] * 4,
        WATERS,
        bonds=[(3 * m, 3 * m + k) for m in range(4) for k in (1, 2)],
    )
    node, P = parameters(
        **{"pair cutoff": 100.0, "triple cutoff": 100.0, "triple rule": "compact"}
    )
    frame = node._frame(conf, node._rules(P), P)
    levels_ = {
        "high": level("high"),
        "molecular": level("molecular"),
        "cluster": level("cluster"),
    }
    results = node._evaluate(levels_, [frame])
    import seamm_mbe

    with pytest.raises(seamm_mbe.MissingFragmentsError) as error:
        node._assemble(frame, levels_, results, None)
    text = node._failure_text(frame, error.value, results)
    assert "d00_01 (molecular): SCF did not converge" in text
