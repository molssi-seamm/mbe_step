"""The step's helpers: tables, offsets, the stress layout, the selection rules
from the parameters."""

import numpy as np
import pytest

from mbe_step import labels as labels_

from .helpers import parameters


def test_parse_table():
    table = labels_.parse_table("water water 4.5; Li+ * 3.0\n* * 5")
    assert table == {("water", "water"): 4.5, ("Li+", "*"): 3.0, ("*", "*"): 5.0}
    with pytest.raises(ValueError, match="'type type value'"):
        labels_.parse_table("water 4.5")
    with pytest.raises(ValueError, match="not a number"):
        labels_.parse_table("water water x")
    with pytest.raises(ValueError, match="empty"):
        labels_.parse_table(" ; ")


def test_parse_offsets():
    assert labels_.parse_offsets("none") is None
    assert labels_.parse_offsets("") is None
    assert labels_.parse_offsets("water 2074.69325; Li+ 200") == {
        "water": 2074.69325,
        "Li+": 200.0,
    }
    with pytest.raises(ValueError, match="'type value'"):
        labels_.parse_offsets("water")


def test_voigt_order_and_symmetry():
    t = np.array([[1.0, 6.0, 5.0], [6.2, 2.0, 4.0], [5.0, 4.0, 3.0]])
    v = labels_.voigt(t)
    assert v[:3] == [1.0, 2.0, 3.0]
    assert v[3] == 4.0 and v[4] == 5.0 and v[5] == pytest.approx(6.1)


def test_rules_from_parameters():
    node, P = parameters()
    rules = node._rules(P)
    assert rules.max_order == 3
    assert rules.cutoffs == {2: 4.5, 3: 3.5}
    assert rules.rules == {3: "connected"}
    assert rules.criterion == "designated"

    node, P = parameters(
        **{
            "maximum order": "2",
            "cutoffs": "table by type pair",
            "pair cutoff table": "water water 4.0; * * 5.0",
            "distance criterion": "closest heavy-atom contact",
        }
    )
    rules = node._rules(P)
    assert rules.max_order == 2
    assert rules.cutoffs == {2: {("water", "water"): 4.0, ("*", "*"): 5.0}}
    assert rules.criterion == "heavy contact"

    node, P = parameters(**{"pair cutoff": 5.0, "maximum order": "1"})
    assert node._rules(P).cutoffs == {}


def test_description_text():
    node, P = parameters(**{"molecular low level": "ORCA:DFT@r2SCAN/def2-SVP"})
    node._id = (1,)
    text = " ".join(node.description_text(node.parameters.values_to_dict()).split())
    assert "monomers, pairs and triples" in text
    assert "ORCA:DFT@r2SCAN/def2-SVP" in text


def test_extxyz_rewrites_a_configuration_instead_of_duplicating(tmp_path):
    from types import SimpleNamespace

    system = SimpleNamespace(
        periodic=False,
        symbols=["O", "H", "H"],
        coordinates=np.zeros((3, 3)),
    )

    def labels(e):
        return SimpleNamespace(reference_energy=e, forces=np.zeros((3, 3)))

    path = tmp_path / "labels.extxyz"
    labels_.write_extxyz(path, system, labels(-1.0), name="a", model="m", identifier=1)
    labels_.write_extxyz(path, system, labels(-2.0), name="b", model="m", identifier=2)
    labels_.write_extxyz(path, system, labels(-3.0), name="a", model="m", identifier=1)
    frames = labels_._read_frames(path)
    assert len(frames) == 2
    assert "REF_energy=-2.00000000" in frames[0][1]
    assert "REF_energy=-3.00000000" in frames[1][1]
