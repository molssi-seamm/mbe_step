"""Energy offsets onto the DfE0 scale from the thermochemistry database."""

import sys
import types

import pytest

import seamm_mbe
from mbe_step import thermo

from .helpers import parameters

MC = {
    "level": "ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-QZVPPD",
    "step": "ORCA",
    "method": "REVDSD-PBEP86-D4_2021",
    "basis": "def2-QZVPPD",
}


def ec():
    return seamm_mbe.System(
        ["O", "C", "O", "C", "C", "O", "H", "H", "H", "H"],
        [
            [0, 0, 0],
            [1.2, 0, 0],
            [1.9, 1.1, 0],
            [3.3, 0.8, 0],
            [3.3, -0.7, 0],
            [1.9, -1.1, 0],
            [3.8, 1.2, 0.9],
            [3.8, 1.2, -0.9],
            [3.8, -1.1, 0.9],
            [3.8, -1.1, -0.9],
        ],
        bonds=[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 1)]
        + [(3, 6), (3, 7), (4, 8), (4, 9)],
    )


def test_automatic_is_the_default():
    node, P = parameters()
    assert thermo.is_automatic(P["energy offsets"])
    assert not thermo.is_automatic("none")
    assert not thermo.is_automatic("water 2074.69325")
    node._id = (1,)
    assert "DfE0" in node.description_text(P)


def test_the_reference_of_a_model_chemistry():
    assert thermo.reference(MC) == ("orca", "REVDSD-PBEP86-D4_2021", "def2-QZVPPD")
    assert thermo.reference({**MC, "basis": "bse:def2-QZVPPD"})[2] == "def2-QZVPPD"
    hook = types.SimpleNamespace(thermochemistry_reference=lambda mc: ("x", "y", "z"))
    assert thermo.reference(MC, hook) == ("x", "y", "z")


@pytest.fixture
def fake_thermo(monkeypatch):
    """A stand-in for seamm_thermochemistry: DfE0 - E = -100 kJ/mol per atom,
    and no atoms at all for 'missing'."""
    calls = []

    class DB:
        def __init__(self, read_only=True):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def formation_energy(composition, energy, db, code, method, **kw):
        calls.append((dict(composition), code, method, kw["settings"]))
        if method == "missing":
            raise KeyError("no atom")
        assert kw["anchor"] and kw["anchor_at_0K"] and kw["units"] == "kJ/mol"
        return energy - 100.0 * sum(composition.values())

    module = types.SimpleNamespace(ThermoDB=DB, formation_energy=formation_energy)
    monkeypatch.setitem(sys.modules, "seamm_thermochemistry", module)
    return calls


def test_offsets_per_type(fake_thermo):
    offsets = thermo.dfe0_offsets(MC, [ec()])
    assert list(offsets) == ["EC"]
    assert offsets["EC"] == pytest.approx(
        -1000.0 / seamm_mbe.units.EV_TO_KJ_PER_MOL, rel=1e-12
    )
    composition, code, method, settings = fake_thermo[0]
    assert composition == {"C": 3, "O": 3, "H": 4}
    assert (code, method, settings) == ("orca", MC["method"], "def2-QZVPPD")


def test_a_level_without_atoms_is_refused(fake_thermo):
    with pytest.raises(ValueError, match="No DfE0 energy offset for EC"):
        thermo.dfe0_offsets({**MC, "method": "missing"}, [ec()])


def test_without_seamm_thermochemistry(monkeypatch):
    monkeypatch.setitem(sys.modules, "seamm_thermochemistry", None)
    with pytest.raises(ValueError, match="need seamm_thermochemistry"):
        thermo.dfe0_offsets(MC, [ec()])


def test_the_real_database():
    """science's values (2026-10-06): EC 9300.36250 eV per molecule at
    revDSD/def2-QZVPPD, 9299.16652 at def2-TZVPPD."""
    try:
        offsets = thermo.dfe0_offsets(MC, [ec()])
    except ValueError as e:
        pytest.skip(f"No thermochemistry database here: {e}")
    assert offsets["EC"] == pytest.approx(9300.3625, abs=1e-4)
    tz = thermo.dfe0_offsets({**MC, "basis": "def2-TZVPPD"}, [ec()])
    assert tz["EC"] == pytest.approx(9299.16652, abs=1e-4)


def test_without_the_database_file(monkeypatch):
    """seamm_thermochemistry installed by pip alone has no database file."""

    class DB:
        def __init__(self, read_only=True):
            raise FileNotFoundError("/somewhere/thermochemistry.db")

    module = types.SimpleNamespace(ThermoDB=DB, formation_energy=None)
    monkeypatch.setitem(sys.modules, "seamm_thermochemistry", module)
    with pytest.raises(ValueError, match="database, which is not installed"):
        thermo.dfe0_offsets(MC, [ec()])
