"""The pilot-frame regression through the step's own path.

The prototype's stored fragment results are fed to the step as the
Evaluator-shaped results it consumes (kJ/mol, gradients in kJ/mol/Å, the
cell's stress in GPa in the pressure convention), keyed as the step keys them.
The step builds the frame from a molsystem configuration, submits the
structures, converts the units, assembles, and writes the labels; the numbers
must match the prototype (seamm_mbe's regression, tests/data/README.md).
"""

import json
from pathlib import Path

import numpy as np
import pytest

import seamm_mbe
from seamm_mbe import units
from mbe_step import labels as labels_
from mbe_step import levels

from .helpers import configuration, level, parameters, result

DATA = Path(__file__).parent / "data"
PROTOTYPE_KB_ATM = 986.923
EXACT_KB_ATM = 1e8 / 101325


@pytest.fixture(scope="module")
def pilot():
    data = np.load(DATA / "water64_pilot.npz")
    meta = json.loads(str(data["meta"]))
    offsets = np.cumsum([0] + [3 * len(m["mols"]) for m in meta])
    slices = {m["name"]: slice(offsets[k], offsets[k + 1]) for k, m in enumerate(meta)}
    index = {m["name"]: k for k, m in enumerate(meta)}
    expected = json.load(open(DATA / "water64_pilot_expected.json"))
    return data, meta, slices, index, expected


def stored_results(pilot):
    """A fake Level.evaluate serving the stored results."""
    data, meta, slices, index, _ = pilot
    box = float(data["box"])
    volume = box**3
    w = units.vasp_stress_to_virial(
        data["cell_vasp_stress_kB"], volume
    ) + units.dftd4_virial_to_virial(data["cell_d4_dEdstrain"])
    cell_pressure = units.tensor_from_virial(
        w, volume, convention="pressure", units="GPa"
    )

    def evaluate(self, node, structures, stress=False):
        out = {}
        for key in structures:
            name = key.split("-", 1)[1]
            if name == "cell":
                e = float(data["cell_vasp_energy"] + data["cell_d4_energy"])
                f = data["cell_vasp_forces"] + data["cell_d4_forces"]
                out[key] = result(key, e, f, cell_pressure)
                continue
            k, s = index[name], slices[name]
            if self.name == "high":
                e, f = data["high_E"][k], data["high_F"][s]
            elif self.name == "molecular":
                e, f = data["molecular_E"][k], data["molecular_F"][s]
            else:
                e = data["periodic_E"][k] + data["d4_E"][k]
                f = data["periodic_F"][s] + data["d4_F"][s]
            out[key] = result(key, float(e), f)
        return out

    return evaluate


@pytest.fixture
def db():
    from molsystem import SystemDB

    db = SystemDB(filename="file:mbe_pilot?mode=memory&cache=shared")
    yield db
    db.close()


def test_pilot_through_the_step(pilot, db, monkeypatch, tmp_path):
    data, meta, _, _, expected = pilot
    monkeypatch.setattr(levels.Level, "evaluate", stored_results(pilot))
    box = float(data["box"])
    X = data["X"]
    conf = configuration(
        db,
        ["O", "H", "H"] * 64,
        X,
        cell=np.eye(3) * box,
        bonds=[(3 * m, 3 * m + k) for m in range(64) for k in (1, 2)],
    )
    node, P = parameters(**{"periodic low level": "fake periodic"})
    frame = node._frame(conf, node._rules(P), P)
    assert frame["fragments"].names == [m["name"] for m in meta]
    calcs = frame["calculations"]
    assert (
        len(calcs["periodic"]) + 1,
        len(calcs["high"]) + len(calcs["molecular"]),
    ) == (
        219,
        2564,
    )
    levels_ = {
        "high": level("high"),
        "molecular": level("molecular"),
        "periodic": level("periodic"),
        "cell": level("cell", convention="pressure"),
    }
    results = node._evaluate(levels_, [frame])
    offsets = labels_.parse_offsets(P["energy offsets"])
    labels = node._assemble(frame, levels_, results, offsets)

    debug = expected["debug"]
    shift = debug["P_vasp_atm"] * (EXACT_KB_ATM / PROTOTYPE_KB_ATM - 1)
    assert labels.breakdown["MBE"]["pressure"] == pytest.approx(
        expected["P_correction_atm"], abs=1e-5
    )
    assert labels.pressure == pytest.approx(
        expected["P_atom_conf_atm"] + shift, abs=1e-5
    )
    assert labels.molecular_pressure == pytest.approx(
        expected["P_mol_conf_atm"] + shift, abs=1e-5
    )
    for order in (1, 2, 3):
        assert labels.per_body[order]["pressure"] == pytest.approx(
            debug["P_per_body_atm"][order - 1], abs=1e-5
        )
    assert labels.reference_energy == pytest.approx(expected["REF_energy_eV"], abs=1e-7)
    assert labels.max_increment_net_force * 1000 == pytest.approx(
        expected["max_increment_net_force_meV_A"], abs=1e-8
    )

    # SEAMM's units and the labels file
    seamm = labels_.to_seamm(labels, box**3)
    assert len(seamm["stress"]) == 6
    p_from_voigt = -np.mean(seamm["stress"][:3]) * 1e9 / 101325
    assert p_from_voigt == pytest.approx(labels.pressure, abs=1e-3)
    path = tmp_path / "labels.extxyz"
    labels_.write_extxyz(path, frame["system"], labels, name="pilot", model="test")
    lines = path.read_text().splitlines()
    assert lines[0] == "192"
    header = lines[1]
    assert f"REF_energy={labels.reference_energy:.8f}" in header
    stress = [float(v) for v in header.split('REF_stress="')[1].split('"')[0].split()]
    assert np.allclose(np.array(stress).reshape(3, 3), labels.stress, atol=1e-12)
    row = lines[2].split()
    assert row[0] == "O" and float(row[4]) == pytest.approx(
        labels.forces[0, 0], abs=1e-8
    )


def test_periodic_level_must_declare_stress_convention():
    with pytest.raises(ValueError, match="stress_convention"):
        levels.stress_convention({"level": "X:DFT@PBE", "options": {}})
    assert levels.stress_convention({"options": {"stress_convention": "stress"}}) == (
        "stress"
    )
    assert seamm_mbe.units.CONVENTIONS == ("pressure", "stress")
