# -*- coding: utf-8 -*-

"""Helpers of the MBE step that need no flowchart: parsing the cutoff and
offset tables, the stress layout, and the extended XYZ labels file."""

import numpy as np

import seamm_mbe


def parse_shells(text):
    """The ion-shell cutoffs.

    Parameters
    ----------
    text : str
        Entries 'ion-element partner-element cutoff' separated by ';' or new
        lines, e.g. 'Li O 2.6; Li F 2.6'.

    Returns
    -------
    {str: {str: float}}
        Per ion element, the cutoff (Å) to each partner element.
    """
    table = {}
    for (ion, partner), value in parse_table(text, "shell cutoff").items():
        if "*" in (ion, partner):
            raise ValueError(
                f"The shell cutoffs need explicit elements, not '*' ({ion} {partner})."
            )
        table.setdefault(ion, {})[partner] = value
    return table


def parse_table(text, what="cutoff"):
    """A table of values by molecule-type pair.

    Parameters
    ----------
    text : str
        Entries 'type type value' separated by ';' or new lines, e.g.
        'water water 4.5; Li+ * 3.0'. '*' matches any type.

    Returns
    -------
    {(str, str): float}
    """
    table = {}
    for entry in text.replace("\n", ";").split(";"):
        entry = entry.strip()
        if not entry:
            continue
        fields = entry.split()
        if len(fields) != 3:
            raise ValueError(
                f"The {what} table entry '{entry}' should be 'type type value', "
                "e.g. 'water water 4.5'."
            )
        try:
            value = float(fields[2])
        except ValueError:
            raise ValueError(f"The {what} in '{entry}' is not a number.") from None
        table[(fields[0], fields[1])] = value
    if not table:
        raise ValueError(f"The {what} table is empty.")
    return table


def parse_offsets(text):
    """Energy offsets in eV per molecule by type, from entries 'type value'
    separated by ';' or new lines, or None for 'none' or an empty text."""
    text = text.strip()
    if text.lower() in ("", "none"):
        return None
    offsets = {}
    for entry in text.replace("\n", ";").split(";"):
        entry = entry.strip()
        if not entry:
            continue
        fields = entry.split()
        if len(fields) != 2:
            raise ValueError(
                f"The energy offset '{entry}' should be 'type value', e.g. "
                "'water 2074.69325'."
            )
        try:
            offsets[fields[0]] = float(fields[1])
        except ValueError:
            raise ValueError(f"The offset in '{entry}' is not a number.") from None
    return offsets


def voigt(tensor):
    """A (3, 3) tensor as Voigt [6]: xx yy zz yz xz xy (symmetrized)."""
    t = np.asarray(tensor, dtype=float)
    t = (t + t.T) / 2
    return [t[0, 0], t[1, 1], t[2, 2], t[1, 2], t[0, 2], t[0, 1]]


def write_extxyz(
    path,
    system,
    labels,
    *,
    name,
    model,
    identifier=None,
    counterpoise=False,
    extra="",
):
    """Add one configuration's labels to an extended XYZ file.

    The file collects configurations (several in one step, or across a loop's
    iterations when it is a job-level file). Writing is idempotent: an earlier
    frame with the same ``identifier`` (SEAMM/configuration_id) is replaced, so
    rerunning a job does not duplicate labels.

    The columns are species, pos and REF_forces (eV/Å); the header carries
    REF_energy (eV, the reference scale), REF_stress (the nine values of
    sigma = -P in eV/Å³, for a cell), the Lattice and pbc, the model, whether
    counterpoise was used and the configuration's name. Note that ASE reads the
    nine-value REF_stress as a flat vector: a reader must reshape it.
    """
    fields = []
    if system.periodic:
        lattice = " ".join(f"{v:.8f}" for v in system.cell.ravel())
        fields.append(f'Lattice="{lattice}"')
    fields.append("Properties=species:S:1:pos:R:3:REF_forces:R:3")
    fields.append(f"REF_energy={labels.reference_energy:.8f}")
    if system.periodic:
        stress = " ".join(f"{v:.10e}" for v in labels.stress.ravel())
        fields.append(f'REF_stress="{stress}"')
        fields.append('pbc="T T T"')
    else:
        fields.append('pbc="F F F"')
    fields.append(f'model="{model}"')
    fields.append(f"counterpoise={'T' if counterpoise else 'F'}")
    fields.append(f'SEAMM/configuration_name="{name}"')
    if identifier is not None:
        fields.append(f"SEAMM/configuration_id={identifier}")
    if extra:
        fields.append(extra)
    lines = [str(len(system.symbols)), " ".join(fields)]
    for symbol, xyz, force in zip(system.symbols, system.coordinates, labels.forces):
        lines.append(
            f"{symbol:2s} {xyz[0]:14.8f} {xyz[1]:14.8f} {xyz[2]:14.8f} "
            f"{force[0]:14.8f} {force[1]:14.8f} {force[2]:14.8f}"
        )
    frames = _read_frames(path)
    if identifier is not None:
        tag = f"SEAMM/configuration_id={identifier}"
        frames = [f for f in frames if tag not in f[1].split()]
    frames.append(lines)
    with open(path, "w") as fd:
        for frame in frames:
            fd.write("\n".join(frame) + "\n")


def _read_frames(path):
    """The frames of an extended XYZ file as lists of lines ([] if none)."""
    try:
        text = open(path).read().splitlines()
    except FileNotFoundError:
        return []
    frames = []
    i = 0
    while i < len(text):
        if not text[i].strip():
            i += 1
            continue
        try:
            n = int(text[i])
        except ValueError:
            raise ValueError(
                f"{path}, line {i + 1}: expected the number of atoms of a frame, "
                f"found '{text[i][:40]}'. Is it an extended XYZ file?"
            ) from None
        frames.append(text[i : i + n + 2])
        i += n + 2
    return frames


def counts_text(fragments):
    """'64 monomers, 391 pairs (+268 outer), 559 triples'."""
    names = {1: "monomers", 2: "pairs", 3: "triples", 4: "4-bodies"}
    parts = []
    for order, count in sorted(fragments.counts().items()):
        text = f"{count['selected']} {names.get(order, f'{order}-bodies')}"
        if count["auxiliary"]:
            extra = "outer" if order == 2 else "auxiliary"
            text += f" (+{count['auxiliary']} {extra})"
        parts.append(text)
    return ", ".join(parts)


def to_seamm(labels, volume):
    """The labels in SEAMM's units for store_results: energies kJ/mol,
    gradients kJ/mol/Å, stress Voigt [6] in GPa as sigma = -P, pressures in
    atm."""
    units = seamm_mbe.units
    data = units.to_seamm(energy=labels.energy, forces=labels.forces)
    data["gradients"] = data["gradients"].tolist()
    data["reference energy"] = float(units.ev_to_kj_per_mol(labels.reference_energy))
    if labels.virial is not None:
        stress = units.tensor_from_virial(
            labels.virial, volume, convention="stress", units="GPa"
        )
        data["stress"] = voigt(stress)
        data["pressure"] = labels.pressure
        data["molecular pressure"] = labels.molecular_pressure
        data["MBE pressure"] = labels.breakdown["MBE"]["pressure"]
        data["MBE molecular pressure"] = labels.breakdown["MBE"]["molecular pressure"]
    data["MBE energy"] = float(
        units.ev_to_kj_per_mol(labels.breakdown["MBE"]["energy"])
    )
    data["per-body energies"] = {
        str(k): float(units.ev_to_kj_per_mol(v["energy"]))
        for k, v in labels.per_body.items()
    }
    if labels.virial is not None:
        data["per-body pressures"] = {
            str(k): v["pressure"] for k, v in labels.per_body.items()
        }
        data["per-body molecular pressures"] = {
            str(k): v["molecular pressure"] for k, v in labels.per_body.items()
        }
    data["maximum increment net force"] = float(
        units.ev_to_kj_per_mol(labels.max_increment_net_force)
    )
    forces = np.asarray(labels.forces)
    data["maximum force"] = float(units.ev_to_kj_per_mol(np.abs(forces).max()))
    data["rms force"] = float(
        units.ev_to_kj_per_mol(np.sqrt((forces**2).sum(1).mean()))
    )
    return data
