# -*- coding: utf-8 -*-

"""Energy offsets onto the DfE0 scale from SEAMM's thermochemistry database.

A configuration's labels are on the high level's absolute energy scale (the
VASP cell's pseudopotential reference is replaced by the high-level monomers').
The training sets use the energy of formation at 0 K instead, DfE0 = E -
sum(atomic reference energies) + sum(atomic heats of formation at 0 K), with the
atomic references at the same level. For a molecule the difference DfE0 - E is a
constant, so it enters as an offset per molecule type, exactly the number the
code steps' own DfE0 uses (seamm_thermochemistry.formation_energy, anchored at
0 K).
"""

from collections import Counter

from seamm_mbe.units import EV_TO_KJ_PER_MOL

#: The parameter value that asks for the offsets from the database
AUTOMATIC = "DfE0 (from the thermochemistry database)"


def is_automatic(text):
    """Whether the 'energy offsets' parameter asks for DfE0 offsets."""
    return (text or "").strip().lower() in (AUTOMATIC.lower(), "dfe0")


def reference(model_chemistry, provider=None):
    """(code, method, settings) under which the database holds the atomic
    references for a model chemistry. A program can say with a
    ``thermochemistry_reference(model_chemistry)`` hook on its step; otherwise
    the code is the step's name in lower case, the method the model
    chemistry's (in its model-chemistry spelling, as the reference runs were
    tagged) and the settings its basis, without a ``bse:`` prefix."""
    hook = getattr(provider, "thermochemistry_reference", None)
    if hook is not None:
        return hook(model_chemistry)
    code = str(model_chemistry.get("step") or model_chemistry.get("owner")).lower()
    method = model_chemistry.get("method") or ""
    basis = model_chemistry.get("basis") or ""
    if basis.lower().startswith("bse:"):
        basis = basis[4:]
    return code, method, basis


def compositions(systems):
    """{type name: Counter of element symbols} for the molecule types of
    seamm_mbe Systems."""
    result = {}
    for system in systems:
        for molecule in system.molecules:
            if molecule.type not in result:
                result[molecule.type] = Counter(
                    system.symbols[i] for i in molecule.atoms
                )
    return result


def dfe0_offsets(model_chemistry, systems, provider=None):
    """The offsets DfE0 - E in eV per molecule, by type, at a model chemistry.

    Raises
    ------
    ValueError
        If seamm_thermochemistry is not installed or has no atomic references
        for the elements at this level.
    """
    try:
        from seamm_thermochemistry import ThermoDB, formation_energy
    except ImportError:
        raise ValueError(
            "The energy offsets onto the DfE0 scale need seamm_thermochemistry, "
            "which is not installed. Set 'energy offsets' to 'none', or give them."
        ) from None
    code, method, settings = reference(model_chemistry, provider)
    level = model_chemistry.get("level", f"{code}:{method}/{settings}")
    offsets = {}
    try:
        database = ThermoDB(read_only=True)
    except (FileNotFoundError, OSError) as e:
        raise ValueError(
            "The energy offsets onto the DfE0 scale need the thermochemistry "
            f"database, which is not installed here ({e}). Install it with "
            "seamm-thermochemistry's installer, or set 'energy offsets' to 'none' "
            "or give them."
        ) from None
    with database as db:
        for name, composition in sorted(compositions(systems).items()):
            try:
                kj = formation_energy(
                    dict(composition),
                    0.0,
                    db,
                    code,
                    method,
                    settings=settings,
                    anchor=True,
                    anchor_at_0K=True,
                    units="kJ/mol",
                )
            except Exception as e:
                raise ValueError(
                    f"No DfE0 energy offset for {name} at {level}: the "
                    f"thermochemistry database has no atomic references for it "
                    f"({type(e).__name__}: {e}). Set 'energy offsets' to 'none' "
                    "for absolute energies, or give the offsets as 'type value'."
                ) from None
            offsets[name] = kj / EV_TO_KJ_PER_MOL
    return offsets
