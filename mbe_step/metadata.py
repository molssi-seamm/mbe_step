# -*- coding: utf-8 -*-

"""This file contains metadata describing the results from Mbe"""

metadata = {}

"""Description of the computational models for Mbe.

Hamiltonians, approximations, and basis set or parameterizations,
only if appropriate for this code. For example::

    metadata["computational models"] = {
        "Hartree-Fock": {
            "models": {
                "PM7": {
                    "parameterizations": {
                        "PM7": {
                            "elements": "1-60,62-83",
                            "periodic": True,
                            "reactions": True,
                            "optimization": True,
                            "code": "mopac",
                        },
                        "PM7-TS": {
                            "elements": "1-60,62-83",
                            "periodic": True,
                            "reactions": True,
                            "optimization": False,
                            "code": "mopac",
                        },
                    },
                },
            },
        },
    }
"""
# metadata["computational models"] = {
# }

"""Description of the Mbe keywords.

(Only needed if this code uses keywords)

Fields
------
description : str
    A human readable description of the keyword.
takes values : int (optional)
    Number of values the keyword takes. If missing the keyword takes no values.
default : str (optional)
    The default value(s) if the keyword takes values.
format : str (optional)
    How the keyword is formatted in the MOPAC input.

For example::
    metadata["keywords"] = {
        "0SCF": {
            "description": "Read in data, then stop",
        },
        "ALT_A": {
            "description": "In PDB files with alternative atoms, select atoms A",
            "takes values": 1,
            "default": "A",
            "format": "{}={}",
        },
    }
"""
# metadata["keywords"] = {
# }

"""Properties that Mbe produces.
`metadata["results"]` describes the results that this step can produce. It is a
dictionary where the keys are the internal names of the results within this step, and
the values are a dictionary describing the result. For example::

    metadata["results"] = {
        "total_energy": {
            "calculation": [
                "energy",
                "optimization",
            ],
            "description": "The total energy",
            "dimensionality": "scalar",
            "methods": [
                "ccsd",
                "ccsd(t)",
                "dft",
                "hf",
            ],
            "property": "total energy#Psi4#{model}",
            "type": "float",
            "units": "E_h",
        },
    }

Fields
______

calculation : [str]
    Optional metadata describing what subtype of the step produces this result.
    The subtypes are completely arbitrary, but often they are types of calculations
    which is why this is name `calculation`. To use this, the step or a substep
    define `self._calculation` as a value. That value is used to select only the
    results with that value in this field.

description : str
    A human-readable description of the result.

dimensionality : str
    The dimensions of the data. The value can be "scalar" or an array definition
    of the form "[dim1, dim2,...]". Symmetric tringular matrices are denoted
    "triangular[n,n]". The dimensions can be integers, other scalar
    results, or standard parameters such as `n_atoms`. For example, '[3]',
    [3, n_atoms], or "triangular[n_aos, n_aos]".

methods : str
    Optional metadata like the `calculation` data. `methods` provides a second
    level of filtering, often used for the Hamiltionian for *ab initio* calculations
    where some properties may or may not be calculated depending on the type of
    theory.

property : str
    An optional definition of the property for storing this result. Must be one of
    the standard properties defined either in SEAMM or in this steps property
    metadata in `data/properties.csv`.

type : str
    The type of the data: string, integer, or float.

units : str
    Optional units for the result. If present, the value should be in these units.
"""


def _result(description, units="", dimensionality="scalar", _type="float", **kw):
    return {
        "description": description,
        "dimensionality": dimensionality,
        "type": _type,
        "units": units,
        **kw,
    }


metadata["results"] = {
    "energy": _result(
        "The MBE-corrected total energy",
        "kJ/mol",
        property="energy#MBE#{model}",
        format=".3f",
    ),
    "reference energy": _result(
        "The energy on the reference scale of the labels (with the offsets)",
        "kJ/mol",
        format=".3f",
    ),
    "gradients": _result(
        "The MBE-corrected gradients",
        "kJ/mol/Å",
        "[n_atoms, 3]",
        "json",
        property="gradients#MBE#{model}",
    ),
    "stress": _result(
        "The MBE-corrected stress, Voigt xx yy zz yz xz xy, sigma = -P",
        "GPa",
        "[6]",
        "json",
        property="stress#MBE#{model}",
        format=".4f",
    ),
    "pressure": _result("The atomic configurational pressure", "atm", format=".1f"),
    "molecular pressure": _result(
        "The molecular configurational pressure", "atm", format=".1f"
    ),
    "MBE pressure": _result(
        "The many-body correction's contribution to the pressure", "atm", format=".1f"
    ),
    "MBE energy": _result("The many-body correction to the energy", "kJ/mol"),
    "per-body energies": _result(
        "The correction to the energy by order, in kJ/mol", "", "json", "json"
    ),
    "per-body pressures": _result(
        "The correction to the pressure by order, in atm", "", "json", "json"
    ),
    "fragment counts": _result(
        "The selected and auxiliary fragments by order", "", "json", "json"
    ),
    "maximum increment net force": _result(
        "The largest component of any increment's net force (ideally zero)",
        "kJ/mol/Å",
    ),
    "counterpoise correction": _result(
        "The pairwise counterpoise correction to the energy", "kJ/mol", format=".3f"
    ),
    "counterpoise fallbacks": _result(
        "Pairs whose counterpoise gradient was unphysical (forces left uncorrected)",
        "",
        _type="integer",
    ),
    "maximum force": _result(
        "The largest Cartesian component of the forces", "kJ/mol/Å"
    ),
    "rms force": _result("The RMS force on the atoms", "kJ/mol/Å"),
    "configuration name": _result("The configuration", "", _type="string"),
    "model chemistry": _result("The high level", "", _type="string"),
}
