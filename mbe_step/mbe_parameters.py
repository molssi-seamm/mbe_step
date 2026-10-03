# -*- coding: utf-8 -*-
"""
Control parameters for the MBE step in a SEAMM flowchart
"""

import logging
import pprint  # noqa: F401

import seamm

logger = logging.getLogger(__name__)

#: The distance criteria, as shown, and their names in seamm_mbe
CRITERIA = {
    "designated atoms": "designated",
    "closest contact": "contact",
    "closest heavy-atom contact": "heavy contact",
    "center of mass": "com",
    "center of geometry": "cog",
}


class MbeParameters(seamm.Parameters):
    """
    The control parameters for the MBE step.

    The step corrects a cheap calculation of a whole periodic cell (or cluster)
    with many-body increments of [high - low] computed on small isolated
    fragments: monomers, selected pairs and triples. The bookkeeping is in the
    ``seamm_mbe`` library; this step runs the calculations through each level's
    Model Chemistry.

    See Also
    --------
    Mbe, TkMbe, MbeStep
    """

    parameters = {
        # Which structures: SEAMM's standard block
        **seamm.standard_parameters.structure_selection_parameters,
        # ------------------------------------------------------------- levels
        "high level": {
            "default": "current model chemistry",
            "kind": "string",
            "default_units": "",
            "enumeration": ("current model chemistry",),
            "format_string": "",
            "description": "High level:",
            "help_text": (
                "The model chemistry the labels approximate, run on every "
                "fragment, e.g. 'ORCA:DFT@revDSD-PBEP86-D4/2021/def2-TZVPPD'. By "
                "default the one chosen by the Model Chemistry step before this "
                "step."
            ),
        },
        "molecular low level": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Molecular low level:",
            "help_text": (
                "The cheap model chemistry run on the fragments whose increments "
                "are referenced to a molecular code, e.g. "
                "'ORCA:DFT@r2SCAN-D4/def2-TZVPPD'."
            ),
        },
        "periodic low level": {
            "default": "none",
            "kind": "string",
            "default_units": "",
            "enumeration": ("none",),
            "format_string": "",
            "description": "Periodic low level:",
            "help_text": (
                "The periodic code used for the compact fragments (monomers and "
                "close pairs), the same code as the cell's low level, which "
                "cancels its errors best. 'none' references every increment to "
                "the molecular low level. For now the fragments are run as "
                "isolated molecules, so any model chemistry is accepted here; "
                "the registered periodic boxes for VASP come in a later "
                "version."
            ),
        },
        "cell low level": {
            "default": "automatic",
            "kind": "string",
            "default_units": "",
            "enumeration": ("automatic",),
            "format_string": "",
            "description": "Whole-system low level:",
            "help_text": (
                "The cheap model chemistry run on the whole cell (or cluster). "
                "'automatic' uses the periodic low level for a periodic system "
                "and the molecular low level for a cluster."
            ),
        },
        # ---------------------------------------------------------- fragments
        "maximum order": {
            "default": "3",
            "kind": "integer",
            "default_units": "",
            "enumeration": ("1", "2", "3"),
            "format_string": "",
            "description": "Highest order:",
            "help_text": (
                "The largest fragments: 1 = monomers, 2 = pairs, 3 = triples."
            ),
        },
        "distance criterion": {
            "default": "designated atoms",
            "kind": "enum",
            "default_units": "",
            "enumeration": tuple(CRITERIA),
            "format_string": "",
            "description": "Distance between molecules:",
            "help_text": (
                "How the distance between two molecules is measured: between "
                "their designated atoms (water's O, a carbonate's carbonyl C, an "
                "ion), the closest contact of any atoms or of the heavy atoms, or "
                "the centres of mass or geometry."
            ),
        },
        "cutoffs": {
            "default": "single value",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("single value", "table by type pair"),
            "format_string": "",
            "description": "Cutoffs:",
            "help_text": (
                "One cutoff for every pair of molecules, or a table by the types "
                "of the two molecules (needed for mixtures)."
            ),
        },
        "pair cutoff": {
            "default": 4.5,
            "kind": "float",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Pair cutoff:",
            "help_text": "Pairs of molecules closer than this are selected.",
        },
        "pair cutoff table": {
            "default": "* * 4.5",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Pair cutoffs (Å):",
            "help_text": (
                "The pair cutoff in Å by molecule types, as entries "
                "'type type cutoff' separated by ';', e.g. 'water water 4.5; Li+ "
                "* 3.0; * * 5.0'. '*' matches any type and the most specific "
                "entry wins. Every pair of types present must be covered."
            ),
        },
        "triple rule": {
            "default": "connected",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("none", "connected", "compact"),
            "format_string": "",
            "description": "Triples:",
            "help_text": (
                "Which triples are selected: 'connected' if at least two of the "
                "three pairs are within the triple cutoff (a hub bonded to both "
                "others), 'compact' if all three are."
            ),
        },
        "triple cutoff": {
            "default": 3.5,
            "kind": "float",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Triple cutoff:",
            "help_text": (
                "The distance within which two molecules of a triple count as "
                "bonded, for the triple rule."
            ),
        },
        "triple cutoff table": {
            "default": "* * 3.5",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Triple cutoffs (Å):",
            "help_text": (
                "The triple cutoff in Å by molecule types, in the same form as "
                "the pair cutoffs."
            ),
        },
        # ------------------------------------------------- low-level assignment
        "periodic monomers": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Periodic low level for monomers:",
            "help_text": (
                "Reference the monomer increments to the periodic low level."
            ),
        },
        "periodic pair cutoff": {
            "default": 3.5,
            "kind": "float",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Periodic low level for pairs closer than:",
            "help_text": (
                "Reference the increments of pairs closer than this to the "
                "periodic low level; extended fragments pick up interactions "
                "with their images in affordable periodic boxes, so they use the "
                "molecular low level. 0 for none."
            ),
        },
        # ------------------------------------------------------ energy scale
        "energy offsets": {
            "default": "none",
            "kind": "string",
            "default_units": "",
            "enumeration": ("none",),
            "format_string": "",
            "description": "Energy offsets (eV per molecule):",
            "help_text": (
                "Added to the energy per molecule of each type, to put the labels "
                "on the scale of other training data, as entries 'type offset' "
                "separated by ';', e.g. 'water 2074.69325' (the water training "
                "sets' formation-energy scale). Every molecule type present needs "
                "one. 'none' keeps the absolute energy."
            ),
        },
        # ------------------------------------------------------------ output
        "extxyz file": {
            "default": "mbe_labels.extxyz",
            "kind": "string",
            "default_units": "",
            "enumeration": ("none",),
            "format_string": "",
            "description": "Labels file:",
            "help_text": (
                "The extended XYZ file the labels are appended to (energy, forces "
                "and, for cells, the stress), in the step's directory; "
                "'job:NAME' puts it in the job's directory, so a loop gathers "
                "every configuration into one file. 'none' writes no file."
            ),
        },
        # --------------------------------------------------------- execution
        "molecular ranks": {
            "default": 4,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Cores per molecular calculation:",
            "help_text": "MPI ranks for each molecular fragment calculation.",
        },
        "molecular memory": {
            "default": 1500,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Memory per core (MB):",
            "help_text": "Memory per rank of a molecular calculation, in MB.",
        },
        "molecular bundle": {
            "default": 240,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Molecular calculations per job:",
            "help_text": (
                "How many molecular fragment calculations share one batch job on "
                "a queue."
            ),
        },
        "periodic ranks": {
            "default": 8,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Cores per periodic calculation:",
            "help_text": "MPI ranks for each periodic fragment calculation.",
        },
        "periodic memory": {
            "default": 2000,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Memory per core (MB):",
            "help_text": "Memory per rank of a periodic calculation, in MB.",
        },
        "periodic bundle": {
            "default": 8,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Periodic calculations per job:",
            "help_text": (
                "How many periodic fragment calculations share one batch job on "
                "a queue."
            ),
        },
        "cell ranks": {
            "default": 16,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Cores for the whole system:",
            "help_text": "MPI ranks for the low-level calculation of the cell.",
        },
        "cell memory": {
            "default": 2000,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Memory per core (MB):",
            "help_text": "Memory per rank of the cell's calculation, in MB.",
        },
        "bundle walltime": {
            "default": 4.0,
            "kind": "float",
            "default_units": "h",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Time limit per job:",
            "help_text": "The walltime of each batch job on a queue.",
        },
        "archive": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Archive the calculations:",
            "help_text": (
                "Pack each finished batch of fragment calculations into a tar "
                "file, keeping a frame to a few files instead of thousands."
            ),
        },
        # ------------------------------------------------------------ results
        "results": {
            "default": {},
            "kind": "dictionary",
            "default_units": None,
            "enumeration": tuple(),
            "format_string": "",
            "description": "results",
            "help_text": "The results to save to variables or in tables.",
        },
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters
            above are used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units
            for updating the current, default values.

        Returns
        -------
        None
        """

        logger.debug("MbeParameters.__init__")

        super().__init__(defaults={**MbeParameters.parameters, **defaults}, data=data)
