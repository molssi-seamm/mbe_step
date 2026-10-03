.. _user-guide:

**********
User Guide
**********

The MBE step estimates high-level energies, forces and stress for a periodic cell
or a large cluster that the high-level method can't treat directly. It corrects a
cheap calculation of the whole system with many-body increments, each the
difference [high − low] computed on a small isolated fragment::

    E(system) ≈ E_low(system) + Σ_F dE_F

where dE_F is the F-body part of [high − low] for each selected fragment F:
monomers, pairs and triples. Forces and the virial are assembled the same way. The
bookkeeping is done by the `seamm_mbe <https://molssi-seamm.github.io/seamm_mbe/>`_
library.

The levels of theory
====================

The step uses up to four Model Chemistries:

**High level**
    The level the labels approximate, run on every fragment. By default it is
    the one chosen by the Model Chemistry step before this step.
**Molecular low level**
    The cheap level run on the fragments referenced to a molecular code. It
    is required.
**Periodic low level**
    Optional. A periodic code for the compact fragments: monomers, and pairs
    closer than a distance r1. It should be the code of the cell's low level,
    since that cancels the cell's own low-level error best. Extended
    fragments use the molecular low level, because they pick up interactions
    with their periodic images in affordable boxes. This is the "mixed" scheme.
**Whole-system low level**
    The cheap level run on the whole cell or cluster. *automatic* uses the
    periodic low level for a cell and the molecular low level for a cluster.

Each field accepts any model chemistry the installed programs offer, typed or
chosen, and ``$variables`` in its parts. A periodic level must come from a program
that declares the sign convention of its stress. MOPAC does so from mopac_step
2026.10.3.2.

Each increment is built entirely from calculations at its own low level. A triple
referenced to the molecular level subtracts the molecular-level increments of its
pairs, even where those pairs are referenced to the periodic level in the sum.

Choosing the fragments
======================

**Maximum order**
    1 (monomers), 2 (pairs) or 3 (triples).
**Distance between molecules**
    Measured between the molecules' designated atoms (water's O, a
    carbonate's carbonyl C, an ion), between their closest atoms or closest
    heavy atoms, or between their centres of mass or geometry.
**Cutoffs**
    One value, or a table by molecule type for mixtures, e.g. ``water water
    4.5; Li+ * 3.0; * * 5.0``. ``*`` matches any type, and the most specific
    entry wins. Every pair of types present must be covered.
**Triples**
    *connected* selects a triple when at least two of its three pairs are
    within the triple cutoff (a hub bonded to both others). *compact*
    requires all three. A connected triple's third pair is computed as well,
    since its increment needs it, even when it is longer than the pair cutoff.
    It does not enter the 2-body sum.

The molecules are found from the bonds (they are perceived if the structure has
none) and typed by formula and topology. Charges come from the structure's formal
charges, or from the type catalog: Li⁺, BF₄⁻, PF₆⁻ and the common monatomic ions.
The whole system runs with the charge its molecules add up to. If the
configuration's own charge is set and disagrees, the step refuses to run.
A configuration charge of 0 counts as not set (plain xyz input has no charge),
so a structure holding a Li⁺ runs at +1 even if its charge field says 0. To make
a molecule neutral, give the structure formal charges or give the charge of its
type.

A cell too small for the cutoffs is refused, because the same molecules could
form two different fragments.

The labels
==========

For each configuration the step stores:

- the energy (kJ/mol), the gradients (kJ/mol/Å) and, for a cell, the stress (a
  Voigt [6] vector in GPa, σ = −P), as properties ``energy#MBE#<model>`` and so
  on, with the gradients also on the atoms;
- the atomic and molecular pressures, the correction and its breakdown by body
  order, the fragment counts and the largest increment net force, as results
  for variables and tables;
- an extended XYZ file (default ``mbe_labels.extxyz`` in the step's directory;
  ``job:NAME`` puts it in the job's directory) with ``REF_energy``,
  ``REF_forces`` and, for a cell, the nine-value ``REF_stress``. Rerunning
  replaces a configuration's frame rather than adding a second copy.

**Energy offsets** (eV per molecule, by type) put the labels on the scale of
other training data. The default, ``water 2074.69325``, is the water training
sets' formation-energy scale. Use ``none`` for absolute energies.

A configuration with a failed or missing calculation gets no labels. The step
reports which calculations failed, stores the complete configurations, and
stops with an error. Rerunning the job reuses every finished calculation and
retries the failed ones.

Running the calculations
========================

Each level runs through SEAMM's task layer: as batch tasks wherever the job's
target sends them (a local pool, or bundled jobs on a cluster), or on a warm MDI
engine for cheap codes run locally. The *Execution* settings give the cores and
memory of each calculation and how many share a batch job. Finished bundles are
archived, keeping a configuration to a few files. In this version the levels run
one after the other.

Counterpoise corrections, the registered-box VASP fragments, and 4-body terms
come in later versions.
