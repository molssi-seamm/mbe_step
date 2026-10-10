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

The step uses up to six Model Chemistries:

**High level**
    The level the labels approximate, run on every fragment. By default it is
    the one chosen by the Model Chemistry step before this step.
**Triples' high level**
    Optional: a high level for the triples' increments when it differs from the
    pairs', e.g. revDSD at def2-QZVPPD for the monomers and pairs and at
    def2-TZVPPD for the triples. Each triple's increment is built entirely at its
    own level (the triple, its three pairs and its three monomers), so those
    monomers and pairs are computed at both levels. Shown only when there are
    triples; by default it is the high level.
**Molecular low level**
    The cheap level run on the fragments referenced to a molecular code. It
    is required.
**Triples' low level**
    Optional: a molecular low level for the triples' increments when it differs
    from the pairs', e.g. r2SCAN-D4 at def2-QZVPPD for the pairs and at
    def2-TZVPPD for the triples. As with the triples' high level, each triple's
    increment is built entirely at its own levels, so its monomers and pairs are
    computed at both. Shown only when there are triples; by default it is the
    molecular low level.

    Different levels for different orders are not a shortcut: an order's
    increment is a difference within its own ladder (a triple minus its pairs
    plus its monomers, all at the same levels), so its levels only need to be
    consistent within that ladder, never across orders.
**Periodic low level**
    Optional. A periodic code for the compact fragments: monomers, and pairs
    closer than a distance r1. It should be the code of the cell's low level,
    since that cancels the cell's own low-level error best. Extended
    fragments use the molecular low level, because they pick up interactions
    with their periodic images in affordable boxes. This is the "mixed" scheme.

    Pairs are more sensitive to their images than this suggests. In boxes the
    size of a 32-molecule cell of ethylene carbonate (15.2 Å), VASP pair
    increments differed from ORCA's by 15 meV/Å rms, and by 6.6 meV/Å even for
    pairs 8 Å or more from their images, while the monomers were converged to
    1 meV/Å. Extrapolated to an infinite box, VASP and ORCA agreed to about
    1 meV/Å. So for polar molecules in small cells, use the periodic low level
    for the monomers only: set **periodic pair cutoff** to 0.
**Whole-system low level**
    The cheap level run on the whole cell or cluster. *automatic* uses the
    periodic low level for a cell and the molecular low level for a cluster.
    For a cell with a periodic low level it must be that same level, so that the
    compact fragments and the cell are the same calculation.

Each field accepts any model chemistry the installed programs offer, typed or
chosen, and ``$variables`` in its parts. A periodic level must come from a program
that declares the sign convention of its stress: MOPAC does from mopac_step
2026.10.3.2, VASP from vasp_step 2026.10.3. A periodic low level applies only to
periodic cells; for a cluster set it to *none*.

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
  order (energy, and both the atomic and the molecular pressure), the fragment
  counts and the largest increment net force, as results for variables and
  tables. Each order's share of the molecular pressure is its virial less its
  own intramolecular part, so the shares add up and the monomers' share is zero:
  it is the breakdown to compare with the virial errors of a cluster validation;
- an extended XYZ file (default ``mbe_labels.extxyz`` in the step's directory;
  ``job:NAME`` puts it in the job's directory) with ``REF_energy``,
  ``REF_forces`` and, for a cell, the nine-value ``REF_stress``. Rerunning
  replaces a configuration's frame rather than adding a second copy.

Each selected fragment's increment is also written to
``increments_c<configuration id>.json`` in the step's directory, for later
analysis such as a correction by triple topology. Each record holds the
fragment's molecules, periodic images and pair distances; its low and high
levels; and its energy (eV) and forces (eV/Å). Triples are marked "chain" or
"closed" by how many of their pairs are within the triple cutoff.

**Energy offsets** (eV per molecule, by type) put the labels on the scale of
other training data. The labels come out on the high level's absolute scale: the
monomers' high-level energies replace the periodic code's pseudopotential
reference in the cell's energy. An offset per molecule type then shifts them;
forces and stress are unchanged.

* **DfE0 (from the thermochemistry database)**, the default, gives the energy of
  formation at 0 K, the scale of the code steps' own ``DfE0``. Each type's offset
  is DfE0 − E for one molecule, from the atomic reference energies at the high
  level (the monomers' level) and the atoms' heats of formation at 0 K. The step
  prints the offsets it used, e.g. ``EC 9300.36250`` eV per molecule at
  revDSD-PBEP86-D4/def2-QZVPPD. If the database has no atoms at that level the
  step stops before any calculation.
* **none** gives absolute energies.
* Or give the offsets as ``type value`` entries separated by ``;``, e.g.
  ``water 2074.69325``. Every molecule type present then needs one.

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
archived, keeping a configuration to a few files.

The levels that run as batch tasks run at the same time, and the long periodic
calculations (the whole cell and the periodic fragments) are started first, so
they overlap the many molecular calculations instead of following them. For the
32-molecule ethylene carbonate pilot this roughly halves the time to the labels.
Levels that use an MDI engine run one at a time.

If a level fails outright (an error, or every one of its calculations failing),
no configuration can be labelled, so the other levels are stopped at once rather
than left to run for hours: their calculations in the queue are cancelled. A
rerun of the job keeps every result already finished and submits the rest
afresh.

Periodic fragments on the cell's grid
=====================================

A plane-wave code's energy depends slightly on where each atom sits relative to
its FFT grid. The increments subtract fragments from each other and from the
cell, so this error cancels only if every atom keeps the same position relative to
the grid. With a periodic low level the step therefore gives the cell an explicit
grid, with a spacing no larger than **Largest FFT grid spacing** (default
0.0829 Å), and puts each periodic fragment in a box that is a whole number of
those grid steps, at least 12 Å and the fragment's extent plus **Box padding**
(default 7.5 Å), shifted by whole grid steps. The cell must be orthorhombic.
These settings appear only when a periodic low level is chosen.

The step runs the cell at the Γ point alone, so with VASP the cell must be at least
10 Å across; VASP refuses a smaller cell, and the step reports it as a failed
calculation.

Ion shells
==========

Around a small, strongly polarizing ion such as Li⁺, the per-molecule expansion
converges slowly. In a test cluster of Li⁺, BF₄⁻, six waters and four EC molecules,
the Li⁺ 3-body [revDSD − r2SCAN] increments summed to about +20 kJ/mol and the
4-body ones to −5 to −18 kJ/mol.

**Ion shells** *Li+ first shell* fixes this. Each Li⁺ and the molecules in its
first shell become one unit of the expansion, and pairs and triples are then
built from units. In the test cluster this brought the error at triples to +1–2
kJ/mol and the forces to 3.6 meV/Å RMS.

- **Who joins a shell:** a molecule with an atom of a listed element within that
  element's cutoff of the ion. **Shell cutoffs** defaults to ``Li O 2.6; Li F
  2.6``, just inside the first minima of g(r) in carbonate electrolytes.
- **When:** the shells are rebuilt for every configuration.
- **Anions:** an anion in contact with Li⁺ joins its shell. Anions get no shells
  of their own.
- **Shared molecules:** a molecule within reach of two ions joins the nearer one,
  so shells never overlap.
- **Crowded shells:** **Most molecules in a shell** (default 5) keeps a crowded
  shell's nearest molecules. The others stay units of their own.
- **Cheaper truncation:** a shell has 40–70 atoms, so its triples hold most of the
  cost. **Fragments with a shell** *pairs* selects a shell's pairs but none of its
  triples.

What stays per molecule:

- the energy offsets;
- the charges, which a shell's charge sums;
- the molecular virial.

The job output lists the shells and their sizes for each configuration.
In ``increments_c<id>.json``, each increment's ``members`` lists the real
molecules in each of its units.

**A contact criterion.** Ion shells need **Distance between molecules** set to a
contact criterion. A shell's designated atom is its ion, so measuring from it
would miss partners near the shell's members. The dialog offers only the contact
criteria when shells are on.

**Larger cells.** The selection needs the cell's smallest width to exceed the pair
cutoff plus twice the largest unit's radius, so shells need larger cells than
single molecules (about 17 Å for carbonate shells). A cell that is too small is
refused, with the reason.

**Not yet available:**

- **counterpoise:** it is hidden in the dialog, and refused at run time;
- **a periodic low level for a charged shell's monomer term:** a charged unit in a
  box picks up a Madelung energy that depends on the box. Validate it before
  relying on it.

Counterpoise
============

**Counterpoise** *pairwise* corrects each selected pair for the basis-set
superposition error (Boys–Bernardi) at the high level. *pairwise, both levels*
also corrects the molecular low level, for the pairs referenced to it.

- **Cost:** two calculations per pair at each corrected level, each monomer in the
  pair's basis with the other monomer's atoms as ghosts.
- **Where it enters:** the 2-body sum only. Triples still subtract the uncorrected
  pairs, so the pairs' error does not move into the 3-body terms.
- **Unphysical ghost gradients:** where they occur, the energy is still corrected
  but the forces are left uncorrected, and the step counts these pairs. The total
  correction and that count are stored as results.

The high level alone is usually right. In tests on Li⁺/BF₄⁻/water/EC clusters in
def2-TZVPPD:

- at revDSD the corrections were 3–4 kJ/mol per molecule;
- at r2SCAN-D4 they were tiny (about −0.02 kJ/mol per pair), and about 5% of its
  ghost calculations ended at spurious energies.

**Forces with a double-hybrid high level:** ORCA's gradients on ghost atoms are
unreliable for double hybrids such as revDSD-PBEP86-D4 (orca_step#42). The step
then keeps the uncorrected forces only for pairs whose corrected gradient fails the
net-force check, so treat counterpoise-corrected forces from a double hybrid with
caution.

Counterpoise needs a level that runs as batch calculations with ghost atoms, and
is not available for periodic cells (it is hidden in the dialog when a periodic
low level is chosen).

4-body terms come in a later version.
