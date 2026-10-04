=======
History
=======

2026.10.4.1 -- Pairwise counterpoise corrects the high level only by default
    * 'pairwise' counterpoise now corrects each pair at the high level only. The new
      choice 'pairwise, both levels' also corrects the molecular low level, as
      'pairwise' did before. In tests on Li⁺/BF₄⁻/water/EC clusters the r2SCAN-D4
      correction was tiny where it could be computed, and its ghost calculations
      sometimes ended at spurious energies, while the revDSD correction was 3–4
      kJ/mol per molecule.
    * The user guide notes that ORCA's ghost-atom gradients are unreliable for
      double hybrids such as revDSD, so counterpoise-corrected forces from them need
      caution.

2026.10.4 -- Requires seamm-exec 2026.10.3.2 for running on a cluster
    * Requires seamm-exec 2026.10.3.2: each bundle of calculations sent to a cluster's
      queue now runs on one node, a calculation runs inside the job itself only when it
      fits there, and bundles get a time limit from their calculations' estimated cost.
      Before, a VASP or ORCA calculation could be spread over several nodes and fail, or
      be cut off by the queue's default time limit.
    * Checked end to end on TinkerCliffs: the 64-water pilot frame reproduces the
      prototype's labels (energy to 0.001 meV, forces to 0.2 meV/Å RMS, pressures to
      about 1 atm).

2026.10.3 -- Initial release: many-body expansion (MBE) corrections
    * A step that estimates high-level energies, forces and, for periodic cells, the
      stress of a cell or a large cluster: a low-level calculation of the whole system
      plus [high - low] corrections computed on monomers, pairs and triples of
      molecules, using the seamm_mbe library.
    * The levels are model chemistries: the high level, a molecular low level, and
      optionally a periodic low level (e.g. VASP) for the fragments of periodic cells,
      run as batches of calculations through seamm_exec.
    * Pair and triple cutoffs, a rule for which triples to include, and a check that
      the cell is large enough for the chosen cutoffs.
    * For a periodic low level, the fragments are placed on the cell's FFT grid so that
      the grid errors cancel between the cell and its fragments.
    * Optional pairwise counterpoise correction of the basis-set superposition error.
    * Optional energy offsets per molecule type, and the results stored as properties
      and written to an extended XYZ file for training machine-learned force fields.
