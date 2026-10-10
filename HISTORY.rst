=======
History
=======

2026.10.10 -- Ion shells
    * 'Ion shells' *Li+ first shell* makes each Li⁺ and the molecules in its first
      shell one unit of the expansion, rebuilt for every configuration. Around Li⁺
      the per-molecule expansion converges slowly; with the shell as one unit it
      converges as for neutral molecules. 'Shell cutoffs' (default Li-O and Li-F
      2.6 Å) say who joins, a molecule near two ions joins the nearer one, and 'Most
      molecules in a shell' (default 5) caps crowded shells. 'Fragments with a
      shell' *pairs* selects a shell's pairs but not its triples, which hold most
      of the cost. The job output lists the shells per configuration, and the
      increments file lists each unit's molecules. Shells need a contact
      criterion; the dialog offers only those when shells are on.
    * The energy offsets, the charges and the molecular virial stay per molecule.
      Counterpoise is not available with ion shells yet (hidden in the dialog,
      refused at run time).
    * Requires seamm-mbe 2026.10.10.

2026.10.7 -- Energies on the DfE0 scale by default
    * The new default for 'energy offsets', "DfE0 (from the thermochemistry
      database)", puts the labels on the energy-of-formation scale (DfE0, at 0 K)
      that the code steps use: each molecule type's offset comes from the atomic
      reference energies at the high level in SEAMM's thermochemistry database. The
      step prints the offsets it used, and stops before any calculation if the
      database has no atoms at the high level. 'none' and offsets given by hand work
      as before; existing flowcharts keep the value they were saved with.
    * Requires seamm-thermochemistry 2026.9.28.

2026.10.6.3 -- The molecular pressure by body order
    * The summary and the results give the many-body correction's share of the
      molecular pressure, in total and by body order, alongside the atomic one. The
      molecular breakdown is the one to compare with the virial errors of a cluster
      validation; the monomers' share of it is zero.
    * Requires seamm-mbe 2026.10.6.1.

2026.10.6.2 -- A failed level stops the others
    * If a level fails outright (an error, or every one of its calculations failing),
      the other levels are stopped at once instead of being left to run for hours:
      their queued and running calculations are cancelled. A rerun keeps everything
      already finished and submits the rest afresh.
    * Requires seamm-exec 2026.10.6.2.

2026.10.6.1 -- The levels run at the same time
    * The levels that run as batch tasks now run at the same time, with the whole cell
      and the periodic fragments started first. Their long calculations overlap the
      molecular ones instead of following them, roughly halving the time a frame takes.
      Levels that use an MDI engine still run one at a time.

2026.10.6 -- A low level of their own for the triples
    * The new 'triples low level' gives the triples' increments their own molecular low
      level, e.g. r2SCAN-D4 at def2-TZVPPD with the molecular low level at def2-QZVPPD
      for the pairs. Each triple's increment is built entirely at that level, so the
      triples' monomers and pairs are computed at both levels. It is shown only when
      there are triples; by default it is the molecular low level. A triples' low level
      that is the molecular low level under another name is refused.
    * The increments file records each increment's low level.
    * The user guide explains when to use the periodic low level for the monomers only.
    * Requires seamm-mbe 2026.10.6.

2026.10.5 -- A high level of their own for the triples; the increments in a file
    * The new 'triples high level' gives the triples' increments their own high
      level, e.g. revDSD at def2-TZVPPD with the high level at def2-QZVPPD for the
      monomers and pairs. Each triple's increment is built entirely at that level, so
      the triples' monomers and pairs are computed at both levels. It is shown only
      when there are triples; by default it is the high level. A triples' level that
      is the high level under another name is refused, since every monomer and pair
      would be computed twice.
    * Each configuration's increments are written to
      increments_c<configuration id>.json in the step's directory, for later
      analysis: each fragment's molecules, images, pair distances, levels, energy and
      forces, with triples marked "chain" or "closed".
    * Requires seamm-mbe 2026.10.5.

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
