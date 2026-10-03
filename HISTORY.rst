=======
History
=======

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
