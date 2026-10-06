2026-10-05: preparing the 32-EC pilot
=====================================

The ethylene carbonate (EC) reference, as decided on 2026-10-05: QZ pairs and TZ
triples, no counterpoise, with criterion "contact", a pair cutoff of 4.5 Å and
connected triples within 3.0 Å. The low level is mixed: VASP for the monomers and
close pairs, ORCA for the rest, and the cell at the VASP level.

What it needed
--------------

- **A high level per order:** seamm-mbe 2026.10.5 (``high_by_order``) and mbe_step
  2026.10.5 (the "triples high level" and the ``increments_c<id>.json`` file).
- **Fragments through different images:** seamm-mbe 2026.10.5.1. The 32-EC cell is
  narrower than the uniqueness bound, which is 20.07 Å for these pairs and
  25.61 Å for these triples, so the same molecules can form several triples.

  - Each of these triples is enumerated. The minimum-image one keeps the plain name
    and the others carry their images.
  - A molecule bonded to two images of a partner only skips the placement that
    would hold the partner twice.
  - A molecule within the pair cutoff of its own image is refused.

Expected counts for the pilot frame
-----------------------------------

The frame is from ChemAI job 5170: 32 EC, a cubic cell of 15.2423 Å, with the
cell held fixed in the minimization. Enumerated with seamm-mbe 2026.10.5.1 at the
settings above, it gives:

- 2 warnings, one for the pairs' bound and one for the triples';
- 32 monomers;
- 241 pairs, 15.1 contacts per molecule;
- 878 connected triples. Fifteen molecule sets each form three different triples,
  e.g. (0, 1, 25) gives ``t00_01_25``, ``t00_01_25_xm0mm0m`` and
  ``t00_01_25_xm00m00``;
- 443 auxiliary pairs, the triples' outer pairs and image pairs;
- 1,594 fragments in all, with unique names.

A pilot run on this frame should list exactly these. A different count means the
enumeration or the settings changed.

The staged pre-pilot (TinkerCliffs job 7859415)
-----------------------------------------------

The pre-pilot ran before the image fix, so its 48 fragments were sampled from the
2×2×2 supercell, where images are not ambiguous. Each fragment was mapped back
to the 32-EC cell's grid. The run computes:

- all 48 fragments with VASP in a box the size of the cell;
- all 48 with ORCA at TZ, and the 43 monomers and pairs at QZ;
- the tightest stacked pair, d00_26, and its two monomers with VASP in a 1.5×
  box.

It compares the VASP and ORCA increments by body order and contact class, and the
pair's box-image effect at 1.5×.

The probe calibrated the VASP cost:

- an EC monomer in the cell-sized box took 1,543 s on 8 ranks, 3.4 core-h, with a
  peak of about 2.3 GB per rank, so 4 GB per rank is needed;
- that is 0.656 of the model's time.
