2026-10-08: the 32-FEC pilot
============================

The fluoroethylene carbonate (FEC) pilot uses option A, as decided for EC on
2026-10-05:

- **Monomers:** VASP r2SCAN-D4 (hard PAW, 1,200 eV) in registered, cell-sized
  boxes.
- **Pairs:** ORCA r2SCAN-D4/def2-QZVPPD as the low level and
  revDSD-PBEP86-D4(2021)/def2-QZVPPD as the high level.
- **Triples:** the same two methods at def2-TZVPPD, with each triple's pairs and
  monomers recomputed at those levels.
- **Selection:** criterion "contact", a pair cutoff of 4.5 Å, connected triples
  within 3.0 Å, no counterpoise.
- **The cell:** VASP r2SCAN-D4.
- **Energy offsets:** automatic, DfE0 from the thermochemistry database (mbe_step
  2026.10.7), FEC 11997.11204 eV per molecule at def2-QZVPPD.

It finished on TinkerCliffs as job 7903005 (2026-10-08, return code 0) with
mbe-step 2026.10.7, seamm-mbe 2026.10.6.1, seamm-exec 2026.10.7.1 and vasp-step
2026.10.8. It took six attempts; see `The runs`_.

The files
---------

All are under ``/projects/seamm/psaxe/fec_pilot``:

- ``fec32.extxyz``: the frame, 32 FEC in a cubic cell of 15.70789 Å, with the
  velocities stripped;
- ``fec_pilot.flow``: the flowchart, with cell ranks 64 (``fec_pilot_cell32.flow``
  is the 32-rank version);
- ``run/3/mbe_labels.extxyz``: the labels (energy, forces, stress), already on
  the DfE0 scale: REF_energy = −286.50801 eV;
- ``run/3/increments_c1.json``: every increment's energy and forces;
- ``scf_tests/``: the SCF tests of the cell (below).

Results
-------

========================================  ==============  ==============
                                          32 FEC          32 EC
========================================  ==============  ==============
monomers                                  32              32
pairs (+ outer and image pairs)           230 (+408)      241 (+443)
connected triples                         764             878
— closed / chain                          86 / 678
2-body, kJ/mol per molecule               −3.043          −4.238
3-body, kJ/mol per molecule               +0.356          +1.189
P, atomic virial (atm)                    +4,824.3        +13,612.4
P, molecular virial (atm)                 −6,180.1        −6,800.3
MBE correction to P, atomic (atm)         −3,633.7        −2,059.8
— 1-body / 2-body / 3-body                −4,011.7 /      −2,606.7 /
                                          +72.9 /         −33.3 /
                                          +305.1          +580.2
MBE correction to P, molecular (atm)      −77.2           −71.5
— 1-body / 2-body / 3-body                −0.6 /          −1.0 /
                                          −357.1 /        −756.2 /
                                          +280.5          +685.7
largest increment net force (meV/Å)       9.15            7.25
cost (core-h)                             ~4,600          5,180
========================================  ==============  ==============

Notes:

- A triple is "closed" when all three of its pairs are within 3.0 Å, and a
  "chain" when two are.
- The EC molecular breakdown was computed afterwards from its increments
  (``ec_pilot/molecular_by_body.py``); the step prints it since mbe_step
  2026.10.6.3.
- **Net forces by order:** monomers at most 10.4 meV/Å (m27, the ORCA − VASP
  monomer increment), pairs 1.9 meV/Å, triples 4.7 meV/Å.
- **The many-body terms are smaller than EC's:** the 2-body term is 72% of EC's
  and the 3-body term 30%. In the molecular virial, the 2- and 3-body terms
  again largely cancel, leaving a correction of −77 atm, close to EC's −72 atm.
- **FEC's cost:** about 4,600 core-h for the final run, which recomputed the cell
  and 728 high-level tasks and reused the rest. The six attempts used about
  7,230 core-h in all:

  ==================  =========
  level               core-h
  ==================  =========
  high (triples)      5,280
  high (pairs)        825
  cell                425
  molecular triples   413
  periodic monomers   163
  molecular pairs     126
  ==================  =========

  These are sacct CPU times of the bundles whose logs remain; bundles cancelled
  in the first attempts are partly missing.

The cell's SCF
--------------

VASP's minimizer for the whole cell needed three changes to converge both
carbonates:

- **ALGO = All with the new line search (ISEARCH = 1)**, vasp-step up to
  2026.10.6, converged EC in 30 steps. On FEC it aborted at SCF step 6, the
  first self-consistent step, with "EDWAV: internal error, the gradient is not
  orthogonal", on 24 and 32 ranks alike.
- **ALGO = Normal (Davidson) with VASP's default mixing**, vasp-step 2026.10.7,
  diverged on FEC. The density residual grew from 1.5 to 7.6 over steps 10–30,
  and the energy swung by up to 500 eV a step. It was cancelled after 6.5 hours.
- **A three-way test on the FEC cell's own inputs** (64 ranks, TinkerCliffs
  7884921):

  ===============================  =======================================
  setting                          FEC
  ===============================  =======================================
  ALGO = All, ISEARCH = 0          converged, 44 steps, −2210.922324 eV
  Normal, AMIX 0.2, BMIX 0.0001,   converged, 24 steps, −2210.922324 eV
  MAXMIX 50
  ALGO = Damped, TIME 0.2          ~100 eV off after 10 steps; cancelled
  ===============================  =======================================

- **The EC check** (7890904, 7899733):

  - Normal with the gentler mixing came within 0.4 eV at step 10, then
    collapsed: by step 23 the energy was −1.4 × 10⁶ eV.
  - ALGO = All with ISEARCH = 0 converged in 33 steps to −2153.9255 eV, the
    energy of EC's original run with ISEARCH = 1 (30 steps), so the legacy line
    search costs EC almost nothing.

So since vasp-step 2026.10.8 a whole cell uses ALGO = All with ISEARCH = 0, and
fragments keep ISEARCH = 1. In the pilot the FEC cell converged in 27 steps, in
2 h 22 min on 64 ranks, to −2210.92232 eV. The cell's time limit is 3× the
hand-fitted estimate; on 32 ranks that is 11:15 against about 8 hours, so the
pilot used 64 ranks (an 8-hour limit against the measured 2.4 hours).

The runs
--------

1. **Attempts 1–2** (32 and 24 cell ranks): the cell aborted with EDWAV on ALGO =
   All. The fix was vasp-step 2026.10.7 (ALGO = Normal for cells).
2. **Attempt 3:** seamm-exec 2026.10.6's cost model was fitted only on science's
   timing rows, and it predicted about 13 s per ORCA task. So it packed about 240
   tasks into each 51-minute bundle. The ORCA bundles were cancelled. The fix was
   seamm-exec 2026.10.7, which refuses to predict outside the fitted range.
3. **Attempt 4** sat idle. The tasks lost in attempt 3's timeouts had used their
   three attempts, so they were never resubmitted. Their manifest entries (state
   lost, attempts ≥ 3) were deleted, keeping backups
   ``*/tasks/manifest.json.before-reset-2026-10-07``.
4. **Attempt 5** (7881654) had two problems:

   - the cell diverged on ALGO = Normal (above);
   - 22 high and 706 high-triples tasks could never run. The bundles cancelled in
     attempt 3 had left them with state cancelled and three attempts, and the
     reset had removed only the lost ones.

   It was cancelled, and those entries were deleted too (backups
   ``*.before-reset2-2026-10-07``).
5. **Attempt 6** (7903005) ran on vasp-step 2026.10.8 with 64 cell ranks:

   - the ORCA triples went one per 4-hour bundle;
   - seven tasks timed out on a contended node, with six triples on tc307 in the
     MP2 gradient at about 26% CPU, and passed on retry;
   - it finished in 5.5 hours.

Lessons
-------

- **Resetting after a cancellation:** delete the manifest entries that are not
  finished and have used their attempts, whatever their state, not only the
  lost ones. Check first for a DONE marker. Tasks that finished before the
  cancellation can still show "cancelled" or no state, but their DONE is reused.
  seamm-exec dev f9aba39 (unreleased) fixes both: cancelling collects tasks
  whose DONE or FAILED exists, and a cancelled task gets its attempt back.
- **Timeouts counted as attempts,** so a low estimate burned a task's three
  attempts. In seamm-exec dev f9aba39 (unreleased) the first ``max_timeouts``
  (default 4) timeouts don't count; each retry still doubles the time.
- **Bundles in a queue wrote no wall time** in their timing records, so none of
  this pilot's ORCA rows could be fitted. Fixed in seamm-exec 2026.10.8.1.
- **A converged cell SCF on one carbonate says little about another.** Check a
  new cell's SCF on its own inputs before a pilot, e.g. 15 steps on 64 ranks,
  about an hour.
