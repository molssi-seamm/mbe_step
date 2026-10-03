=========================================================
2026-10-03: the MBE step, phase 2 (``mbe_step``)
=========================================================

The design is ``~/Work/SEAMM/MBE_correction_step_design.rst`` (its "Decisions
(2026-10-03)" are binding). Phase 1 is the ``seamm_mbe`` library, released as
2026.10.3; see its campaign notes. The plan below was agreed with the "design"
session on 2026-10-03 (questions Q1–Q5).

.. contents::
   :local:

The plan as agreed
==================

Q1, execution
    - **One Evaluator per level** (high, molecular, periodic, cell or
      cluster), each with its own directory, so each has its own manifest and
      bundles. Each passes ``task_set_options`` (archive, bundle size, bundle
      walltime) and the new ``resources``.
    - The facade chooses the path: batch tasks, or a warm MDI engine for
      cheap codes run locally.
    - The levels run **one after the other** in this phase. Running them
      concurrently is a seamm_exec follow-up: a gather over Evaluators sharing
      one LocalPool, owned by the parallel-execution campaign.
    - **seamm_exec change:** ``Evaluator(resources=...)`` is forwarded to the
      provider's ``get_task`` on the batch path. Before this, ORCA batch tasks
      always got one rank. The change rides in the "worker" session's phase-4
      seamm_exec release; mbe_step pins seamm-exec at that version.
Q2, several model chemistries
    - The high level defaults to the flowchart's ``_model_chemistry``.
    - The three low levels are parameters holding level strings, resolved with
      model_chemistry_step's ``resolve_level`` (component-wise ``$variables``),
      ``discover_model_chemistries`` (periodic-only for the periodic and cell
      levels) and ``match_model_chemistry``.
Q3, the sign of the stress
    - Providers declare ``options["stress_convention"]`` ("pressure" or
      "stress") in ``get_model_chemistry_options``, and the step refuses a
      periodic level without it. It is never a GUI choice.
    - The contract paragraph is in seamm_exec's docs, in the same change.
    - MOPAC's MDI engine sends the pressure convention (it negates MOPAC's
      tensile-positive stress), so mopac_step should declare "pressure".
      vasp_step declares its convention in phase 3.
Q4, citations
    - The many-body expansion (Herbert 2019; Richard, Lao and Herbert 2014;
      Beran 2016, all checked against Crossref) and seamm_mbe.
    - The report lists each level's program. Citing each level's code waits
      for a provider ``citations(model_chemistry)`` hook in the batch contract
      (a campaign follow-up).
Q5, validation
    See below.

Layout
======

- ``mbe.py``: the node. It selects the configurations, builds each frame
  (``seamm_mbe.System.from_configuration``, perceiving bonds if there are
  none), resolves the levels, evaluates each level, assembles with seamm_mbe,
  stores the labels and reports.
- ``levels.py``: resolves the model chemistries; ``Level`` wraps the Evaluator;
  fragments are submitted as ``seamm_exec.Geometry`` (with the fragment's charge
  and multiplicity), the cell as its configuration; the stress convention.
- ``labels.py``: no flowchart needed. The cutoff and offset tables, the Voigt
  layout, SEAMM's units, and the idempotent extxyz writer.
- ``tk_mbe.py``: the dialog. It shows only what applies:

  - triples only at order 3, pairs from order 2;
  - single cutoffs or the type tables, not both;
  - the periodic assignment and resources only with a periodic level;
  - order capped at 3, and no counterpoise until phase 3.

- Task keys are ``c<configuration id>-<fragment name>``, and ``-cell`` for the
  whole system.

Validation
==========

**Unit and regression tests (12):**

- The cluster identity through the step's own code. The levels are analytic
  potentials with at most 3-body terms behind a stand-in Evaluator, and every
  pair and triple is selected. Labels equal the high level on the whole cluster:
  energy to 1e-9 eV, forces to 1e-8 eV/Å.
- A failed calculation gives no labels, and the report names it.
- The pilot frame through the step's path, from the stored fragment results as
  Evaluator results (kJ/mol, gradients, the cell's stress in GPa, pressure
  convention):

  - names identical, and 219 + 2564 calculations;
  - P_MBE, P_atom, P_mol and the per-body pressures and energies match
    seamm_mbe's regression (with its exact kB→atm constant);
  - the Voigt [6] σ in GPa reproduces P;
  - the extxyz header and forces are correct.

- The dialog over every layout-driving choice; the helpers; the selection rules
  from the parameters.

**End to end in a private venv** (released packages + the patched seamm_exec +
mbe_step; ``~/SEAMM`` and ``~/SEAMM_DEV`` untouched), with a scratch root holding
copies of ``orca.ini`` and ``mopac.ini``:

- **ORCA, batch path:** a 3-water cluster, B3LYP/def2-SVP high and HF/def2-SVP
  low, every subset. 15 calculations in 15.5 s in the local pool.

  - MBE labels vs the B3LYP trimer's own ``orca.engrad``: energy differs by
    2.7e-9 eV, forces by at most 4.6e-7 eV/Å (SCF noise between the two HF
    runs).
  - A rerun restored all 15 from the manifest in 0.0 s and rewrote the labels
    file without duplicating the frame.
  - ``%maxcore 1500`` (from the 1500 MB per rank) and, with 2 ranks,
    ``%pal nprocs 2`` reach ORCA's input. Identical labels.

- **MOPAC, MDI path:** the same cluster, PM6 high and AM1 low. It ran in
  18 s on warm engines.

**Periodic, real code (MOPAC over MDI) -- a smoke test, not a validation number:**

- Setup: the 64-water pilot geometry (L = 12.4297 Å), default selection
  (pairs < 4.5 Å, connected triples < 3.5 Å). PM7 high, PM6 molecular and
  periodic low, the cell PM6 periodic with its stress. mopac_step declares
  ``stress_convention = "pressure"`` (2026.10.3.2).
- In the same flowchart, the Energy step ran PM7 on the cell directly.
- MBE ran 2,783 calculations in 22 s on warm engines (1282 high, 1282
  molecular, 218 periodic, the cell). Largest increment net force 0.19 meV/Å.
- MBE (PM7 on PM6) vs direct periodic PM7:

  - P_atom −60,507 vs −66,491 atm;
  - energy +10.7 kJ/mol per molecule;
  - forces rms difference 72 meV/Å on an rms of 1,510 meV/Å.

  The signs agree. With a wrong stress sign MBE would give about +40,000
  atm.
- The gap is the truncation of a long-range difference. PM7 has built-in
  dispersion and PM6 has none, so [PM7 − PM6] carries an attractive tail
  beyond the 4.5 Å pairs and the connected triples. Missing it raises the
  energy and makes the pressure less negative, as seen. It is a property of
  this toy pair of levels, not of the step: the production levels (revDSD on
  r2SCAN-D4) both include D4.
- Read Structure needed two workarounds for the frame file (an ASE-style
  extxyz). It needs a ``Properties=`` key, and ``indices = 1:end`` became an
  empty range when it found no frames.

Bugs found on the way
=====================

- **MDI engine name.** The name was "MBE high", and ``MDI_Init`` rejects a
  space in ``-name``. It is now ``MBE_high``.
- **Duplicate labels.** Appending to the extxyz duplicated a frame on a rerun.
  Frames are now replaced by ``SEAMM/configuration_id``.
- **Not this step's bug: MOPAC's MDI engine under an activated venv.** The
  engine is launched with ``conda run -n seamm-mopac python``, but with a
  virtual environment *activated*, ``python`` resolves to the venv's
  interpreter, which lacks mopactools ("mopactools not found"). Without
  activation it works. If the JobServer runs with its venv activated, MOPAC
  over MDI would fail the same way; worth checking in mopac_step
  (``conda run`` with the env's full python path would avoid it).

Review (design, 2026-10-03)
===========================

The review found the core sound. These were fixed, each with a test:

#. **Release order.** ``Evaluator(resources=...)`` is in no released
   seamm_exec yet. It is applied on seamm_exec dev by the "worker" session
   (56ef993; the contract sentence, that the declaration covers both paths,
   is 5f3139e) and rides in its phase-4 release. mbe_step pins that version
   once it is on PyPI. ``test_real_evaluator_with_resources`` runs the real
   Evaluator, so it fails until then. A periodic MOPAC level needs mopac_step
   2026.10.3.2 (PR molssi-seamm/mopac_step#160).
#. **Energy offsets.** The default is now ``none``. Offsets that miss a
   molecule type are refused before any calculation; they used to fail
   after everything ran, storing nothing.
#. **The whole system's charge.** It was the configuration's own charge and
   multiplicity, so a plain-xyz Li⁺ + 2 H₂O ran at charge 0, doublet, while
   its fragments ran Li⁺ at +1. The whole system is now a Geometry with the
   molecules' summed charge, singlet. A configuration charge that disagrees
   is refused.
#. **Help text.** The labels file's help gives ``job:NAME``; a leading ``/``
   is an absolute path.
#. **Configuration charge 0 counts as not set** (plain xyz has none), so the
   molecules' charges win: a Li⁺ structure runs at +1 even if its charge field
   says 0. The user guide says so.
#. **Cell without a stress.** It is marked failed with a reason, instead of
   failing in the assembly.
#. **Order field.** A typed ``$variable`` in it no longer breaks the dialog.
#. **No configurations** is an error.
#. **Contract.** The stress declaration covers both paths (seamm_exec docs;
   mopac_step's docstring and user guide).
#. **Periodic low level.** It accepts any model chemistry in phase 2; the
   help text says so.
#. **Open-shell molecules** in fragments of several molecules are refused
   before running.

Nits:

- the per-body results have no units on the dictionaries;
- "maximum force" is the largest Cartesian component;
- the extxyz reader gives a clear message on a malformed file.

Issues filed for what is not this step's:

- the MOPAC MDI engine under an activated venv,
  `mopac_step#161 <https://github.com/molssi-seamm/mopac_step/issues/161>`_;
- the Energy step's nested stress,
  `energy_step#2 <https://github.com/molssi-seamm/energy_step/issues/2>`_;
- the extxyz reader without ``Properties=``,
  `read_structure_step#81 <https://github.com/molssi-seamm/read_structure_step/issues/81>`_.

The ORCA identity run, repeated with the whole system as a Geometry, gave
identical labels.

Phase-2 limitations (open items)
================================

- The levels run one after the other (see Q1).
- The only real-code periodic check is the MOPAC toy above. VASP has no
  ``get_task`` until phase 3; the pilot regression covers the production
  numbers.
- Periodic fragments are submitted as isolated molecules. The registered-box
  options (n·h boxes, explicit NG, hard PAW, the dipole correction) are the
  phase-3 vasp_step adapter.
- The run-time check that the periodic fragment level and the cell level share
  code, POTCAR, ENCUT and grid is phase 3.
- No counterpoise: phase 3, through seamm_bsse's ghost options.
- The default energy offsets are the water value. Later they should come from
  seamm_mbe's catalog or the thermochemistry database.
- Code citations wait for the provider hook (Q4).
- The Energy step stores its stress as a nested 3×3, which the Write Structure
  extxyz writer (a flat list) would mangle. MBE stores Voigt [6], which works.


Phase 3 (2026-10-03, approved by Paul)
======================================

Order: vasp_step (the batch contract, its own notes), then seamm_mbe (the
``corrections`` hook, local dev 07dd78d), then this step, then the C.1 check on
TinkerCliffs.

- **Grid options.**

  - With a periodic low level, the periodic fragments get
    ``options = {"grid": {"max_spacing", "padding", "reference_cell": cell}}``.
  - The cell gets ``{"grid": {"max_spacing"}}``, so vasp_step registers every
    fragment on the cell's explicit FFT grid.
  - New parameters: *grid spacing* (0.0829 Å) and *box padding* (7.5 Å). The
    dialog shows them only with a periodic level.

- **Consistency.** The cell's low level must be the periodic fragments' (the same
  level string, so the same code, potentials and cutoff, and through the shared
  grid options the same grid). The step refuses otherwise.
- **Pairwise counterpoise** (``counterpoise.py``, the new *counterpoise*
  parameter, shown from order 2):

  - For each selected pair, the two monomers in the pair's basis run as ghost-atom
    tasks at the high level, and at the molecular low level when the pair is
    referenced to it.
  - ``seamm_bsse.combine`` makes the corrected pair; the correction
    (E_high^CP − E_high) − (E_low^CP − E_low), with its forces, enters the sum
    only, through seamm_mbe's ``corrections``.
  - Refused for periodic cells, since the pressure would need the
    energy-scaling virial.
  - Refused for a molecular level that runs through MDI (no ghost atoms there).
  - ``seamm_bsse.combine``'s ghost-gradient fallback is counted and reported.
    ``counterpoise=T`` is written in the extxyz.

- **Tests (22).**

  - A planted BSSE: E^CP = E + 2b per level, so the 2-body sum shifts by exactly
    2(b_high − b_low) per pair, the triples are unchanged, and the forces are
    unchanged with zero ghost gradients.
  - The grid options reach the periodic and cell levels.
  - The three refusals.
  - The dialog's new rules.

- **Real ORCA run** (3 waters, B3LYP high and HF low, def2-SVP, every subset):

  - 13 calculations per molecular level (7 fragments plus 6 ghost jobs);
  - the 2-body term goes from −5.302 to −1.818 kJ/mol per molecule;
  - the 3-body term is unchanged at +0.308;
  - no gradient fallbacks; 20.5 s.
