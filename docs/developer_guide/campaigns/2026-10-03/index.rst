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

Phase-2 limitations (open items)
================================

- The levels run one after the other (see Q1).
- No real periodic run yet. No installed provider declares
  ``stress_convention``, and VASP has no ``get_task`` until phase 3. The cell
  path is covered by the pilot regression. A MOPAC-periodic toy needs mopac_step
  to declare "pressure".
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
