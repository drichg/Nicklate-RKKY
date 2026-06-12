# Project Reorganization Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Reorganize the research folder into modular Julia source, scripts, data, results, paper, and reference directories without changing the intended numerical calculations.

**Architecture:** Add a lightweight `NickelateRKKY` loader that includes model and susceptibility source files, while keeping optional Sunny excitation work outside the default load path. Make every runtime path independent of the current working directory, then move artifacts by role and remove reproducible build products.

**Tech Stack:** Julia, Test, JLD2, Plots, PowerShell.

---

### Task 1: Add Structure And Path Tests

**Files:**
- Create: `tests/test_project_structure.jl`
- Modify: existing tests after moves

1. Add tests for the main module loader and project-root/data path constants.
2. Run the new test and confirm it fails because the new loader does not exist.
3. Add the source hierarchy and split the model from susceptibility code.
4. Run the focused test and existing chemical-potential test.

### Task 2: Move Executable Scripts

**Files:**
- Move calculation entry points to `scripts/susceptibility/`
- Move plot entry points to `scripts/plotting/`
- Move the Sunny draft to `src/excitations/`

1. Update tests to target the planned script paths and confirm failure before
   moving the scripts.
2. Move scripts and replace fragile relative includes with the main loader.
3. Anchor default input and output paths at the project root.
4. Run focused script tests.

### Task 3: Sort Data And Research Artifacts

**Files:**
- Move model tables to `data/model/`
- Move tabular RKKY inputs to `data/rkky/`
- Move JLD2 files to `data/results/`
- Move PNG outputs to `results/`
- Move TeX/PDF/PPT files to `paper/`
- Move downloaded references to `references/`

1. Move preserved artifacts without deleting scientific data.
2. Update source defaults to the new locations.
3. Confirm all expected files exist in their destination directories.

### Task 4: Add Project Documentation And Ignore Rules

**Files:**
- Create: `README.md`
- Create: `.gitignore`
- Create: `src/ground_state/README.md`

1. Document module boundaries, primary commands, data locations, and future
   extension points.
2. Ignore Julia caches, LaTeX build products, and generated result patterns
   while retaining existing results.

### Task 5: Remove Reproducible Clutter

**Files:**
- Delete LaTeX intermediate files
- Delete `Manifest.toml.bak`
- Delete empty `rslt/`

1. Verify every deletion candidate matches the approved cleanup policy.
2. Remove only approved intermediates and empty directories.
3. List the final root directory and check no loose generated artifacts remain.

### Task 6: Full Verification

1. Run all Julia test files from the project root.
2. Load `src/NickelateRKKY.jl` from a different working directory to verify
   path independence.
3. Report pre-existing scientific test disagreements separately from
   organization failures.

