# Simple High-Symmetry Susceptibility Plot Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Simplify `plot_chi_highsym.jl` so the user edits `U` and `component = 11/12/...` at the top and runs the script.

**Architecture:** Keep the script include-safe for tests, but make direct use simple: top-level variables call one plotting function only when the file is executed directly. The only supported component format is a two-digit integer matrix element.

**Tech Stack:** Julia, Test, Plots, existing `susceptibility.jl`.

---

### Task 1: Replace Tests

**Files:**
- Modify: `tests/test_plot_chi_highsym.jl`

**Step 1:** Test `component_indices(12) == (1, 2)` and invalid components throw.

**Step 2:** Test the path endpoints.

**Step 3:** Test a small `Nk=4`, `n_per_segment=1` calculation and PNG save.

**Step 4:** Run `julia --project=. tests/test_plot_chi_highsym.jl` and verify it fails against the current complex script.

### Task 2: Replace Script

**Files:**
- Modify: `plot_chi_highsym.jl`

**Step 1:** Put editable variables at the top.

**Step 2:** Implement only the helpers needed by the tests.

**Step 3:** Run `julia --project=. tests/test_plot_chi_highsym.jl` and verify it passes.

### Task 3: Verify Related Test

Run `julia --project=. tests/test_mu.jl`.
