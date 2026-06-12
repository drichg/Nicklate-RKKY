# High-Symmetry Susceptibility Plot Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a Julia utility that plots user-selected susceptibility along `(0,0) -> (0,pi) -> (pi,pi) -> (0,0)`.

**Architecture:** Add a reusable script with pure helpers for path generation, component parsing, calculation, plotting, and interactive execution. Tests include the script without triggering long computations and verify helper behavior on a small grid.

**Tech Stack:** Julia, Test, Plots, existing `susceptibility.jl`.

---

### Task 1: Write Tests

**Files:**
- Create: `tests/test_plot_chi_highsym.jl`
- Create: `plot_chi_highsym.jl`

**Step 1: Write the failing tests**

Test that:

- the high-symmetry path starts at `(0,0)`, passes `(0,pi)` and `(pi,pi)`, and ends at `(0,0)`
- component inputs `11`, `1,1`, `intra`, and `inter` parse correctly
- a small `Nk=4`, `n_per_segment=1` run returns finite values and can save a plot into a temp directory

**Step 2: Run test to verify it fails**

Run: `julia --project=. tests/test_plot_chi_highsym.jl`

Expected: FAIL because `plot_chi_highsym.jl` does not exist yet.

### Task 2: Implement Helpers

**Files:**
- Modify: `plot_chi_highsym.jl`

**Step 1: Add minimal implementation**

Implement:

- `high_symmetry_chi_path`
- `parse_chi_component`
- `component_label`
- `component_value`
- `chi_highsymmetry_data`
- `plot_chi_highsymmetry`
- `main`

**Step 2: Run test to verify it passes**

Run: `julia --project=. tests/test_plot_chi_highsym.jl`

Expected: PASS.

### Task 3: Verify Existing Tests

**Files:**
- No edits expected.

**Step 1: Run the focused and related tests**

Run: `julia --project=. tests/test_plot_chi_highsym.jl`

Run: `julia --project=. tests/test_mu.jl`

Expected: both PASS.

## Notes

No commit step is included because `G:\Nicklate-RKKY` is not a git repository.
