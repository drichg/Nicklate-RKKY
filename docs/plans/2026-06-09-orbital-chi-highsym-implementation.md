# Orbital High-Symmetry Susceptibility Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Support `"x2"` and `"z2"` inputs in `plot_chi_highsym.jl` for plotting `chi11+chi13` and `chi22+chi24`.

**Architecture:** Replace the integer component API with an orbital selector and a single projection helper. Preserve the cached momentum-grid calculation and high-symmetry path, while guarding direct execution so the module is safe to include in tests.

**Tech Stack:** Julia, Test, Plots, existing `susceptibility.jl`.

---

### Task 1: Add Failing Orbital Tests

**Files:**
- Modify: `tests/test_plot_chi_highsym.jl`

1. Test selector normalization and invalid input.
2. Test `orbital_chi_value` on a known matrix.
3. Test `compute_chi_path(..., "x2")` and `"z2"` against direct matrix sums.
4. Run `julia --project=. tests/test_plot_chi_highsym.jl`.
5. Confirm failure because the string orbital API does not exist.

### Task 2: Implement the Orbital Plot Utility

**Files:**
- Modify: `plot_chi_highsym.jl`

1. Add `normalize_orbital`, `orbital_chi_value`, and `orbital_label`.
2. Change `compute_chi_path` and `plot_chi_path` to accept an orbital string.
3. Use automatic output names `chi_x2_U*.png` and `chi_z2_U*.png`.
4. Guard the executable configuration with `abspath(PROGRAM_FILE) == @__FILE__`.
5. Run the focused test and confirm it passes.

### Task 3: Regression Verification

1. Run `julia --project=. tests/test_plot_chi_highsym.jl`.
2. Run `julia --project=. tests/test_mu.jl`.
3. Confirm both test files pass without starting a large calculation during include.
