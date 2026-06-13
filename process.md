# Calculation Process

All energies, temperatures, frequencies, and exchange constants are expressed
in eV unless explicitly stated otherwise.

## 2026-06-13

- Created and switched to branch `codex/rkky-sunny-workflow`.
- Preserved the pre-existing untracked file `rkky_bilayer_nickelate_codex.md`.
- Baseline tests passed: `test_mu.jl`, `test_run2.jl`,
  `test_plot_rkky_jld2.jl`, `test_plot_total.jl`, and
  `test_project_structure.jl`.
- Baseline `test_plot_chi_highsym.jl` failed two assertions because the
  implementation returned diagonal-only orbital susceptibility while its
  labels and tests specified the layer even combinations.
- Started test-driven implementation of the projected bilayer x2 RKKY
  workflow. The new test initially failed with `x2_layer_channels` undefined,
  confirming the expected RED state.
- Implemented the projected x2 kernel, five Fourier components, strict
  metadata validation, atomic JLD2 replacement, resumable scans, and the
  include-safe production script.
- Corrected `orbital_chi_value` to match its documented and tested layer-even
  combinations. All existing regression tests now pass.
- Twenty-thread projected-kernel benchmark after warmup:
  - `N=20`: 1.5545 s per U.
  - `N=40`: 6.0523 s per U.
  - `N=60`: 14.2192 s per U.
- The observed benchmark suggests that the 61-point `N=100` scan is feasible
  as a multi-hour resumable local calculation.
