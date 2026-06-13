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
- The production `Nk=Nq=100`, `U=0:0.1:6` scan completed all 61 points in
  approximately 40 minutes and produced
  `data/results/rkky_x2_bilayer_U0-6_Nk100.jld2`.
- Production integrity diagnostics:
  - maximum filling error: `9.8843e-11`;
  - maximum intralayer symmetry residual: `9.8810e-15`;
  - maximum interlayer symmetry residual: `1.2490e-16`;
  - maximum Fourier imaginary residual: `1.4097e-15`.
- Sunny 0.8.0 was added as a direct project dependency. This version has no
  native eV unit, so exchange and frequency values are converted at the Sunny
  boundary while all persisted and reported values remain in eV.
- A real `U=0` initial 4x4 attempt minimized to a period-4 state with
  antialigned layers and `E/N=-0.00437425` eV, but Sunny reported an unstable
  mode near the magnetic Gamma point.
- The full nine-attempt retry protocol at `U=0` remained unstable for five
  4x4 seeds, two 6x6 seeds, and two 8x8 seeds. The lowest-energy 4x4 and 8x8
  states agree to numerical precision and retain period-4 order; this point
  requires physical interpretation after the full scan.
- Fixed a direct-script namespace bug found before production: the script
  referenced `Sunny` from `Main` although the dependency was imported inside
  `NickelateRKKY`. A regression-tested `sunny_runtime_version()` API now keeps
  the script independent of that module-loading detail.
- The formal Sunny scan completed all 61 U points and wrote
  `data/results/sunny_bilayer_U0-6.jld2`. Every point remained classified as
  `unstable_spin_wave` after the full retry protocol.
- The selected finite-cell states are period-4 throughout the scan. Their
  layer correlation changes from -1 for `U <= 0.9` eV to +1 for
  `U >= 1.0` eV, locating the sampled transition in `0.9 < U < 1.0` eV.
- Fifty-six of the 61 selected states did not trigger a Sunny minimizer
  convergence warning. The repeated instability of converged states shows
  that minimizer iteration count alone is not the dominant cause.
- A continuous exchange-matrix minimization gives an incommensurate diagonal
  ordering vector between approximately `(0.273, 0.273)` and
  `(0.289, 0.289)` RLU. This is incompatible with the quarter-grid momenta of
  a 4x4 cell and explains the near-Gamma negative-curvature modes reported by
  Sunny. The interpretation is consistent with the incommensurate LSWT
  treatment of Toth and Lake (2015) and Sunny's requirement that the supplied
  magnetic structure be a true energy minimum.
- RKKY convergence checks were run at `U = 0.0, 2.9, 3.0, 3.6, 6.0`.
  The `Nk=80` versus 100 threshold triggered `Nk=120` diagnostics at
  `U=0.0, 2.9, 3.6`. The 120-point values remain close to the 100-point
  results; the largest relative sensitivities occur in exchange components
  that are close to zero.
- Generated 61 instability diagnostic images, four summary figures, a compact
  report-data JLD2 file, and the five-page Chinese report
  `paper/rkky_sunny_results.pdf`.
- GitHub pushes were attempted after the first milestone. The sandboxed and
  escalated execution environments reached GitHub but did not expose the SSH
  identity accepted in the user's interactive PowerShell session, so local
  commits remain the authoritative backup until the user runs `git push`.
