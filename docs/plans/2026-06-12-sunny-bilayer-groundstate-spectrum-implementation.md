# Sunny Bilayer Ground-State And Spectrum Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Read the completed five-exchange RKKY scan and produce a resumable 61-point Sunny calculation of deterministic randomized ground-state candidates, spin-wave dispersions, `ssf_perp` weights, stability records, and per-\(U\) plots.

**Architecture:** Add one module for the AA bilayer spin model and ground-state minimization, and one module for spin waves, stability diagnostics, plotting, persistence, and scan orchestration. Keep the default 61-point run in an include-safe script and use a single atomically replaced JLD2 result file.

**Tech Stack:** Julia 1.11, Sunny 0.9.x, JLD2, Test, plotting backend compatible with Sunny, existing NickelateRKKY project.

---

## Prerequisite

Do not begin this plan until
`docs/plans/2026-06-12-x2-bilayer-rkky-implementation.md` is complete and the
file

```text
data/results/rkky_x2_bilayer_U0-6_Nk100.jld2
```

passes its metadata and completeness validation.

### Task 1: Add And Verify Sunny

**Files:**
- Modify: `Project.toml`
- Modify: `Manifest.toml`
- Create: `tests/test_sunny_bilayer.jl`

**Step 1: Add a failing dependency smoke test**

Create a focused test that imports Sunny and checks the APIs required by the
design are defined:

```julia
using Test
using Sunny

@test isdefined(Sunny, :System)
@test isdefined(Sunny, :set_exchange!)
@test isdefined(Sunny, :randomize_spins!)
@test isdefined(Sunny, :minimize_energy!)
@test isdefined(Sunny, :SpinWaveTheory)
@test isdefined(Sunny, :ssf_perp)
@test isdefined(Sunny, :intensities_bands)
```

**Step 2: Run the test to verify RED**

Run:

```powershell
julia --project=. tests/test_sunny_bilayer.jl
```

Expected: fail because Sunny is not a project dependency.

**Step 3: Add Sunny through Julia Pkg**

Run:

```powershell
julia --project=. -e "using Pkg; Pkg.add(name=\"Sunny\", version=\"0.9\")"
```

Add the plotting dependency required by the selected Sunny 0.9 plotting path
only if it is not already provided by the project. Do not install packages
globally.

**Step 4: Run the smoke test**

Expected: all API assertions pass. Record the exact installed Sunny version.

**Step 5: Commit**

```powershell
git add Project.toml Manifest.toml tests/test_sunny_bilayer.jl
git commit -m "build: add Sunny spin-wave dependency"
```

### Task 2: Build The AA Bilayer Spin Model

**Files:**
- Create: `src/ground_state/sunny_bilayer.jl`
- Modify: `src/NickelateRKKY.jl`
- Modify: `tests/test_sunny_bilayer.jl`

**Step 1: Write failing crystal and system tests**

Define a wished-for API:

```julia
model = build_sunny_bilayer(
	J0=0.01,
	J1=0.02,
	J2=-0.003,
	J3=0.004,
	J1p=0.005,
	dims=(4, 4, 1),
)
```

Assert:

- energy units are eV and length units are angstrom;
- lattice constants are `(3.85, 3.85, 20.0)`;
- basis positions are `[0,0,0.25]` and `[0,0,0.75]`;
- the system has 32 sites;
- both moments have `s=1/2` and `g=2`;
- the model records all five exchange values.

Also inspect Sunny reference bonds or symmetry tables to verify representative
bonds for `J0`, `J1`, `J2`, `J3`, and `J1p` exist without duplicate
assignment.

**Step 2: Run the test to verify RED**

Expected: fail because `build_sunny_bilayer` does not exist.

**Step 3: Implement the model builder**

Use:

```julia
units = Units(:eV, :angstrom)
latvecs = lattice_vectors(3.85, 3.85, 20.0, 90, 90, 90)
positions = [[0.0, 0.0, 0.25], [0.0, 0.0, 0.75]]
cryst = Crystal(latvecs, positions; types=["Ni1", "Ni2"])
sys = System(
	cryst,
	[1 => Moment(s=1/2, g=2), 2 => Moment(s=1/2, g=2)],
	:dipole;
	dims,
	units,
)
```

Assign:

```julia
Bond(1, 2, [0, 0, 0]) # J0
Bond(1, 1, [1, 0, 0]) # J1, layer 1
Bond(2, 2, [1, 0, 0]) # J1, layer 2
Bond(1, 1, [1, 1, 0]) # J2, layer 1
Bond(2, 2, [1, 1, 0]) # J2, layer 2
Bond(1, 1, [2, 0, 0]) # J3, layer 1
Bond(2, 2, [2, 0, 0]) # J3, layer 2
Bond(1, 2, [1, 0, 0]) # J1p
```

If Sunny symmetry makes explicit layer-2 assignments redundant or rejects
them, use `reference_bonds` to choose the minimal non-duplicating set and
update the tests to assert the resulting physical bonds, not incidental call
count.

**Step 4: Run the focused test**

Expected: crystal, site-count, moment, and exchange-bond tests pass.

**Step 5: Commit**

```powershell
git add src/NickelateRKKY.jl src/ground_state/sunny_bilayer.jl tests/test_sunny_bilayer.jl
git commit -m "feat: build Sunny AA bilayer spin model"
```

### Task 3: Add Deterministic Ground-State Minimization

**Files:**
- Modify: `src/ground_state/sunny_bilayer.jl`
- Modify: `tests/test_sunny_bilayer.jl`

**Step 1: Write failing minimization tests**

Call:

```julia
result1 = minimize_sunny_bilayer(
	J0=...,
	J1=...,
	J2=...,
	J3=...,
	J1p=...,
	seed=20260613,
)
```

Assert:

- `initial_energy`, `final_energy`, and `energy_per_site` are finite;
- the final energy does not exceed the initial energy beyond tolerance;
- `spins` has shape `(32, 3)`;
- every spin vector has the Sunny dipole magnitude expected for \(S=1/2\);
- running again with the same seed and parameters reproduces energies and
  spins within tight tolerance.

Add a simple ferromagnetic or antiferromagnetic parameter case whose energy
sign can be predicted independently, to verify exchange sign convention.

**Step 2: Run the test to verify RED**

Expected: fail because `minimize_sunny_bilayer` does not exist.

**Step 3: Implement one-start minimization**

The function must:

1. build a fresh model;
2. create a local deterministic RNG or seed the exact RNG used by Sunny;
3. call `randomize_spins!`;
4. record the initial energy;
5. call `minimize_energy!` exactly once;
6. serialize the final dipoles in stable site order;
7. validate finite energies and spin vectors.

Do not add multiple restarts or previous-\(U\) inheritance.

**Step 4: Run the focused test**

Expected: deterministic minimization tests pass.

**Step 5: Commit**

```powershell
git add src/ground_state/sunny_bilayer.jl tests/test_sunny_bilayer.jl
git commit -m "feat: minimize deterministic Sunny bilayer states"
```

### Task 4: Build The 301-Point Spin-Wave Path

**Files:**
- Create: `src/excitations/sunny_spin_wave.jl`
- Modify: `src/NickelateRKKY.jl`
- Modify: `tests/test_sunny_bilayer.jl`

**Step 1: Write failing path tests**

Define:

```julia
path_data = sunny_bilayer_q_path(cryst; n_points=301)
```

Assert:

- exactly 301 sampled points;
- endpoints and vertices are
  `Gamma=[0,0,0]`, `X=[1/2,0,0]`, `M=[1/2,1/2,0]`, `Gamma`;
- tick labels are `["Gamma", "X", "M", "Gamma"]`;
- distances are monotone;
- all points have zero out-of-plane component.

**Step 2: Run the test to verify RED**

Expected: fail because the path helper does not exist.

**Step 3: Implement using Sunny `q_space_path`**

Keep both the Sunny path object and serializable arrays of sampled \(q\)
points, distances, ticks, and labels.

**Step 4: Run the focused test**

Expected: path tests pass.

**Step 5: Commit**

```powershell
git add src/NickelateRKKY.jl src/excitations/sunny_spin_wave.jl tests/test_sunny_bilayer.jl
git commit -m "feat: add Sunny bilayer spin-wave path"
```

### Task 5: Calculate Dispersion And ssf_perp Weights

**Files:**
- Modify: `src/excitations/sunny_spin_wave.jl`
- Modify: `tests/test_sunny_bilayer.jl`

**Step 1: Write a failing small-spectrum integration test**

Use a simple exchange set and a minimized system. Call:

```julia
spectrum = compute_sunny_bilayer_spectrum(sys; n_q=31)
```

Assert:

- mode energies and weights share the same \(q\) and branch dimensions;
- all weights are finite and nonnegative within tolerance;
- energies retain their raw sign;
- metadata says `measure = "ssf_perp"`;
- units are eV.

Use Sunny 0.9's `SpinWaveTheory`, `ssf_perp`, and `intensities_bands`. Adapt
only to the installed public return type, and convert it immediately to
plain serializable arrays.

**Step 2: Run the test to verify RED**

Expected: fail because `compute_sunny_bilayer_spectrum` does not exist.

**Step 3: Implement the spectrum function**

Create:

```julia
measure = ssf_perp(sys)
swt = SpinWaveTheory(sys; measure)
```

Compute the band intensities on the path, then extract branch energies and
weights without broadening.

**Step 4: Run the focused test**

Expected: spectrum shape and finite-value tests pass.

**Step 5: Commit**

```powershell
git add src/excitations/sunny_spin_wave.jl tests/test_sunny_bilayer.jl
git commit -m "feat: calculate Sunny dispersion and spectral weights"
```

### Task 6: Classify Stable And Unstable Spectra

**Files:**
- Modify: `src/excitations/sunny_spin_wave.jl`
- Modify: `tests/test_sunny_bilayer.jl`

**Step 1: Write failing pure classification tests**

Use synthetic arrays:

```julia
stable = classify_spin_wave([0.0 0.1; 0.2 0.3]; negative_tolerance=1e-8)
unstable = classify_spin_wave([0.0 -0.01; 0.2 0.3]; negative_tolerance=1e-8)
```

Assert:

- roundoff-level negatives are `"ok"`;
- significant negatives are `"unstable_spin_wave"`;
- minimum frequency and negative count are correct;
- exact \(q\) and branch indices are returned.

**Step 2: Run the test to verify RED**

Expected: fail because the classifier does not exist.

**Step 3: Implement classification without modifying energies**

Never replace negative values by absolute values or zero. Return a pure,
serializable diagnostic record.

**Step 4: Integrate classification into the spectrum result**

Ensure real Sunny output is classified and retained even when unstable.

**Step 5: Commit**

```powershell
git add src/excitations/sunny_spin_wave.jl tests/test_sunny_bilayer.jl
git commit -m "feat: record unstable Sunny spin-wave modes"
```

### Task 7: Plot Weighted Dispersions

**Files:**
- Modify: `src/excitations/sunny_spin_wave.jl`
- Modify: `tests/test_sunny_bilayer.jl`

**Step 1: Write failing plot tests**

With synthetic stable and unstable spectra, write PNG files into a temp
directory and assert:

- files exist and are nonempty;
- stable title includes \(U\);
- unstable title includes `UNSTABLE`;
- negative branches remain visible in unstable plot input;
- color or marker intensity represents the saved mode weights.

**Step 2: Run the test to verify RED**

Expected: fail because the plot helper does not exist.

**Step 3: Implement plotting**

Use a noninteractive backend. Plot every branch against path distance, encode
`ssf_perp` weight, draw high-symmetry separators, label energy in eV, and use
the exact saved ticks.

**Step 4: Run the focused test**

Expected: PNG tests pass without opening a GUI.

**Step 5: Commit**

```powershell
git add src/excitations/sunny_spin_wave.jl tests/test_sunny_bilayer.jl
git commit -m "feat: plot weighted Sunny dispersions"
```

### Task 8: Define Sunny Result Metadata And Input Validation

**Files:**
- Modify: `src/excitations/sunny_spin_wave.jl`
- Modify: `tests/test_sunny_bilayer.jl`

**Step 1: Write failing metadata tests**

Construct synthetic complete RKKY input and assert Sunny metadata contains:

- source path and source metadata;
- all source exchange arrays or deterministic fingerprint;
- Sunny version;
- units, lattice, positions, spin, \(g\), and dimensions;
- base seed;
- path vertices, labels, and 301 points;
- negative-frequency tolerance;
- bond definitions.

Test rejection of incomplete \(U\) data, changed exchange values, changed
supercell dimensions, and changed tolerance.

**Step 2: Run the test to verify RED**

Expected: metadata helpers are missing.

**Step 3: Implement strict validation**

Require the exact 61-point \(U=0:0.1:6\) contract for the default script.
Reusable lower-level scan functions may accept smaller test inputs but must
still validate prefix and exchange-array consistency.

**Step 4: Run the focused test**

Expected: metadata tests pass.

**Step 5: Commit**

```powershell
git add src/excitations/sunny_spin_wave.jl tests/test_sunny_bilayer.jl
git commit -m "feat: validate Sunny scan input metadata"
```

### Task 9: Add Atomic Persistence And Resume

**Files:**
- Modify: `src/excitations/sunny_spin_wave.jl`
- Modify: `tests/test_sunny_bilayer.jl`

**Step 1: Write failing persistence tests**

Using synthetic point calculators and plot writers:

1. save a one-point `"ok"` result;
2. save a one-point `"unstable_spin_wave"` result;
3. resume and verify both statuses are skipped;
4. verify the next missing \(U\) is calculated;
5. reject changed RKKY input or Sunny metadata before calculation;
6. confirm no temporary file remains after success;
7. delete a PNG and confirm resume redraws it without invoking minimization;
8. confirm `overwrite=true` recalculates every point.

**Step 2: Run the test to verify RED**

Expected: scan persistence helpers do not exist.

**Step 3: Implement persistence and orchestration**

Provide:

```julia
save_sunny_result_atomic(filepath, result)
load_sunny_result(filepath)
scan_sunny_bilayer(rkky_result; output_path, image_dir, overwrite=false, ...)
```

Persist plain arrays and named tuples, not Sunny internal objects.

**Step 4: Run the focused test**

Expected: resume, unstable continuation, overwrite, and missing-image tests
pass.

**Step 5: Commit**

```powershell
git add src/excitations/sunny_spin_wave.jl tests/test_sunny_bilayer.jl
git commit -m "feat: add resumable Sunny bilayer scans"
```

### Task 10: Add The Default 61-Point Script

**Files:**
- Create: `scripts/sunny/run_bilayer_groundstate_spectrum.jl`
- Modify: `README.md`
- Modify: `tests/test_sunny_bilayer.jl`

**Step 1: Write a failing include-safety test**

Include the script and assert its default configuration:

```julia
input_path = data/results/rkky_x2_bilayer_U0-6_Nk100.jld2
output_path = data/results/sunny_bilayer_U0-6.jld2
image_dir = results/excitations
dims = (4,4,1)
base_seed = 20260612
n_q = 301
```

Including the script must not load the large RKKY file or launch Sunny
calculations.

**Step 2: Run the test to verify RED**

Expected: script is missing.

**Step 3: Implement the guarded executable**

On direct execution:

1. load and validate the completed RKKY file;
2. print Sunny version and full configuration;
3. print the number of resumed and pending points;
4. process every pending \(U\);
5. report exchange values, seed, final energy per site, status, minimum
   frequency, and output image;
6. continue after unstable points;
7. print the final JLD2 path and status counts.

Update `README.md` with execution and resume instructions.

**Step 4: Run script parse and include tests**

Expected: script parses and includes without starting the scan.

**Step 5: Commit**

```powershell
git add scripts/sunny/run_bilayer_groundstate_spectrum.jl README.md tests/test_sunny_bilayer.jl
git commit -m "feat: add default Sunny bilayer spectrum scan"
```

### Task 11: Full Verification

**Files:**
- No production edits expected

**Step 1: Run focused Sunny tests**

```powershell
julia --project=. tests/test_sunny_bilayer.jl
```

Expected: all unit and small integration tests pass.

**Step 2: Run RKKY prerequisite tests**

```powershell
julia --project=. tests/test_rkky_x2_bilayer.jl
julia --project=. tests/test_mu.jl
julia --project=. tests/test_run2.jl
```

Expected: all pass.

**Step 3: Run remaining passing regressions**

```powershell
julia --project=. tests/test_plot_rkky_jld2.jl
julia --project=. tests/test_plot_total.jl
julia --project=. tests/test_project_structure.jl
```

Expected: all pass.

**Step 4: Record the known unrelated plotting baseline**

Run:

```powershell
julia --project=. tests/test_plot_chi_highsym.jl
```

The two pre-existing diagonal-versus-layer-combination assertions may remain
failed. Do not modify that code during Sunny implementation.

**Step 5: Perform one real end-to-end smoke point**

Use a copied one-\(U\) RKKY fixture in a temporary directory. Build the real
`(4,4,1)` system, minimize it, calculate a reduced-path spectrum, save and
reload JLD2, and produce one PNG. Do not launch the full 61-point production
scan as a test.

**Step 6: Verify repository scope**

```powershell
git status --short
git diff --check
```

Keep `rkky_bilayer_nickelate_codex.md` untouched unless separately requested.

