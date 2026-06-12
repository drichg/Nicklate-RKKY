# Sunny Bilayer Ground-State And Spectrum Design

## Goal

Consume the completed bilayer \(x^2-y^2\) RKKY result and, for every
\(U=0.0:0.1:6.0\), calculate:

- one minimized classical ground-state candidate on a \(4\times4\) bilayer;
- the linear spin-wave dispersion;
- the `ssf_perp` spectral weight;
- a weighted dispersion plot;
- complete numerical data and stability diagnostics.

This phase begins only after the RKKY extraction plan has produced a complete,
validated JLD2 input.

## Input Contract

Read

```text
data/results/rkky_x2_bilayer_U0-6_Nk100.jld2
```

and require all 61 requested \(U\) points. For each point, use the physical
exchange constants in eV:

\[
J_0,\quad J_1,\quad J_2,\quad J_3,\quad J'_1.
\]

The RKKY metadata, requested \(U\) list, basis, units, and exchange arrays are
copied into the Sunny output metadata. Malformed, incomplete, or incompatible
RKKY input must be rejected before any Sunny calculation starts.

## Spin Hamiltonian

Use the isotropic bilayer Heisenberg model

\[
\begin{aligned}
H={}&J_0\sum_i\mathbf S_{i1}\cdot\mathbf S_{i2}\\
&+J_1\sum_{\langle ij\rangle,l}\mathbf S_{il}\cdot\mathbf S_{jl}
+J_2\sum_{\langle\langle ij\rangle\rangle,l}
\mathbf S_{il}\cdot\mathbf S_{jl}\\
&+J_3\sum_{\langle ij\rangle_3,l}\mathbf S_{il}\cdot\mathbf S_{jl}
+J'_1\sum_{\langle ij\rangle,l\ne l'}
\mathbf S_{il}\cdot\mathbf S_{jl'}.
\end{aligned}
\]

Sunny uses the same sign convention: \(J>0\) is antiferromagnetic for
\(H=\sum J_{ij}\mathbf S_i\cdot\mathbf S_j\). Do not alter the signs read from
the RKKY file.

## Crystal And Spin System

Use eV and angstrom throughout:

```julia
Units(:eV, :angstrom)
```

The AA-stacked bilayer chemical cell has

\[
a=b=3.85\ {\rm \AA},\qquad c=20.0\ {\rm \AA},
\]

with two Ni positions

\[
\tau_1=(0,0,0.25),\qquad
\tau_2=(0,0,0.75).
\]

Each site carries

\[
S=\frac12,\qquad g=2.
\]

Create a dipole system with

```text
dims = (4,4,1)
```

which gives 32 spins.

The representative Sunny bonds are:

- \(J_0\): `Bond(1, 2, [0, 0, 0])`;
- \(J_1\): `Bond(1, 1, [1, 0, 0])` and the equivalent layer-2 bond;
- \(J_2\): `Bond(1, 1, [1, 1, 0])` and the equivalent layer-2 bond;
- \(J_3\): `Bond(1, 1, [2, 0, 0])` and the equivalent layer-2 bond;
- \(J'_1\): `Bond(1, 2, [1, 0, 0])`.

Sunny symmetry propagation supplies equivalent in-plane directions. Tests
must inspect the resulting reference bonds so accidental missing or duplicate
interactions are detected.

## Ground-State Candidate

Each \(U\) is independent.

1. Create a fresh 32-spin system.
2. Seed the RNG with
   `base_seed + U_index`, where `base_seed = 20260612`.
3. Call `randomize_spins!`.
4. Record the initial energy.
5. Call `minimize_energy!` once.
6. Record the final total energy, energy per site, and all 32 spin vectors.

Do not inherit the previous \(U\) configuration and do not perform multiple
random restarts. This workflow provides one deterministic, reproducible local
minimum per \(U\), not a proof of the global classical ground state.

Non-finite energies or malformed spin vectors are processing failures and
must not be marked complete.

## Spin-Wave Spectrum

Create conventional linear spin-wave theory using

```julia
measure = ssf_perp(sys)
swt = SpinWaveTheory(sys; measure)
```

Use the path

\[
\Gamma(0,0,0)
\rightarrow X(1/2,0,0)
\rightarrow M(1/2,1/2,0)
\rightarrow\Gamma(0,0,0)
\]

in reciprocal lattice units of the chemical cell, sampled at 301 total
points.

Save:

- the sampled \(q\) points and path distance;
- every magnon branch energy in eV;
- the corresponding `ssf_perp` mode weight;
- high-symmetry tick positions and labels.

The primary plot is a dispersion plot whose color or marker intensity
represents spectral weight. No artificial broadening or energy-bin intensity
map is required in the first implementation.

## Stability Classification

Inspect the raw spin-wave frequencies without taking absolute values.

For a stable result, save:

```text
status = "ok"
```

If a significant negative frequency occurs, save:

```text
status = "unstable_spin_wave"
```

and record:

- minimum frequency;
- negative-mode count;
- unstable \(q\)-point indices;
- unstable mode indices;
- complete raw frequency and weight arrays;
- ground-state energy and spin configuration.

An unstable point remains a processed point and does not stop the 61-point
scan. Generate a diagnostic plot with `UNSTABLE` in the title and retain the
negative branches.

The negative-frequency tolerance must be explicit metadata and tested. Small
roundoff-level negatives inside the tolerance are not classified as
instability.

## Output And Resume

Save all numerical data in one file:

```text
data/results/sunny_bilayer_U0-6.jld2
```

Write one plot per \(U\):

```text
results/excitations/dispersion_weight_U0.0.png
...
results/excitations/dispersion_weight_U6.0.png
```

After each processed \(U\), atomically replace the formal JLD2 result using a
temporary file in the same directory.

On restart:

- validate the entire Sunny metadata and source RKKY contract;
- treat both `"ok"` and `"unstable_spin_wave"` points as processed;
- continue from the first unprocessed \(U\);
- if a numerical result exists but its PNG is missing, regenerate only the
  plot;
- recompute everything only when `overwrite=true`.

## JLD2 Schema

Metadata includes:

- schema version;
- source RKKY path and copied metadata;
- source exchange arrays or a content fingerprint;
- Sunny package version;
- energy and length units;
- lattice constants and basis positions;
- spin and \(g\);
- supercell dimensions;
- base random seed;
- \(q\)-path vertices, labels, and 301-point count;
- negative-frequency tolerance;
- Hamiltonian and bond definitions.

Per-\(U\) results include:

- \(U\) and five exchanges;
- deterministic seed;
- initial and final total energies;
- final energy per site;
- 32 final spin vectors;
- sampled \(q\) points and path distances;
- mode energies and `ssf_perp` weights;
- status and stability diagnostics.

All completed arrays must have equal lengths, use the requested \(U\) prefix,
and contain finite values except that negative finite frequencies are allowed
and explicitly classified.

## Code Boundaries

Reusable functionality is separated into:

```text
src/ground_state/sunny_bilayer.jl
src/excitations/sunny_spin_wave.jl
```

The executable orchestration is:

```text
scripts/sunny/run_bilayer_groundstate_spectrum.jl
```

The ground-state module owns crystal construction, bond assignment, seeded
initialization, minimization, and spin serialization. The excitation module
owns path construction, spin-wave calculation, mode weights, stability
classification, plotting data, persistence, and resume orchestration.

Sunny must be added as an explicit project dependency and pinned through the
Manifest. Plotting may use the existing plotting dependency or the Makie
extension recommended by the installed Sunny version, but GUI calls such as
`view_crystal` must not appear in the batch workflow.

## Validation

Tests must cover:

1. the AA bilayer cell, 32-site supercell, and \(S=1/2,g=2\);
2. all five representative exchange bonds and propagated equivalents;
3. a simple FM/AFM parameter set with analytically predictable energy sign;
4. repeatability under the same deterministic seed;
5. a 301-point \(\Gamma-X-M-\Gamma\) path with correct endpoints and ticks;
6. extraction and storage of dispersion and `ssf_perp` weights;
7. stable and unstable classification using synthetic frequency arrays;
8. continuation after an unstable point;
9. one-file atomic resume and incompatible-input rejection;
10. missing-image regeneration without repeating minimization;
11. script inclusion without launching the full 61-point scan.

Implementation tests should use a smaller supercell or minimal synthetic data
where possible. At least one Sunny integration test must build the real
\(4\times4\) bilayer and calculate a small spectrum.
