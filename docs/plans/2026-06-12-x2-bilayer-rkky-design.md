# Bilayer x2 RKKY Extraction Design

## Goal

Compute the five exchange constants

\[
J_0,\quad J_1,\quad J_2,\quad J_3,\quad J'_1
\]

for the bilayer nickelate model from the itinerant
\(d_{x^2-y^2}\)-projected susceptibility. This phase does not include Sunny
ground-state or excitation calculations.

## Physical Model

The tight-binding basis is

\[
(x_1,z_1,x_2,z_2),
\]

where the subscript labels the layer. The \(z^2\) orbital supplies the local
moment and the \(x^2-y^2\) orbital is the only itinerant RKKY mediator.
All-orbital susceptibility sums must not be used for this extraction.

Use the layer-averaged projected channels

\[
\chi_{\rm intra}^{x}(\mathbf q)
=
\frac{\chi_{11}(\mathbf q)+\chi_{33}(\mathbf q)}{2},
\]

\[
\chi_{\rm inter}^{x}(\mathbf q)
=
\frac{\chi_{13}(\mathbf q)+\chi_{31}(\mathbf q)}{2}.
\]

Layer-exchange symmetry is expected but the two components are averaged for
numerical stability. Their mismatch is recorded as a diagnostic.

## Exchange Definitions

The spin convention to be used later is based on exchange constants extracted
from

\[
J(\mathbf R)=-J_K^2\chi(\mathbf R).
\]

Using the Fourier transform

\[
\chi(\mathbf R)
=
\frac{1}{N_q^2}
\sum_{\mathbf q}
e^{+i\mathbf q\cdot\mathbf R}\chi(\mathbf q),
\]

define

\[
J_0=-J_K^2\chi_{\rm inter}^{x}(0,0),
\]

\[
J_1=-J_K^2\chi_{\rm intra}^{x}(1,0),
\qquad
J_2=-J_K^2\chi_{\rm intra}^{x}(1,1),
\]

\[
J_3=-J_K^2\chi_{\rm intra}^{x}(2,0),
\qquad
J'_1=-J_K^2\chi_{\rm inter}^{x}(1,0).
\]

\(J'_1\) represents the four equivalent interlayer displaced nearest-neighbor
bonds \((\pm1,0)\) and \((0,\pm1)\). Only one representative direction is
calculated. The same \(C_4\) and inversion assumptions apply to \(J_1\),
\(J_2\), and \(J_3\). The extracted values are per-bond exchange constants;
equivalent directions are not summed into them.

## Units And Parameters

All energies use eV.

- \(N_k=N_q=100\)
- filling \(n=1.5\)
- temperature \(T=0.01\ {\rm eV}\)
- \(U=0.0:0.1:6.0\ {\rm eV}\)
- \(J_K=1.0\ {\rm eV}\)

The susceptibility has units \({\rm eV}^{-1}\). Save both the bare
coefficients

\[
\widetilde J=-\chi
\]

in \({\rm eV}^{-1}\), and physical exchanges

\[
J=J_K^2\widetilde J
\]

in eV. They are numerically equal for the current \(J_K\), but their meanings
and units remain distinct.

## Computation

For each \(U\):

1. Find the chemical potential for filling \(n=1.5\).
2. Build the cached \(100\times100\) eigensystem.
3. Traverse the full \(q\) grid once.
4. Compute the existing \(4\times4\) susceptibility matrix at each \(q\).
5. Form the projected intra- and interlayer \(x^2-y^2\) channels.
6. Accumulate all five Fourier components in the same loop.
7. Divide by \(N_q^2\), retain the real part, and record imaginary residuals.

This initial implementation deliberately reuses the verified full
susceptibility routine. A later optimization may calculate only matrix
elements \(11,33,13,31\), but only after the projected implementation is
validated.

## Interface

Provide:

- a reusable function that computes the five exchanges for one \(U\);
- a scan function that accepts an arbitrary \(U\) collection;
- an executable script with the confirmed default scan.

The new implementation must not reuse the all-orbital extraction in
`run_total_rkky.jl`. Keep old scripts as historical comparison tools.

Suggested script:

```text
scripts/susceptibility/run_x2_bilayer_rkky.jl
```

Suggested result:

```text
data/results/rkky_x2_bilayer_U0-6_Nk100.jld2
```

## JLD2 Schema

The single output file stores metadata and all completed results.

Metadata includes:

- basis convention;
- \(N_k,N_q,n,T,J_K\);
- Fourier sign convention;
- requested \(U\) list;
- channel definitions and units.

Results include:

- completed `U_values`;
- chemical potentials and actual fillings;
- bare `J0`, `J1`, `J2`, `J3`, `J1p`;
- physical `J0`, `J1`, `J2`, `J3`, `J1p`;
- layer-symmetry residuals;
- Fourier imaginary residuals;
- filling errors.

## Incremental Save And Resume

Use one JLD2 file only. After each completed \(U\):

1. construct a result containing all completed points;
2. write it to a temporary file in the destination directory;
3. replace the formal output only after the temporary write succeeds.

When the script starts and the result exists:

1. load and validate the schema and metadata;
2. require exact agreement of basis, grid, filling, temperature, \(J_K\),
   Fourier convention, and requested \(U\) list;
3. require all completed arrays to have the same length and finite values;
4. continue from the first uncomputed \(U\);
5. do no numerical work when all requested points are complete.

Parameter mismatch or malformed partial results must raise an error instead of
mixing incompatible data.

## Validation

Tests use a small grid and verify:

- basis indices select \(1\) and \(3\) for the itinerant orbitals;
- projected intra/inter channels use the stated layer averages;
- all five exchanges match an independent direct Fourier sum;
- the \(J_K^2\) scaling and eV units are represented correctly;
- layer and imaginary residuals are finite and small for a symmetric case;
- one-file resume skips completed \(U\) values;
- incompatible metadata is rejected;
- the formal result is replaced only after a temporary save succeeds.

The tight-binding Hamiltonian and susceptibility formula are outside the scope
of this change.
