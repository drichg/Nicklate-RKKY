# Bilayer x2 RKKY Extraction Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a tested Julia workflow that extracts \(J_0,J_1,J_2,J_3,J'_1\) from the layer-averaged \(x^2-y^2\) susceptibility for one \(U\) or a resumable \(U\) scan.

**Architecture:** Put reusable projection, Fourier accumulation, metadata validation, and resume logic in a new source module loaded by `NickelateRKKY.jl`. Keep the executable default scan in a thin script, and update one JLD2 result atomically after each completed \(U\).

**Tech Stack:** Julia 1.11, LinearAlgebra, JLD2, Test, existing cached susceptibility implementation.

---

### Task 1: Define The x2 Layer Projection

**Files:**
- Create: `src/susceptibility/rkky_x2_bilayer.jl`
- Modify: `src/NickelateRKKY.jl`
- Create: `tests/test_rkky_x2_bilayer.jl`

**Step 1: Write failing projection tests**

Create a test matrix with distinct entries and assert:

```julia
chi = reshape(ComplexF64.(1:16), 4, 4)
channels = x2_layer_channels(chi)

@test X1_INDEX == 1
@test X2_INDEX == 3
@test channels.intra == (chi[1, 1] + chi[3, 3]) / 2
@test channels.inter == (chi[1, 3] + chi[3, 1]) / 2
@test channels.intra_layer_residual == abs(chi[1, 1] - chi[3, 3])
@test channels.inter_layer_residual == abs(chi[1, 3] - chi[3, 1])
```

**Step 2: Run the test to verify RED**

Run:

```powershell
julia --project=. tests/test_rkky_x2_bilayer.jl
```

Expected: fail because `x2_layer_channels` and the named indices do not exist.

**Step 3: Add the minimal projection implementation**

Implement named constants and:

```julia
function x2_layer_channels(chi::AbstractMatrix)
	size(chi) == (4, 4) || throw(DimensionMismatch("chi must be 4x4"))

	chi11 = chi[X1_INDEX, X1_INDEX]
	chi22 = chi[X2_INDEX, X2_INDEX]
	chi12 = chi[X1_INDEX, X2_INDEX]
	chi21 = chi[X2_INDEX, X1_INDEX]

	return (
		intra = (chi11 + chi22) / 2,
		inter = (chi12 + chi21) / 2,
		intra_layer_residual = abs(chi11 - chi22),
		inter_layer_residual = abs(chi12 - chi21),
	)
end
```

Include and export the new file from `src/NickelateRKKY.jl`.

**Step 4: Run the focused test**

Expected: projection test passes.

**Step 5: Commit**

```powershell
git add src/NickelateRKKY.jl src/susceptibility/rkky_x2_bilayer.jl tests/test_rkky_x2_bilayer.jl
git commit -m "feat: add bilayer x2 susceptibility projection"
```

### Task 2: Extract Five Exchanges For One U

**Files:**
- Modify: `src/susceptibility/rkky_x2_bilayer.jl`
- Modify: `tests/test_rkky_x2_bilayer.jl`

**Step 1: Write a failing direct-sum comparison**

On `Nk=4`, calculate a reference by explicitly traversing \(q\), calling
`chi_matrix_q_index`, applying `x2_layer_channels`, and separately summing:

```julia
phase_J1 = exp(im * qx)
phase_J2 = exp(im * (qx + qy))
phase_J3 = exp(im * 2qx)

reference = (
	bare_J0 = -real(inter_R0),
	bare_J1 = -real(intra_R10),
	bare_J2 = -real(intra_R11),
	bare_J3 = -real(intra_R20),
	bare_J1p = -real(inter_R10),
)
```

Assert that:

```julia
actual = compute_x2_bilayer_rkky(0.0; Nk=4, Nq=4, temperature=0.01, JK=2.0)
```

matches the independent reference, and that every physical `J*` equals
`JK^2 * bare_J*`.

**Step 2: Run the test to verify RED**

Expected: fail because `compute_x2_bilayer_rkky` does not exist.

**Step 3: Implement one-pass Fourier accumulation**

The function must:

- require `Nk == Nq`;
- find `mu` for the requested filling;
- create one cached eigensystem;
- visit every \(q\) once;
- accumulate `inter R=(0,0)`, `intra R=(1,0)`, `intra R=(1,1)`,
  `intra R=(2,0)`, and `inter R=(1,0)`;
- divide all sums by `Nq^2`;
- return bare and physical exchanges, `mu`, actual filling, filling error,
  maximum layer residuals, and maximum Fourier imaginary residual.

Use field names:

```julia
bare_J0, bare_J1, bare_J2, bare_J3, bare_J1p
J0, J1, J2, J3, J1p
```

**Step 4: Run focused and regression tests**

Run:

```powershell
julia --project=. tests/test_rkky_x2_bilayer.jl
julia --project=. tests/test_mu.jl
julia --project=. tests/test_run2.jl
```

Expected: all tests pass.

**Step 5: Commit**

```powershell
git add src/susceptibility/rkky_x2_bilayer.jl tests/test_rkky_x2_bilayer.jl
git commit -m "feat: compute five bilayer x2 RKKY exchanges"
```

### Task 3: Define The Result Schema And Metadata Validation

**Files:**
- Modify: `src/susceptibility/rkky_x2_bilayer.jl`
- Modify: `tests/test_rkky_x2_bilayer.jl`

**Step 1: Write failing schema tests**

Test a metadata constructor with:

```julia
U_values = collect(0.0:0.1:0.2)
metadata = x2_rkky_metadata(
	U_values;
	Nk=4,
	Nq=4,
	n_target=1.5,
	temperature=0.01,
	JK=1.0,
)
```

Assert the basis, units, Fourier convention, channel descriptions, and exact
requested \(U\) list. Test that `validate_x2_rkky_resume` accepts matching
metadata and rejects a changed `Nk`, `temperature`, `JK`, or \(U\) list.

**Step 2: Run the test to verify RED**

Expected: fail because the metadata helpers do not exist.

**Step 3: Implement schema helpers**

Use a schema version field and explicit strings:

```julia
schema_version = 1
basis = "x1,z1,x2,z2"
fourier_convention = "chi(R)=sum_q exp(+i q.R) chi(q)/Nq^2"
susceptibility_unit = "eV^-1"
exchange_unit = "eV"
```

Do not use approximate metadata matching. Resume only when all scalar values
and the complete requested \(U\) list match.

**Step 4: Run the focused test**

Expected: schema and mismatch tests pass.

**Step 5: Commit**

```powershell
git add src/susceptibility/rkky_x2_bilayer.jl tests/test_rkky_x2_bilayer.jl
git commit -m "feat: define resumable x2 RKKY result schema"
```

### Task 4: Add Atomic Save And Resume

**Files:**
- Modify: `src/susceptibility/rkky_x2_bilayer.jl`
- Modify: `tests/test_rkky_x2_bilayer.jl`

**Step 1: Write failing persistence tests**

Using `mktempdir()`:

1. Save a one-point partial result.
2. Load it and validate every array length.
3. Run a scan over two \(U\) values with an injected lightweight point
   calculator that records calls.
4. Confirm only the missing point is calculated.
5. Confirm incompatible metadata throws before the calculator is called.
6. Confirm no `.tmp` file remains after a successful save.

The point calculator injection is part of the scan API only to keep resume
tests fast and deterministic.

**Step 2: Run the test to verify RED**

Expected: fail because persistence and scan functions do not exist.

**Step 3: Implement atomic persistence**

Implement:

```julia
save_x2_rkky_result_atomic(filepath, result)
load_x2_rkky_result(filepath)
scan_x2_bilayer_rkky(U_values; ..., filepath, point_calculator=compute_x2_bilayer_rkky)
```

Write the temporary file beside the destination. Use one formal JLD2 file and
replace it only after `jldsave` succeeds. Keep all result arrays in one named
tuple under a stable key such as `"result"`.

Validate:

- metadata;
- equal completed-array lengths;
- completed length not exceeding requested length;
- completed `U_values` equal the requested prefix;
- finite numeric entries.

**Step 4: Run the focused test**

Expected: persistence and resume tests pass.

**Step 5: Commit**

```powershell
git add src/susceptibility/rkky_x2_bilayer.jl tests/test_rkky_x2_bilayer.jl
git commit -m "feat: add atomic resumable x2 RKKY scans"
```

### Task 5: Add The Default Executable Scan

**Files:**
- Create: `scripts/susceptibility/run_x2_bilayer_rkky.jl`
- Modify: `tests/test_rkky_x2_bilayer.jl`
- Modify: `README.md`

**Step 1: Write a failing script-interface test**

Include the script and assert it defines a configuration constructor or
constants matching:

```julia
U_values = collect(0.0:0.1:6.0)
Nk = 100
Nq = 100
n_target = 1.5
temperature = 0.01
JK = 1.0
```

Assert the default output path is:

```text
data/results/rkky_x2_bilayer_U0-6_Nk100.jld2
```

Including the script must not start the scan.

**Step 2: Run the test to verify RED**

Expected: fail because the script does not exist.

**Step 3: Implement the thin executable**

The guarded executable must call `scan_x2_bilayer_rkky` with the confirmed
defaults and print:

- metadata at startup;
- completed/resumed point count;
- each completed \(U\), chemical potential, filling, and five exchanges;
- final output path.

Update `README.md` with the command:

```powershell
julia --project=. scripts/susceptibility/run_x2_bilayer_rkky.jl
```

Warn in the README that `Nk=100`, 61 \(U\) points, and the full \(q\) grid
make this a long-running calculation.

**Step 4: Run focused tests and parse the script**

Run:

```powershell
julia --project=. tests/test_rkky_x2_bilayer.jl
julia --project=. -e "Meta.parseall(read(\"scripts/susceptibility/run_x2_bilayer_rkky.jl\", String))"
```

Expected: both commands succeed without launching the default scan.

**Step 5: Commit**

```powershell
git add scripts/susceptibility/run_x2_bilayer_rkky.jl tests/test_rkky_x2_bilayer.jl README.md
git commit -m "feat: add default bilayer x2 RKKY scan"
```

### Task 6: Full Verification

**Files:**
- No production edits expected

**Step 1: Run the new focused test**

```powershell
julia --project=. tests/test_rkky_x2_bilayer.jl
```

Expected: all new tests pass.

**Step 2: Run existing numerical regression tests**

```powershell
julia --project=. tests/test_mu.jl
julia --project=. tests/test_run2.jl
julia --project=. tests/test_plot_rkky_jld2.jl
julia --project=. tests/test_plot_total.jl
julia --project=. tests/test_project_structure.jl
```

Expected: all listed tests pass.

**Step 3: Run the known inconsistent test separately**

```powershell
julia --project=. tests/test_plot_chi_highsym.jl
```

Expected baseline: the two pre-existing orbital projection assertions may
still fail because that plot currently uses diagonal-only values. Do not alter
that physical plotting behavior as part of this task.

**Step 4: Verify repository scope**

```powershell
git status --short
git diff --check
```

Confirm the untracked `rkky_bilayer_nickelate_codex.md` remains untouched
unless the user separately requests that it be committed.

