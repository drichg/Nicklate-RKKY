# Project Reorganization Design

## Goal

Turn the current mixed research folder into a maintainable Julia project that
keeps susceptibility calculations working and leaves explicit locations for
future ground-state and excitation modules.

## Directory Structure

```text
src/
  NickelateRKKY.jl
  models/
    nickelate_bilayer.jl
  susceptibility/
    susceptibility.jl
  ground_state/
    README.md
  excitations/
    spin_excitation.jl
scripts/
  susceptibility/
  plotting/
data/
  model/
  rkky/
  results/
results/
  susceptibility/
  rkky/
paper/
references/
tests/
docs/plans/
```

`src/NickelateRKKY.jl` is the stable loader for the model and susceptibility
code. The Sunny-based excitation draft remains separate because Sunny is not
currently declared in `Project.toml`; loading the main module must not require
optional future dependencies.

## Code Boundaries

- `src/models/nickelate_bilayer.jl` owns the tight-binding Hamiltonian,
  quasiparticle renormalization, and model-table lookup.
- `src/susceptibility/susceptibility.jl` owns Fermi functions, filling,
  chemical-potential search, cached eigensystems, and susceptibility
  transforms.
- `src/excitations/spin_excitation.jl` stores the current Sunny experiment.
- `src/ground_state/README.md` documents the intended location of future
  ground-state implementations.
- Executable calculations and plot generation live under `scripts/`.

## Data And Output Policy

- Immutable or reusable inputs go under `data/model/` and `data/rkky/`.
- JLD2 calculation results go under `data/results/`.
- Generated PNG files go under `results/susceptibility/` or `results/rkky/`.
- TeX sources, final PDFs, and presentations go under `paper/`.
- Downloaded papers and extracted text go under `references/`.

## Cleanup Policy

Delete reproducible LaTeX build products (`.aux`, `.log`, `.fls`,
`.fdb_latexmk`, `.out`, `.synctex.gz`, `.xdv`), `Manifest.toml.bak`, and the
empty `rslt/` directory. Preserve source files, input data, JLD2 data, final
PDFs, presentations, and PNG results.

## Compatibility And Verification

All file lookup is anchored to `@__DIR__` or a project-root constant instead
of the caller's working directory. Tests load `src/NickelateRKKY.jl` and plot
scripts by absolute paths. Existing numerical formulas are not intentionally
changed. Existing implementation/test disagreements are reported separately
rather than silently changing scientific definitions during folder cleanup.

