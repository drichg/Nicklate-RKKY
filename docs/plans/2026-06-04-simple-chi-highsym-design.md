# Simple High-Symmetry Susceptibility Plot Design

## Goal

Replace the previous interactive/component-parser script with a minimal Julia script controlled by variables near the top of the file.

## User-Facing Behavior

The user edits:

```julia
U = 2.0
component = 22
```

Then runs:

```powershell
julia --project=. plot_chi_highsym.jl
```

The script interprets `component = 22` as `chi[2,2]`, `component = 12` as `chi[1,2]`, computes `Re chi[alpha,beta](q)` along `(0,0) -> (0,pi) -> (pi,pi) -> (0,0)`, and saves a PNG.

## Shape

Keep only small helpers for path generation, `component` to `(alpha,beta)`, data calculation, and plotting. Do not support `intra`, `inter`, or interactive input.
