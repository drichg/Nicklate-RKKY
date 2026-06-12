# Orbital High-Symmetry Susceptibility Design

## Goal

Make `plot_chi_highsym.jl` a dedicated utility for plotting one of two
layer-resolved orbital susceptibility combinations:

- `x2`: `Re(chi[1,1] + chi[1,3])`
- `z2`: `Re(chi[2,2] + chi[2,4])`

## Interface

The user edits an `orbital` variable near the bottom of the script and runs:

```powershell
julia --project=. plot_chi_highsym.jl
```

Accepted values are `"x2"` and `"z2"`. Any other value raises an
`ArgumentError`. Output files use names such as `chi_x2_U0.0.png`.

## Structure

Keep the existing high-symmetry path and cached susceptibility calculation.
Add small helpers that validate the orbital name, return its display label,
and evaluate the requested combination from a susceptibility matrix.

Guard the executable section with:

```julia
if abspath(PROGRAM_FILE) == @__FILE__
```

so tests can include the file without starting the full calculation.

## Testing

Test the orbital selector directly with a known matrix, test invalid inputs,
verify a small `Nk=4` path against `chi_matrix_q_index`, and verify plot output.
