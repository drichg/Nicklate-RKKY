# High-Symmetry Susceptibility Plot Design

## Goal

Add a Julia program that plots the real part of a user-selected bare susceptibility component along the path `(0,0) -> (0,pi) -> (pi,pi) -> (0,0)`.

## User-Facing Behavior

The program asks for:

- `U`
- susceptibility component

Supported component inputs:

- `11` or `1,1` for `chi[1,1]`
- any valid pair from `1..4`, for example `22`, `13`, `24`, or `2,4`
- `intra` for `sum(chi[1:2, 1:2])`
- `inter` for `sum(chi[1:2, 3:4])`

The program computes `mu` from `find_chemical_potential(1.5, U)`, evaluates the susceptibility on the requested path, and saves a PNG such as `chi_11_U0.0.png`.

## Architecture

Create `plot_chi_highsym.jl` with small reusable functions:

- build the high-symmetry path and distance axis
- parse component input into a selector
- compute the selected susceptibility along the path
- plot and save the result
- run an interactive prompt only when the file is executed directly

Use the existing `susceptibility.jl` functions and `Plots.jl`.

## Constraints

This directory is not a git repository, so the design document cannot be committed here.
