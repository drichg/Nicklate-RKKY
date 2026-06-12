# Add Pi-Over-2 Point Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add `(pi/2, pi/2)` as an explicit high-symmetry tick in `plot_chi_highsym.jl`.

**Architecture:** Keep the cached q-index implementation. Require `Nk` to be divisible by 4 so `(pi/2, pi/2)` lies exactly at q-index `(Nk/4, Nk/4)`.

**Tech Stack:** Julia, Test, Plots, existing `susceptibility.jl`.

---

### Task 1: Test Path Tick

Modify `tests/test_plot_chi_highsym.jl` so `high_symmetry_path(4)` expects labels `Gamma, Y, M, P, Gamma` and ticks at `(0,0)`, `(0,pi)`, `(pi,pi)`, `(pi/2,pi/2)`, `(0,0)`.

### Task 2: Implement Path

Modify `plot_chi_highsym.jl` to require `Nk % 4 == 0`, add `P` as a tick when q-index equals `(Nk/4, Nk/4)`, and keep the user's bottom parameters.

### Task 3: Verify

Run:

```powershell
julia --project=. tests/test_plot_chi_highsym.jl
julia --project=. tests/test_mu.jl
```
