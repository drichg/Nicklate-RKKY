using Sunny
using LinearAlgebra

# ----------------------------
# Parameters
# ----------------------------
units = Units(:meV, :angstrom)

J1 = 1.0      # meV, replace by your fitted value
J2 = 0.0
J3 = 5.0
Jp = 10.0     # interlayer exchange
S  = 1/2

# Lattice constants
a = 3.85
c = 20.0      # effective bilayer unit-cell c, can be arbitrary for 2D plot

# ----------------------------
# Crystal: bilayer square lattice
# basis:
#   site 1 = A layer
#   site 2 = B layer
# ----------------------------
latvecs = lattice_vectors(a, a, c, 90, 90, 90)

positions = [
	[0.0, 0.0, 0.25],   # A layer
	[0.0, 0.0, 0.75],   # B layer
]

types = ["NiA", "NiB"]

cryst = Crystal(latvecs, positions; types)

# ----------------------------
# Spin system
# supercell size for classical minimization
# Need large enough cell to allow q = (π/2, π/2), i.e. period 4
# ----------------------------
sys = System(
	cryst,
	[1 => Moment(s = S, g = 2), 2 => Moment(s = S, g = 2)],
	:dipole;
	units,
)
view_crystal(cryst)
# ----------------------------
# Set exchange couplings
# Bond vectors are in lattice-coordinate units.
# Sunny bond specification can vary slightly by version;
# if this syntax errors, use print_symmetry_table(cryst, max_dist)
# or inspect available bonds first.
# ----------------------------

# Intralayer J1: (1,0), (0,1)
# set_exchange!(sys, J1, Bond(1, 1, [1, 0, 0]))

# set_exchange!(sys, J1, Bond(2, 2, [1, 0, 0]))


# # Intralayer J2: (1,1), (1,-1)
# set_exchange!(sys, J2, Bond(1, 1, [1, 1, 0]))

# set_exchange!(sys, J2, Bond(2, 2, [1, 1, 0]))


# # Intralayer J3: (2,0), (0,2)
# set_exchange!(sys, J3, Bond(1, 1, [2, 0, 0]))

# set_exchange!(sys, J3, Bond(2, 2, [2, 0, 0]))


# # Interlayer vertical Jp: A-B in same in-plane unit cell
# set_exchange!(sys, Jp, Bond(1, 2, [0, 0, 0]))
