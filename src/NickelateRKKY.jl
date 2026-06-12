module NickelateRKKY

using DelimitedFiles
using LinearAlgebra

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const MODEL_DATA_DIR = joinpath(PROJECT_ROOT, "data", "model")
const Z_DATA_PATH = joinpath(MODEL_DATA_DIR, "z.txt")

include("models/nickelate_bilayer.jl")
include("susceptibility/susceptibility.jl")

export PROJECT_ROOT, MODEL_DATA_DIR, Z_DATA_PATH
export z_data_row, tight_binding, renormalize, chemical_potential
export fermi, fermi_derivative, susceptibility_kernel
export band_energies_kmesh, filling_number, filling_number_per_orbital
export find_chemical_potential, chi_matrix_q, diagonalized_kmesh
export chi_matrix_q_index, chi_matrix_R

end
