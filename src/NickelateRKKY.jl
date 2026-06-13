module NickelateRKKY

using DelimitedFiles
using JLD2
using LinearAlgebra

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const MODEL_DATA_DIR = joinpath(PROJECT_ROOT, "data", "model")
const Z_DATA_PATH = joinpath(MODEL_DATA_DIR, "z.txt")

include("models/nickelate_bilayer.jl")
include("susceptibility/susceptibility.jl")
include("susceptibility/rkky_x2_bilayer.jl")

export PROJECT_ROOT, MODEL_DATA_DIR, Z_DATA_PATH
export z_data_row, tight_binding, renormalize, chemical_potential
export fermi, fermi_derivative, susceptibility_kernel
export band_energies_kmesh, filling_number, filling_number_per_orbital
export find_chemical_potential, chi_matrix_q, diagonalized_kmesh
export chi_matrix_q_index, chi_matrix_R
export X1_INDEX, X2_INDEX, x2_layer_channels
export compute_x2_bilayer_rkky, x2_rkky_metadata, validate_x2_rkky_resume
export save_x2_rkky_result_atomic, load_x2_rkky_result
export scan_x2_bilayer_rkky

end
