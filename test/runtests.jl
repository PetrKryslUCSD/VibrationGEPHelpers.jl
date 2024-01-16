using Test

using DataDrop
using SparseArrays
using Statistics
using LinearAlgebra
using VibrationGEPHelpers: gep_smallest, check_M_orthogonality, check_K_orthogonality

function __load_pencil(b)
    K = DataDrop.retrieve_matrix(b, "/K")
    M = DataDrop.retrieve_matrix(b, "/M")
    K, M
end

function __load_frequencies(b)
    DataDrop.retrieve_matrix(b, "/frequencies")
end

orthogonality_tol = 1.0e-9
frequency_tol = 1.0e-6
residual_tol = 1.0e-6

b = "unit_cube_tet-16"
if !isfile(joinpath(dirname(@__FILE__()), b * ".h5"))
    success(run(`unzip -qq -d $(dirname(@__FILE__())) $(joinpath(dirname(@__FILE__()), "matrix_files.zip"))`; wait = false))
end


include("test_methods.jl")


true
