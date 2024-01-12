using Test

using DataDrop
using SparseArrays
using LinearAlgebra
using GEPHelpers: gep_smallest, check_M_orthogonality, check_K_orthogonality

function __load_pencil(b)
    K = DataDrop.retrieve_matrix(b, "/K")
    M = DataDrop.retrieve_matrix(b, "/M")
    K, M
end

function __load_frequencies(b)
    DataDrop.retrieve_matrix(b, "/frequencies")
end

b = "unit_cube_tet-16"
if !isfile(joinpath(dirname(@__FILE__()), b * ".h5"))
    success(run(`unzip -qq -d $(dirname(@__FILE__())) $(joinpath(dirname(@__FILE__()), "matrix_files.zip"))`; wait = false))
end

@time @testset "Arpack" begin
include("test_arpack.jl")
end

@time @testset "SubSIt" begin
include("test_subsit.jl")
end

@time @testset "ArnoldiMethod" begin
include("test_arnoldimethod.jl")
end

@time @testset "KrylovKit" begin
include("test_krylovkit.jl")
end

true
