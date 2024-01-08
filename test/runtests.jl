using Test

using DataDrop
using SparseArrays
using LinearAlgebra
using GEPHelpers: gep_smallest

function __load_pencil(b)
    K = DataDrop.retrieve_matrix(b, "/K")
    M = DataDrop.retrieve_matrix(b, "/M")
    K, M
end

function __load_frequencies(b)
    DataDrop.retrieve_matrix(b, "/frequencies")
end

b = "unit_cube_modes-h20-n1=3"

K, M = __load_pencil(b)
fs = __load_frequencies(b)
neigvs = length(fs)

d, v, nconv = gep_smallest(K, M, neigvs)
@test nconv == neigvs
fs1 = sqrt.(abs.(d)) ./ (2*pi)
@test norm(fs1 - fs) / norm(fs) <= 1e-5

b = "unit_cube_modes-h8-n1=3"

K, M = __load_pencil(b)
fs = __load_frequencies(b)
neigvs = length(fs)

d, v, nconv = gep_smallest(K, M, neigvs)
@test nconv == neigvs
fs1 = sqrt.(abs.(d)) ./ (2*pi)
@test norm(fs1 - fs) / norm(fs) <= 1e-5

true
