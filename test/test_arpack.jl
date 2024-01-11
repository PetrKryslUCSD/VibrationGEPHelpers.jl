let

    b = "unit_cube_modes-h20-n1=3"

    K, M = __load_pencil(b)
    fs = __load_frequencies(b)
    neigvs = length(fs)

    d, v, nconv = gep_smallest(K, M, neigvs)
    @test nconv == neigvs
    fs1 = sqrt.(abs.(d)) ./ (2*pi)
    @test norm(fs1 - fs) / norm(fs) <= 1e-5

    r = check_K_orthogonality(d, v, K)
    @test norm(r, Inf) <=  1e-10
    r = check_M_orthogonality(v, M)
    @test norm(r, Inf) <=  1e-10


    b = "unit_cube_modes-h8-n1=3"

    K, M = __load_pencil(b)
    fs = __load_frequencies(b)
    neigvs = length(fs)

    d, v, nconv = gep_smallest(K, M, neigvs)
    @test nconv == neigvs
    fs1 = sqrt.(abs.(d)) ./ (2*pi)
    @test norm(fs1 - fs) / norm(fs) <= 1e-5

    r = check_K_orthogonality(d, v, K)
    @test norm(r, Inf) <=  1e-10
    r = check_M_orthogonality(v, M)
    @test norm(r, Inf) <=  1e-10

end
