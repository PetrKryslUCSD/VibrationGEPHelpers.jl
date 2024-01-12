let

    omega_shift = 2 * pi * 0.1

    b = "unit_cube_modes-h20-n1=3"
    @info "Input $b"
    K, M = __load_pencil(b)
    fs = __load_frequencies(b)
    neigvs = length(fs)

    d, v, nconv = gep_smallest(K + omega_shift^2 * M, M, neigvs)
    d .-= omega_shift^2
    @test nconv == neigvs
    fs1 = sqrt.(abs.(d)) ./ (2*pi)
    @test norm(fs1 - fs) / norm(fs) <= 1e-5

    r = check_K_orthogonality(d, v, K)
    @test norm(r, Inf) <=  1e-10
    r = check_M_orthogonality(v, M)
    @test norm(r, Inf) <=  1e-10


    b = "unit_cube_modes-h8-n1=3"
    @info "Input $b"
    K, M = __load_pencil(b)
    fs = __load_frequencies(b)
    neigvs = length(fs)

    d, v, nconv = gep_smallest(K + omega_shift^2 * M, M, neigvs)
    d .-= omega_shift^2
    @test nconv == neigvs
    fs1 = sqrt.(abs.(d)) ./ (2*pi)
    @test norm(fs1 - fs) / norm(fs) <= 1e-5

    r = check_K_orthogonality(d, v, K)
    @test norm(r, Inf) <=  1e-10
    r = check_M_orthogonality(v, M)
    @test norm(r, Inf) <=  1e-10

    @time d, v, nconv = gep_smallest(K + omega_shift^2 * M, M, neigvs)

    b = "unit_cube_tet-16"
    @info "Input $b"
    if !isfile(joinpath(dirname(@__FILE__()), b * ".h5"))
        success(run(`unzip -qq -d $(dirname(@__FILE__())) $(joinpath(dirname(@__FILE__()), b * ".zip"))`; wait = false))
    end
    K, M = __load_pencil(b)
    fs = __load_frequencies(b)
    neigvs = length(fs)

    d, v, nconv = gep_smallest(K + omega_shift^2 * M, M, neigvs)
    d .-= omega_shift^2
    @test nconv == neigvs
    fs1 = sqrt.(abs.(d)) ./ (2*pi)
    @test norm(fs1 - fs) / norm(fs) <= 1e-5

    r = check_K_orthogonality(d, v, K)
    @test norm(r, Inf) <=  1e-10
    r = check_M_orthogonality(v, M)
    @test norm(r, Inf) <=  1e-10

    @time d, v, nconv = gep_smallest(K + omega_shift^2 * M, M, neigvs)
end
