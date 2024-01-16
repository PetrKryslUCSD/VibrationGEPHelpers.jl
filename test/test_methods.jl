let

    omega_shift = 2 * pi * 0.1

    for m in [:KrylovKit,:Arpack, :SubSIt]# :ArnoldiMethod,
    # for m in [:Arpack, :SubSIt, :KrylovKit, :ArnoldiMethod]
        @info "Method $(m)"


        for b  in  [
            "unit_cube_modes-h20-n1=3", # [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.264197, 0.264197, 0.359590, 0.35959, 0.359590, 0.3634, 0.363401, 0.36340, 0.40699, 0.408478, 0.4084, 0.46259, 0.46259, 0.462592]
            "unit_cube_modes-h8-n1=3", # [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2895884, 0.2895884, 0.403155, 0.4031553, 0.4031553, 0.4031921, 0.403192, 0.4031921, 0.450131, 0.4501314, 0.50321, 0.570422, 0.67125, 0.671251]
            "unit_cube_tet-16" # [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.259392, 0.260642, 0.351863, 0.3524297, 0.353351, 0.3574280, 0.3575518, 0.3579181, 0.399573, 0.405491, 0.405635, 0.4529674, 0.453472, 0.4543342]
            ]
            @info "Input $b"

            K, M = __load_pencil(b)
            fs = __load_frequencies(b)

            for neigvs in [3, 6, 9, 12, 17, 20]
                @info "Number of eigenvalues $(neigvs)"

                d, v, nconv = gep_smallest(K + omega_shift^2 * M, M, neigvs; method = m)
                d .-= omega_shift^2

                @test length(d) == size(v, 2)
                @test length(d) == neigvs

                fs1 = sqrt.(abs.(d)) ./ (2*pi)

                @test nconv >= neigvs

                @test norm(fs1 - fs[1:length(fs1)]) / norm(fs) <= frequency_tol

                r = check_K_orthogonality(d, v, K)
                @test norm(r, Inf) <= max(10000 * eps(1.0), orthogonality_tol * max(maximum(d), 1.0))
                r = check_M_orthogonality(v, M)
                @test norm(r, Inf) <= orthogonality_tol
            end
        end
    end

    nothing
end
