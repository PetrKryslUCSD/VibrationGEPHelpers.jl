let

    omega_shift = 2 * pi * 0.1

    for m in [:Arpack, :SubSIt, :KrylovKit, :ArnoldiMethod]
        @info "Method $(m)"


        for b  in  ["unit_cube_modes-h20-n1=3", "unit_cube_modes-h8-n1=3", "unit_cube_tet-16"]
            @info "Input $b"

            K, M = __load_pencil(b)
            @show fs = __load_frequencies(b)

            for neigvs in [3, 6, 9, 12, 17, 20]
                @info "Number of eigenvalues $(neigvs)"

                d, v, nconv = gep_smallest(K + omega_shift^2 * M, M, neigvs; method = m)
                d .-= omega_shift^2

                @test length(d) == size(v, 2)
                @test length(d) == neigvs

                @show fs1 = sqrt.(abs.(d)) ./ (2*pi)

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
