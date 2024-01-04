module GEPHelpers

using LinearAlgebra

"""
    pwr_largest(K, M, maxit = 30, rtol = 1/10000)

Estimate largest eigenvalue with power iteration.

!!! note
    The mass matrix `M` is assumed to be diagonal.
"""
function pwr_largest(K, M, maxit = 30, rtol = 1/10000)
    invM = fill(0.0, size(M, 1))
    invM .= 1.0 ./ (vec(diag(M)))
    v = rand(size(M, 1))
    w = fill(0.0, size(M, 1))
    everyn = Int(round(maxit / 50)) + 1
    omega = omegap = 0.0
    for i in 1:maxit
        mul!(w, K, v) # ThreadedSparseCSR.bmul!(w, K, v)
        wn = norm(w)
        w .*= (1.0/wn)
        v .= invM .* w
        vn = norm(v)
        v .*= (1.0/vn)
        if i % everyn  == 0
            omega = sqrt(dot(v, mul!(w, K, v)) / (v' * M * v))
            @show i, abs(omega - omegap) / omega
            if abs(omega - omegap) / omega  < rtol
                break
            end
            omegap = omega
        end
    end
    return omega
end

function _pwr(K, M, maxit = 30, rtol = 1/10000)
    invM = fill(0.0, size(M, 1))
    invM .= 1.0 ./ (vec(diag(M)))
    v = rand(size(M, 1))
    w = fill(0.0, size(M, 1))
    everyn = Int(round(maxit / 50)) + 1
    lambda = lambdap = 0.0
    for i in 1:maxit
        mul!(w, K, v) # ThreadedSparseCSR.bmul!(w, K, v)
        wn = norm(w)
        w .*= (1.0/wn)
        v .= invM .* w
        vn = norm(v)
        v .*= (1.0/vn)
        if i % everyn  == 0
            lambda = sqrt(dot(v, mul!(w, K, v)) / (v' * M * v))
            @show i, abs(lambda - lambdap) / lambda
            if abs(lambda - lambdap) / lambda  < rtol
                break
            end
            lambdap = lambda
        end
    end
    return lambda
end

end # module GEPHelpers
