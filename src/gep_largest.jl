"""
    gep_largest(K, M, maxit = 30, rtol = 1/10000)

Estimate largest eigenvalue with power iteration.

Find the largest ``\\omega`` such that
```math
(\\mathbf{K} - \\omega^2 \\mathbf{M})\\mathbf{v} = \\mathbf{0}
```
Here ``\\mathbf{K}`` is a symmetric stiffness matrix, `` \\mathbf{M}`` is a
symmetric positive definite mass matrix, ``\\omega`` is the angular velocity.

!!! note
    The mass matrix `M` is assumed to be diagonal.
"""
function gep_largest(K, M, maxit = 30, rtol = 1/10000)
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
