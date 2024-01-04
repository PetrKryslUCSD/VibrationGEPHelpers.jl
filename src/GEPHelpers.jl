module GEPHelpers

"""
    pwr_largest(K, M)

Estimate largest eigenvalue with power iteration.

!!! note
    The mass matrix `M` is assumed to be diagonal.
"""
function pwr_largest(K, M)
    invM = fill(0.0, size(M, 1))
    invM .= 1.0 ./ (vec(diag(M)))
    v = rand(size(M, 1))
    w = fill(0.0, size(M, 1))
    for i in 1:40
        mul!(w, K, v)
        wn = norm(w)
        w .*= (1.0/wn)
        v .= invM .* w
        vn = norm(v)
        v .*= (1.0/vn)
    end
    # ev = sqrt((v' * K * v) / (v' * M * v)) # EXTREMELY SLOW
    mul!(w, K, v) # FAST
    ev = sqrt((w' * v) / (v' * (M * v))) # FAST
    return ev
end

end # module GEPHelpers
