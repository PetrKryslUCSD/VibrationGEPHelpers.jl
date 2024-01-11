function __mass_orthogonalize!(v, M)
    for i in axes(v, 2)
        v[:, i] /= sqrt(dot(v[:, i], M * v[:, i]))
    end
    return v
end

"""
    gep_smallest(K, M, neigvs; orthogonalize = false,
        which=:SM, tol=0.0, maxiter=300, sigma=nothing, ritzvec=true, v0=zeros((0,))
        )


"""
function gep_smallest(K, M, neigvs; orthogonalize = false,
    which=:SM, tol=0.0, maxiter=300, sigma=nothing, ritzvec=true, v0=zeros((0,))
    )

    @assert which == :SM
    @assert sigma == nothing
    @assert ritzvec == true

    d, v, nconv = eigs(Symmetric(K), Symmetric(M); nev=neigvs, which=:SM, tol=tol, maxiter=maxiter, explicittransform=:none, check = 1) # maxiter = maxiter,

    orthogonalize && __mass_orthogonalize!(v, M)

    return d, v, nconv
end

