"""
    mass_orthogonalize!(v, M)

Orthogonalize the columns of the matrix `v` with respect to the mass matrix `M`.
"""
function mass_orthogonalize!(v, M)
    for i in axes(v, 2)
        v[:, i] ./= sqrt(dot(view(v, :, i), M * view(v, :, i)))
    end
    return v
end

"""
    gep_smallest(
        K,
        M,
        neigvs;
        method = :Arpack,
        tol = 0.0,
        maxiter = 300,
        v0 = fill(zero(eltype(K)), 0, 0),
    )

Solve for the smallest natural frequencies.

# Return

- `d, v, nconv`: vector of squares of angular frequencies `d`, matrix of
  eigenvectors as columns `v`, how many eigenvalues converged `nconv`.
"""
function gep_smallest(
    K,
    M,
    neigvs;
    method = :Arpack,
    tol = 0.0,
    maxiter = 300,
    v0 = fill(zero(eltype(K)), 0, 0),
)
    method in [:KrylovKit,:Arpack, :SubSIt, :ArnoldiMethod,] || error("Unknown method $(method)")

    if method == :Arpack
        d, v, nconv = eigs(
            Symmetric(K),
            Symmetric(M);
            nev = neigvs,
            which = :SM,
            tol = tol,
            maxiter = maxiter,
            explicittransform = :none,
            check = 1,
        )
    elseif method == :ArnoldiMethod
        d, v, nconv = __arnoldimethod_eigs(
            Symmetric(K),
            Symmetric(M);
            nev = neigvs,
            tol = tol,
            maxiter = maxiter,
        )
    elseif method == :KrylovKit
        d, v, nconv = __krylovkit_eigs(
            Symmetric(K),
            Symmetric(M);
            nev = neigvs,
            tol = tol,
            maxiter = maxiter,
        )
    elseif method == :SubSIt
        d, v, nconv = __subsit_eigs(
            Symmetric(K),
            Symmetric(M);
            nev = neigvs,
            tol = tol,
            maxiter = maxiter,
            X = v0,
        )
    else
        error("Unknown method: $(method)")
    end

    return d, v, nconv
end


struct ShiftAndInvert{TA,TB,TT}
    A_factorization::TA
    B::TB
    temp1::TT
    temp2::TT
end

function (SI::ShiftAndInvert)(y, x)
    SI.temp1 .= x
    SI.temp2 .= SI.A_factorization.PtL' \ SI.temp1
    mul!(SI.temp1, SI.B, SI.temp2)
    y .= (SI.A_factorization.PtL \ SI.temp1)
end

function construct_linear_map(A, B)
    a = ShiftAndInvert(cholesky(A), B, Vector{eltype(A)}(undef, size(A, 1)), Vector{eltype(A)}(undef, size(A, 1)))
    LinearMap{eltype(A)}(a, size(A, 1), ismutating = true)
end

function __arnoldimethod_eigs(
    K,
    M;
    nev::Integer = 6,
    ncv::Integer = max(20, 2 * nev + 1),
    tol = 0.0,
    maxiter::Integer = 300,
    v0::Vector = zeros(eltype(K), (0,)),
)
    tol = (tol == 0 ? sqrt(eps(1.0)) : tol)
    si = ShiftAndInvert(cholesky(K), M, Vector{eltype(K)}(undef, size(K, 1)), Vector{eltype(K)}(undef, size(K, 1)))
    m = LinearMap{eltype(K)}(si, size(K, 1), ismutating = true)
    decomp, history = partialschur(
        m,
        nev = nev,
        tol = tol,
        restarts = maxiter,
        which = LM(),
        mindim  = nev + 6,
        # maxdim = max(nev + 8, 2 * nev, 40)
        maxdim = min(max(40, 2 * nev), size(K,1))
    )
    @show history
    d_inv, v = partialeigen(decomp)
    # Recover the solution of the original problem
    d = 1 ./ real.(d_inv)
    # Sort  the angular frequencies by magnitude.  Make sure all imaginary parts
    # of the eigenvalues are removed.
    ix = sortperm(d)
    d, v = d[ix[1:nev]], real.(v[:, ix[1:nev]])
    v = si.A_factorization.PtL'\v
    # Make vectors mass orthogonal
    mass_orthogonalize!(v, M)
    return d, v, history.nconverged
end

function __krylovkit_eigs(
    K,
    M;
    nev::Integer = 6,
    ncv::Integer = max(20, 2 * nev + 1),
    tol = 0.0,
    maxiter::Integer = 300,
)
    _nev = nev + 6
    z = (zero(eltype(M)))
    # Employ invert strategy to accelerate convergence
    Kfactor = cholesky(Symmetric(K))
    PtL = Kfactor.PtL
    di, ww, convinfo = eigsolve(
            x -> PtL \ (M * (PtL' \ x)),
            rand(typeof(z), size(K, 1)),
            _nev,
            :LR;
            maxiter = maxiter,
            krylovdim = 2 * _nev + 6,
            ishermitian = true
        )
    vv = map(w->PtL'\w, ww)
    # Eigen values of the original problem
    d = 1 ./ real.(di)
    # Convert a vector of vectors to a matrix
    v = zeros(size(K, 1), length(d))
    for j in eachindex(vv)
        v[:, j] .= real.(vv[j])
    end
    mass_orthogonalize!(v, M)
    # Sort  the angular frequencies by magnitude. Carve out only the sorted
    # eigenvalues we are interested in.
    ix = sortperm(d)
    return d[ix[1:nev]], v[:, ix[1:nev]], convinfo.converged
end

function __subsit_eigs(
    K,
    M;
    nev::Integer = 6,
    ncv::Integer = max(20, 2 * nev + 1),
    tol = 0.0,
    maxiter::Integer = 300,
    X = fill(zero(eltype(K)), 0, 0),
)
    d, v, nconv = SubSIt.ssit(K, M; nev = nev, maxiter = maxiter, tol = tol, verbose = false)
    return d, v, nconv
end
