function __mass_orthogonalize!(v, M)
    for i in axes(v, 2)
        v[:, i] ./= sqrt(dot(view(v, :, i), M * view(v, :, i)))
    end
    return v
end

"""
    gep_smallest(K, M, neigvs; orthogonalize = false,
        which=:SM, tol=0.0, maxiter=300, sigma=nothing, ritzvec=true, v0=zeros((0,))
        )


"""
function gep_smallest(
    K,
    M,
    neigvs;
    method = :Arpack,
    orthogonalize = false,
    tol = 0.0,
    maxiter = 300,
    v0 = fill(zero(eltype(K)), 0, 0),
)
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

    orthogonalize && __mass_orthogonalize!(v, M)

    return d, v, nconv
end


struct ShiftAndInvert{TA,TB,TT}
    A_factorization::TA
    B::TB
    temp::TT
end

function (M::ShiftAndInvert)(y, x)
    mul!(M.temp, M.B, x)
    # ldiv!(y, M.A_factorization, M.temp)
    y .= M.A_factorization \ M.temp
end

function construct_linear_map(A, B)
    a = ShiftAndInvert(cholesky(A), B, Vector{eltype(A)}(undef, size(A, 1)))
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
    check::Integer = 0,
)
    decomp, history = partialschur(
        construct_linear_map(K, M),
        nev = nev,
        tol = tol,
        restarts = maxiter,
        which = LM(),
        mindim  = nev + 6,
        maxdim = max(nev + 8, 2 * nev)
    )
    d_inv, v = partialeigen(decomp)
    # Invert the eigenvalues
    d = 1 ./ d_inv
    d = real.(d)
    # Order by absolute value
    ix = sortperm(d)
    d, v = d[ix], v[:, ix]
    @show d
    # Compute a set of M-orthogonal vectors
    Mhat = v' * M * v
    R = eigvecs(Mhat)
    __mass_orthogonalize!(R, Mhat)
    v = real.(v * R)
    # Recalculate the eigenvalues so that the order of the eigenvalues and
    # eigenvectors agrees
    d = zeros(size(d))
    for j in 1:nev
        d[j] = view(v, :, j)'  * (K * view(v, :, j))
    end
    @show sort(d)
    # Sort  the angular frequencies by magnitude.  Make sure all imaginary parts
    # of the eigenvalues are removed.
    ix = sortperm(d)
    return d[ix[1:nev]], v[:, ix[1:nev]], history.nconverged
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
    # Employ invert strategy to accelerate convergence
    z = (zero(eltype(M)))
    Kfactor = cholesky(Symmetric(K))
    di, vv, convinfo = eigsolve(
        x -> Kfactor \ (M * x),
        rand(typeof(z), size(K, 1)),
        _nev,
        :LR;
        maxiter = maxiter,
        krylovdim = 2 * _nev + 6
    )
    # Eigen values of the original problem
    d = real.(1 ./ di)
    # Convert a vector of vectors to a matrix
    v = zeros(size(K, 1), length(d))
    for j in 1:length(vv)
        v[:, j] .= real.(vv[j])
    end
    # Order by absolute value
    ix = sortperm(d)
    d, v = d[ix], v[:, ix]
    @show d
    # Compute a set of M-orthogonal vectors
    Mhat = v' * M * v
    R = eigvecs(Mhat)
    __mass_orthogonalize!(R, Mhat)
    v = real.(v * R)
    # Recalculate the eigenvalues so that the order of the eigenvalues and
    # eigenvectors agrees
    d = zeros(size(v, 2))
    for j in 1:size(v, 2)
        d[j] = view(v, :, j)'  * (K * view(v, :, j))
    end
    @show d
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
    d, v, nconv = SubSIt.ssit(K, M; nev = nev, maxiter = maxiter, tol = tol, verbose = true)
    return d, v, nconv
end
