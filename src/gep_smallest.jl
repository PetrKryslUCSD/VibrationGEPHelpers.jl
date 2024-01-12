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
function gep_smallest(K, M, neigvs; method = :Arpack, orthogonalize = false,
    which=:SM, tol=0.0, maxiter=300, sigma=nothing, ritzvec=true, v0=zeros((0,))
    )

    @assert which == :SM
    @assert sigma == nothing
    @assert ritzvec == true

    if method == :Arpack
        d, v, nconv = eigs(Symmetric(K), Symmetric(M); nev=neigvs, which=:SM, tol=tol, maxiter=maxiter, explicittransform=:none, check = 1)
    elseif method == :ArnoldiMethod
        d, v, nconv = arnoldimethod_eigs(Symmetric(K), Symmetric(M); nev=neigvs, which=:SM, tol=tol, maxiter=maxiter, explicittransform=:none, check = 1)
    elseif method == :KrylovKit
        d, v, nconv = krylovkit_eigs(Symmetric(K), Symmetric(M); nev=neigvs, which=:SR, tol=tol, maxiter=maxiter)
    elseif method == :SubSIt
        d, v, nconv = subsit_eigs(Symmetric(K), Symmetric(M); nev=neigvs, which=:SM, tol=tol, maxiter=maxiter)
    else
        error("Unknown method: $(method)")
    end

    orthogonalize && __mass_orthogonalize!(v, M)

    return d, v, nconv
end


struct ShiftAndInvert{TA, TB, TT}
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
    a = ShiftAndInvert(cholesky(A), B, Vector{eltype(A)}(undef, size(A,1)))
    LinearMap{eltype(A)}(a, size(A,1), ismutating=true)
end


function arnoldimethod_eigs(K, M; nev::Integer=6, ncv::Integer=max(20,2*nev+1), which=:SM, tol=0.0, maxiter::Integer=300, sigma=nothing, v0::Vector=zeros(eltype(K),(0,)), ritzvec::Bool=true, explicittransform::Symbol=:auto, check::Integer=0)
    which == :SM || error("Argument which: The only recognized which is :SM")
    sigma == nothing || error("Argument sigma not supported")
    ritzvec == true || error("Argument ritzvec not supported")
    explicittransform == :none || error("Argument explicittransform only supported as :none")

    decomp, history = partialschur(construct_linear_map(K, M), nev=nev, tol=tol, restarts=maxiter, which=LM())
    d_inv, v = partialeigen(decomp)
    d = 1 ./ d_inv
    d = real.(d)
    ix = sortperm(real.(d))

    return d[ix], v[:, ix], history.nconverged
end

function krylovkit_eigs(K, M; nev::Integer=6, ncv::Integer=max(20,2*nev+1), which=:SR, tol=0.0, maxiter::Integer=300)
    which == :SR || error("Argument which: The only recognized which is :SR")
    d, vv, convinfo = geneigsolve((Symmetric(K), Symmetric(M)), nev, :SR; maxiter=maxiter, issymmetric = true, ishermitian = true, isposdef = true)
    v = zeros(size(K, 1), length(d))
    for j in 1:length(d)
        v[:, j] .= vv[j]
    end
    @show convinfo
    nconv = convinfo.converged
    d, v, nconv
end

function subsit_eigs(K, M; nev::Integer=6, ncv::Integer=max(20,2*nev+1), which=:SM, tol=0.0, maxiter::Integer=300)
    which == :SM || error("Argument which: The only recognized which is :SM")
    d, v, nconv = SubSIt.ssit(K, M; nev=nev, maxiter=maxiter, tol = tol)
end
