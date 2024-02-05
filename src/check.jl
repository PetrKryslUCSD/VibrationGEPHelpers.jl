"""
    check_M_orthogonality(v, M)

Check the mass-orthogonality of the eigenvectors.

# Returns

- `max_vMv_diag_error`, `max_vMv_offdiag_error`: absolute deviations of the
  diagonal entries of the reduced mass matrix from unity, absolute deviations
  of the off-diagonal entries of the reduced mass matrix from zero.
"""
function check_M_orthogonality(v, M)
    max_vMv_diag_error = 0.0
    max_vMv_offdiag_error = 0.0
    Mred = v' * M * v
    for i in axes(Mred, 1), j = i:size(Mred, 2)
        p = (Mred[i, j] + Mred[j, i]) / 2
        if i == j
            max_vMv_diag_error = max(max_vMv_diag_error, abs(p - 1))
        else
            max_vMv_offdiag_error = max(max_vMv_offdiag_error, abs(p))
        end
    end
    return max_vMv_diag_error, max_vMv_offdiag_error
end

"""
    check_K_orthogonality(d, v, K)

Check the stiffness-orthogonality of the eigenvectors.

# Returns

- `max_vKv_diag_error`, `max_vKv_offdiag_error`: absolute deviations of the
  diagonal entries of the reduced stiffness matrix from the eigenvalue squared,
  absolute deviations of the off-diagonal entries of the reduced stiffness
  matrix from zero.
"""
function check_K_orthogonality(d, v, K)
    max_vKv_diag_error = 0.0
    max_vKv_offdiag_error = 0.0
    Kred = v' * K * v
    for i in eachindex(d), j = i:length(d)
        p = (Kred[i, j] + Kred[j, i]) / 2
        if i == j
            max_vKv_diag_error = max(max_vKv_diag_error, abs(p - d[i]))
        else
            max_vKv_offdiag_error = max(max_vKv_offdiag_error, abs(p))
        end
    end
    return max_vKv_diag_error, max_vKv_offdiag_error
end
