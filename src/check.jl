"""
    check_M_orthogonality(v, M)

Check the mass-orthogonality of the eigenvectors.

# Returns

- `max_vMv_diag_error`, `max_vMv_offdiag_error`, `max_diag_ij`, `max_i`,
  `max_j`, `Mred`: absolute deviations of the diagonal entries of the reduced
  mass matrix from unity, absolute deviations of the off-diagonal entries of the
  reduced mass matrix from zero, and the row/column diagonal index, the row and
  column index of the off-diagonal entry with the largest error. The last output
  argument is the reduced matrix itself (which in this case is expected to be
  identity).
"""
function check_M_orthogonality(v, M)
    max_vMv_diag_error = 0.0
    max_vMv_offdiag_error = 0.0
    max_diag_ij, max_i, max_j = 0, 0, 0
    Mred = v' * M * v
    for i in axes(Mred, 1), j = i:size(Mred, 2)
        p = (Mred[i, j] + Mred[j, i]) / 2
        if i == j
            if abs(p - 1) > max_vMv_diag_error
                max_vMv_diag_error = abs(p - 1)
                max_diag_ij = i
            end
        else
            if abs(p) > max_vMv_offdiag_error
                max_vMv_offdiag_error = abs(p)
                max_i, max_j = i, j
            end
        end
    end
    return max_vMv_diag_error, max_vMv_offdiag_error, max_diag_ij, max_i, max_j, Mred
end

"""
    check_K_orthogonality(d, v, K)

Check the stiffness-orthogonality of the eigenvectors.

# Returns

- `max_vKv_diag_error`, `max_vKv_offdiag_error`, `max_diag_ij`, `max_i`,
  `max_j`, `Kred`: absolute deviations of the diagonal entries of the reduced stiffness
  matrix from the eigenvalue squared, absolute deviations of the off-diagonal
  entries of the reduced mass matrix from zero, and the row/column diagonal
  index, the row and column index of the off-diagonal entry with the largest
  error. The last output argument is the reduced matrix itself. In this case it
  is expected to be a diagonal matrix with the numbers `d` (squared angular
  frequencies) on the diagonal.

"""
function check_K_orthogonality(d, v, K)
    max_vKv_diag_error = 0.0
    max_vKv_offdiag_error = 0.0
    max_diag_ij, max_i, max_j = 0, 0, 0
    Kred = v' * K * v
    for i in eachindex(d), j = i:length(d)
        p = (Kred[i, j] + Kred[j, i]) / 2
        if i == j
            if abs(p - d[i]) > max_vKv_diag_error
                max_vKv_diag_error = abs(p - d[i])
                max_diag_ij = i
            end
        else
            if abs(p) > max_vKv_offdiag_error
                max_vKv_offdiag_error = abs(p)
                max_i, max_j = i, j
            end
        end
    end
    return max_vKv_diag_error, max_vKv_offdiag_error, max_diag_ij, max_i, max_j, Kred
end
