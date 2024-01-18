# VibrationGEPHelpers.jl

Vibration Generalized Eigenvalue Problem (GEP) helper functions.

Targeted at matrix pencils with real symmetric matrices. The GEP reads

(K - omega^2 * M) * v = 0

Here K is a symmetric, positive semi-definite, stiffness matrix, M is a
symmetric positive definite mass matrix, omega is the angular velocity.

# Usage

- Solve the vibration GEP for the smallest eigenvalues:
```
d, v, nconv = gep_smallest(K + omega_shift^2 * M, M, neigvs; method = :Arpack)
```
The solution is useful in constructing modal expansions in solid dynamics.

- Solve the vibration GEP for the largest eigenvalue:
```
omega_max = gep_largest(K, M)
```
The solution is useful in estimating the stable time step in explicit
integration of the equations of motion.


