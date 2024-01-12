# GEPHelpers.jl

Generalized Eigenvalue Problem helper functions.

Targeted at matrix pencils with real symmetric matrices. The GEP reads

(K - omega^2 * M) * v = 0

Here K is a symmetric, positive semi-definite, stiffness matrix, M is a
symmetric positive definite mass matrix, omega is the angular velocity.

