In matlab one can do:

I = h5read('unit_cube_modes-h8-n1=3.h5', '/K/I');
J = h5read('unit_cube_modes-h8-n1=3.h5', '/K/J');
V = h5read('unit_cube_modes-h8-n1=3.h5', '/K/V');
m = max(I);
n = max(J);
K = sparse(I, J, V, m, n);

I = h5read('unit_cube_modes-h8-n1=3.h5', '/M/I');
J = h5read('unit_cube_modes-h8-n1=3.h5', '/M/J');
V = h5read('unit_cube_modes-h8-n1=3.h5', '/M/V');
m = max(I);
n = max(J);
M = sparse(I, J, V, m, n);

tic;
[v, d] = eigs(K + 0.1 * M, M, 20, 'SM');
toc
d = diag(diag(d) - 0.1);
norm(K * v - M * v * d)

f = h5read('unit_cube_modes-h8-n1=3.h5', '/frequencies/matrix');
