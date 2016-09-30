function L=graph_laplacian(A)
% Compute the laplacian
n = length(A);
D = diag(sum(A));
L = D-A;