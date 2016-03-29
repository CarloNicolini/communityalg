function L=graph_laplacian(A)
% Compute the laplacian
n = length(A);
D = eye(n);
D(1:n+1:n*n) = degrees_und(A);
L = D-A;