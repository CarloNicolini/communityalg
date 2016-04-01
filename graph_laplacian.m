function L=graph_laplacian(A)
% Compute the laplacian
n = length(A);
D = eye(n);
D(1:n+1:n*n) = strengths_und(A);
L = D-A;