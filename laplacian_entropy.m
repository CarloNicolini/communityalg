function S = laplacian_entropy(A)

% Compute the laplacian
L=laplacian(A);
eival = eig(L);
S = -sum(eival.*log(eival));
