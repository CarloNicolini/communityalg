function S = laplacian_entropy(A)

% Compute the degree normalized laplacian
L=quantum_density(A);

% eival = eig(L);
% eival(find(abs(eival<=eps)))=0; % set to zero eigenvalues very close to 0
% [eival log2(eival)]
% 
% S = -sum(eival.*log2(eival));
