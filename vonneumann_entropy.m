function S  = vonneumann_entropy(A)
L = quantum_density(A);

S = -trace(L.*log2(L));
lambda = eig(L);
S2 = - sum(lambda.*log2(lambda));

[S,S2]