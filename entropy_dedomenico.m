function S = entropy_dedomenico(A,beta)
% Compute laplacian
%L=graph_laplacian(A);
L = diag(sum(A)) - A;
% Compute spectral decomposition of laplacian
lambdas=eig(L);

% Compute density matrix rho (eq 3.)
rho = expm(-beta*L);
rho = rho/trace(rho);

% Compute entropy S of density matrix rho in log base two (bits)
lambdas_rho = eig(rho);
Srho = -sum(diag(lambdas_rho).*log2(lambdas_rho));

% Compute partition function (eq 4)
Z = sum(exp(-lambdas.*beta));

% Compute entropy of graph (eq 7)
S = log2(Z) + beta*trace(L*rho);
