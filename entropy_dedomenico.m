function S = entropy_dedomenico(A,beta)
% Compute laplacian
%L=graph_laplacian(A);
L = diag(sum(A)) - A;
% Compute spectral decomposition of laplacian
LambdaL=eig(L);

% Compute density matrix rho (eq 3.)
expminusbetaL = expm(-beta*L);
rho = expminusbetaL/(trace(expminusbetaL));

% Compute entropy S of density matrix rho in log base two (bits)
% LambdaS = eig(rho);
% Srho = -sum(diag(lambdaS).*log2(lambdaS));

% Compute partition function (eq 4)
Z = sum(exp(-LambdaL.*beta));

% Compute entropy of graph (eq 7)
S = log2(Z) + beta*trace(L*rho);
