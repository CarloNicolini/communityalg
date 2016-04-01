function S=rwalkent(A,beta)
% Random walker entropy in a graph as in 
% Estrada et al. Walk entropies in graphs, "Linear Algebra and its Applications 443 (2014) 235–244"

% The parameter beta is the inverse of temperature
eba = expm(beta.*A);

Z = trace(eba); % partition function

pib = diag(eba)/Z;

% random walk entropy of a graph 
S = -sum(diag(eba)/Z.*(diag(log(eba))-log(Z)));
% renormalize it in [0,1] as entropy of a random walker with beta->0 and beta->infinity is
% exactly log(n)
S = S/log(length(A));
