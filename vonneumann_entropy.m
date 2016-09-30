function S  = vonneumann_entropy(A)
% Implementation from Structural reducibility of multilayer networks.
% De Domenico, Manlio Nicosia, Vincenzo Arenas, Alexandre Latora, Vito 
% Nature Communications 6, 6864, (2015)
L = quantum_density(A);
% The Von Neumann entropy is -\Tr[ L \log(L) ] where the product is
% intended to be a matrix product and the logarithm is intended to be a
% matrix logarithm
% S = -trace(L*logm(L));
% Another way to compute it to observe that the Von Neumann Entropy of a
% density matrix corresponds to the Shannon Entropy of its power spectrum
lambda = eig(L); % compute the spectrum
% keep only positive eigenvalues
lambda = lambda(lambda>0);
S = - sum(lambda.*log2(lambda));