function S  = vonneumann_entropy(A)
% Implementation from Structural reducibility of multilayer networks.
% De Domenico, Manlio Nicosia, Vincenzo Arenas, Alexandre Latora, Vito 
% Nature Communications 6, 6864, (2015)

L = quantum_density(A);
S = -trace(L.*log(L));
% Another way to compute it to observe that the Von Neumann Entropy of a
% density matrix corresponds to the Shannon Entropy of its power spectrum
lambda = eig(L); % compute the spectrum
lambda(find(lambda<=0))=eps; % fix for numerical precision
S2 = - sum(lambda.*log(lambda));