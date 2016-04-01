function d = JS_similarity(A1,A2)
% Returns the Jensen-Shannon similarity between the laplacian of two graphs
% De Domenico et al, Mapping multiplex hubs in human functional brain network (2016)
% The similarity of two layers can be calculated in terms of differences in their entropy. Given two rescaled Laplacian matrices L[α] and L[β], it is possible to quantify to which extent layer α is different from layer β by their Kullback-Liebler divergence

L1 = quantum_density(A1);
L2 = quantum_density(A2);

d = 0.5*(KL_similarity(L1,(L1+L2)/2) + KL_similarity(L2,(L1+L2)/2));
%vonneumann_entropy(mu)-0.5*(vonneumann_entropy(rho)+vonneumann_entropy(sigma))