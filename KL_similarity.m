function DKL = KL_similarity(A1, A2)
% Returns the Kullback-Leibler similarity between the laplacian of two graphs
% De Domenico et al, Mapping multiplex hubs in human functional brain network (2016)
% The similarity of two layers can be calculated in terms of differences in their entropy. Given two rescaled Laplacian matrices L[α] and L[β], it is possible to quantify to which extent layer α is different from layer β by their Kullback-Liebler divergence

L1 = graph_laplacian(A1);
L2 = graph_laplacian(A2);

DKL = sum(diag(L1.*(log(L1)-log(L2))));