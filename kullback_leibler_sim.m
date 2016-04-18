function d = kullback_leibler_sim(L1, L2)
% Returns the Kullback-Leibler similarity between the laplacian of two graphs
% De Domenico et al, Mapping multiplex hubs in human functional brain network (2016)
% The similarity of two layers can be calculated in terms of differences in their entropy.
% Given two rescaled Laplacian matrices L[α] and L[β], it is possible to quantify to which 
% extent layer α is different from layer β by their Kullback-Leibler divergence
d = trace(L1*(logm(L1)-logm(L2)));