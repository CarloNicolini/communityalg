function [d,L1,L2] = jensen_shannon_sim(L1,L2)
% Returns the Jensen-Shannon similarity between the laplacian of two graphs
% De Domenico et al, Mapping multiplex hubs in human functional brain network (2016)
% The similarity of two layers can be calculated in terms of differences in their entropy. Given two rescaled Laplacian matrices L[α] and L[β], it is possible to quantify to which extent layer α is different from layer β by their Kullback-Liebler divergence
mu = (L1+L2)/2;
d = KL_similarity(L1,mu)/2 + KL_similarity(L2,mu)/2;