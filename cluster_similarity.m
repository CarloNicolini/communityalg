function [sj,sp] = cluster_similarity(e1, e2)
%CLUSTER_SIMILARITY    Compute the cluster similarity between two clustering as groups
%
%   Input: group1, as compute from membership2groups function
%          group2, as compute from membership2groups function
%
%   Output:
%         sj: Jaccard similarty
%         sp: p-value of the cluster homogeneity
%         http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099755#pone.0099755.e106
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).

sj = zeros([length(e2),length(e1)]);
sp = zeros([length(e2),length(e1)]);
% Count the total number of nodes in the network
N =  sum(cellfun(@length,e1));

for i=1:length(e2)
    c2i = e2{i};
    for j=1:length(e1)
        c1j = e1{j};
        inter_set = length(intersect(c1j,c2i));
        union_set = length(union(c1j,c2i));
        sj(i,j) = inter_set/union_set;
        sp(i,j) = -logHyperProbability(N,length(c1j),length(c2i),inter_set);
    end
end
