function [sj,sp] = cluster_similarity(e1, e2)
    sj = zeros([length(e2),length(e1)]);
    sp = zeros([length(e2),length(e1)]);
    % Count the total number of nodes in the network
    N =  sum(cellfun(@length,e1));
    % Use the p-value of the cluster homogeneity:
    %http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099755#pone.0099755.e106
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
