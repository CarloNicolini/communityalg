function [sj,sp] = membership_similarity(m1, m2)
%MEMBERSHIP_SIMILARITY    Compute the cluster similarity between two membership vectors,
%   Input: m1 vertex membership of the first partition
%          m2 vertex membership of the second partition
%   Output:
%         sj: Jaccard similarty
%         sp: p-value of the cluster homogeneity
%         http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099755#pone.0099755.e106
%         ssta: statistical significance of extracted community structure
%
%   See also cluster_similarity function
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).
    e1 = membership2groups(m1);
    e2 = membership2groups(m2);

    [sj,sp]=cluster_similarity(e1,e2);