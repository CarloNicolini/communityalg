function [S,pars] = surprise(A, ci)
%SURPRISE      Compute surprise of a vertex partition on a binary network.
%
%
%   Inputs      A,  undirected binary network. If weighted, weights are
%                   ignored.
%               ci, membership vector
%
%   Outputs:    S,  value of Surprise (is log10 base, to get natural logarithm,
%                   multiply S by log(10)).
%
%               pars, partition parameters. [intraclusteredges,
%               intracluster_pairs, number of edges, number of pairs]
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2015).
%

if length(unique(A(:))) ~= 2
    warning on;
    warning('Input matrix is not binary {0,1}. Ignoring edge weights to compute Surprise.');
end

[mc, pc, m, p] = partition_params(A,ci);

S=compute_surprise(p, sum(pc), m, sum(mc));