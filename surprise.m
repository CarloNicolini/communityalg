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

B = double(A~=0);
n = length(B);
m = number_of_edges(B);
p = n*(n-1)/2;
groups = membership2groups(ci); % convert membership vector to groups
ncomms = length(groups); % number of communities
intraedges = 0; % number of intracluster edges
intrapairs = 0; % number of intracluster pairs
for i=1:ncomms
    nodes = groups{i};
    % generate subgraph
    g = B(nodes,nodes);
    ec = sum(g(:))/2;
    nc = length(g);
    pc = nc*(nc-1)/2;
    intraedges = intraedges + ec;
    intrapairs = intrapairs + pc;
end
% Clustering parameters.
pars = [intraedges,intrapairs,m,p];

S=compute_surprise(p, intrapairs, m, intraedges);