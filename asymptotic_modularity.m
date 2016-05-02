function q = asymptotic_modularity(W,ci)
%ASYMPTOTIC_MODULARITY      Compute the asymptotic modularity of a vertex partition on a binary network.
%
%   Inputs:      W,  undirected weighted or unweighted network.
%               ci, membership vector
%
%   Outputs:    Q,  value of Asymptotic Modularity (natural logs).
%               pars, partition parameters. [intraclusteredges,
%               intracluster_pairs, number of edges, number of pairs]
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).

n = length(W);
m = sum(nonzeros(triu(W)));
ncomms = length(unique(ci)); % number of communities

q=0; % asympt modularity value
k=degrees_und(W);
for c=1:ncomms
    nodes = find(ci==c);
    g = W(nodes,nodes);
    wc = sum(nonzeros(triu(g))); % intracomm edges
    % compute configuration model
    ks = sum(k(nodes))
    prob = (ks/(2*m))^2;
    q = q + KL(wc/m, prob);
end