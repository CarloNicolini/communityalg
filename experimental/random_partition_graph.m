% Build a random modular graph, given number of modules, and link densities.
%
% INPUTs: number of nodes (n), number of modules (c), total link density (p),
%         and ratio of nodal degree to nodes within the same module
%         to the degree to nodes in other modules (r);
%         if specified, "labels" overwrites the random node
%         assignment to clusters, eg: labels = [1,1,2,3,3,4]
% OUTPUTs: adjacency matrix, modules to which the nodes are assigned
%
% Idea and code about pre-specified labels by Jonathan Hadida, June 12, 2014
% GB: last updated, July 6, 2014

function [adj, membership] = random_partition_graph(n,p_in,p_out, membership)

% n - number of nodes
% p_in - intracluster probability of attachment
% p_out - intercluster probability of attachment
% membership - pre-specified cluster assignments

assert( length(membership) == n )

% DERIVATION of probabilities
% k_in/k_out = r, k_in + k_out = k = p(n-1)
% => (1/r)k_in + k_in = p(n-1), => k_in = rp(n-1)/(r+1), k_out = p(n-1)/(r+1)
% k_in = p_in*(n/c-1) => p_in = rpc(n-1)/((r+1)(n-c))
% k_out = p_out*(n-n/c) => p_out = pc(n-1)/(n(r+1)(c-1))

adj = double(rand(n) < p_out); % already set the outer links (the majority, and fast!)
% then set the inner links
for c=unique(membership)
    nc = sum(membership==c);
    nodes = membership==c;
    adj(nodes,nodes) = double(rand(nc) < p_in);
end
adj(1:n+1:n^2)=0;
