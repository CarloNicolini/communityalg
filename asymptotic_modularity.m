function q = asymptotic_modularity(W,ci)
%ASYMPTOTIC_MODULARITY

%ASYMPTOTIC_MODULARITY      Compute modularity of a vertex partition on a binary network.
%
%
%   Inputs      W,  undirected weighted or unweighted network.
%               ci, membership vector
%
%   Outputs:    Q,  value of Asymptotic Modularity (natural logs).
%               pars, partition parameters. [intraclusteredges,
%               intracluster_pairs, number of edges, number of pairs]
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).
%

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
    ks = sum(k(nodes));
    prob = (ks/(2*m))^2;
    q = q + KL(wc/m, prob);
end


function D = KL(q,p)
if (q==p)
    D=0;
    return;
end

D = 0.0;
if (q > 0.0 && p > 0.0)
    D = D + q*log10(q/p);
end

if (q < 1.0 && p < 1.0)
    D = D + (1.0-q)*log10((1.0-q)/(1.0-p));
end
