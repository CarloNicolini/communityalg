function [w, pars] = wonder(W, ci)
%WONDER
%WONDER      Compute wonder of a vertex partition on a binary network.
%
%
%   Inputs      A,  undirected weighted or unweighted network.
%               ci, membership vector
%
%   Outputs:    S,  value of Wonder (base 10 logs).
%               pars, partition parameters. [intraclusteredges,
%               intracluster_pairs, number of edges, number of pairs]
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).
%
error('===WARNING==== This implementation is not correct');
n = length(W);
m = sum(nonzeros(triu(W)));
p = n*(n-1)/2;
blockmat = comm_mat(W,ci);

intraedges = sum(diag(blockmat));
intrapairs = sum(sum(blockmat,2).*(sum(blockmat,2)-1)/2);

x = intraedges/m;
y = (intrapairs/p);
mxy = (x+y)/2;

w = (KL(x,mxy) + KL(y,mxy))/2;

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
