function [S, pars] = asymptotic_surprise_old(W, ci)
%ASYMPTOTIC_SURPRISE

%ASYMPTOTIC_SURPRISE      Compute surprise of a vertex partition on a binary network.
%
%   If you use this program, please cite:
%   V.A.Traag, R.Aldecoa. 2015. “Detecting Communities Using Asymptotical Surprise.” 
%   http://arxiv.org/abs/1503.0044.
%
%   Inputs      A,  undirected weighted or unweighted network.
%               ci, membership vector
%
%   Outputs:    S,  value of Asymptotic Surprise (base10 logs).
%               pars, partition parameters. [intraclusteredges,
%               intracluster_pairs, number of edges, number of pairs]
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2015).
%

n = length(W);
m = sum(sum(triu(W)));
p = n*(n-1)/2;
groups = membership2groups(ci); % convert membership vector to groups
ncomms = length(groups); % number of communities
intraedges = 0; % number of intracluster edges
intrapairs = 0; % number of intracluster pairs
for i=1:ncomms
    nodes = groups{i};
    % generate subgraph
    g = W(nodes,nodes);
    wc = sum(nonzeros(triu(g)));
    nc = length(nodes);
    pc = nc*(nc-1)/2;
    %fprintf('%d\t%.2f\t%d\t%d\n',i,wc,nc,pc);
    intraedges = intraedges + wc;
    intrapairs = intrapairs + pc;
end
% Clustering parameters.
pars = [intraedges,intrapairs,m,p]
% Use Kullback Leibler divergence to compute asymptotic surprise
S = m*KL(intraedges/m, intrapairs/p);
