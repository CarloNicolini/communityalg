function Zq = zmodularity(W,ci)
%ZMODULARITY      Compute the zmodularity of a vertex partition on a binary network.
% 				  based on the work of Miyauchi and Kawase
%				  Z-Score-Based Modularity for Community Detection in Networks, Plos One
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
[B,C,Kc,~,m] = comm_mat(W,ci);

P = Kc/(2*m); % degree matrix
% zmodularity definition as in 
Zq=sum(diag(B)./m - P.^2)/sqrt(sum(P.^2).*(1-sum(P.^2))); 
