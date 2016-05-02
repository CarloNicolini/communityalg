function [S, pars] = asymptotic_surprise(W, ci)
%ASYMPTOTIC_SURPRISE      Compute asymptotical surprise of a vertex partition on a binary network.
%
%
%   Inputs:     A,  undirected weighted or unweighted network.
%               ci, membership vector
%
%   Outputs:    S,  value of Asymptotic Surprise (base10 logs).
%               pars, partition parameters. [intraclusteredges,
%               intracluster_pairs, number of edges, number of pairs]
%
%   If you use this program, please cite:
%   V.A.Traag, R.Aldecoa. 2015. “Detecting Communities Using Asymptotical Surprise.” 
%   http://arxiv.org/abs/1503.0044.
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).

% Get the block matrix
[B,C]=comm_mat(W,ci);

n=length(W); % number of nodes
m=sum(sum(triu(W))); % number of edges 
p=n*(n-1)/2; % number of pairs
nc = sum(C,2);

mc = sum(diag(B)); % number of intracluster edges
pc = sum(nc.*(nc-1)/2); % number of intracluster pairs

% Use Kullback Leibler divergence to compute asymptotic surprise
S = m*KL(mc/m, pc/p);
