function [S,mcd,K] = degree_corrected_surprise(W,ci)
%DEGREE_CORRECTED_SURPRISE      Compute degree corrected surprise of a vertex partition on a binary network.
%
%
%   Inputs      W,  undirected weighted or unweighted network.
%               ci, membership vector
%
%   Outputs:    Q,  value of Degree Corrected Surprise (natural logs).
%               pars, partition parameters. [intraclusteredges,
%               intracluster_pairs, number of edges, number of pairs]
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).
%

[B,C,~,~,m,p]=comm_mat(W,ci);


% The probabilistic definition
