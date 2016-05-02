function [S] = surprise(A, ci)
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
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).
%

if length(unique(A(:))) ~= 2
    warning on;
    warning('Input matrix is not binary {0,1}. Ignoring edge weights to compute Surprise.');
end

% Get the block matrix
[B,C,~,~,m,p]=comm_mat(A,ci);

nc = sum(C,2); % number of nodes per community

mc = sum(diag(B)); % number of intracluster edges
pc = sum(nc.*(nc-1)/2); % number of intracluster pairs

S=compute_surprise(p, pc, m, mc);