function q = modularity(W,ci)
%MODULARITY      Compute modularity of a vertex partition on a binary network.
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
[B,C,K,~,m] = comm_mat(W,ci);

P = K/(2*m); % degree matrix
Mc = diag(B); % intracluster weights matrix
q=sum(Mc/m-P.^2); % modularity as difference
