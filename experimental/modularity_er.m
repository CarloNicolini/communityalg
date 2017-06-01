function qer = modularity_er(W,ci)
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
[B,C,K,n,m,p,Bnorm,nc] = comm_mat(W,ci);

mc = diag(B);
pc = nc.*(nc-1)/2;
rho = 2*m/(n*(n-1));

qer=1/m * sum(mc - rho*pc ); % modularity as difference
