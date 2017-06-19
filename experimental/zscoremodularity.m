function Zq = zscoremodularity(W,ci)
%ZMODULARITY      Compute the zscore modularity as in Nicolini, Garlaschelli
%
%
%   Inputs      W,  undirected weighted or unweighted network.
%               ci, membership vector
%
%   Outputs:    Zq,  value of zscore Modularity
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2017).
%
[B,C,Kc,~,m] = comm_mat(W,ci);

P = Kc/(2*m); % degree matrix
% zmodularity definition as in 
Zq=sum(diag(B)./m - P.^2)/sqrt(sum(P.^2).*(1-sum(P.^2))); 
