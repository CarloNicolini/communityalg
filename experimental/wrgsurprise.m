function S = wrgsurprise(W,ci)
%weighted random graph surprise
%
%
%   Inputs      W,  undirected weighted network.
%               ci, membership vector
%
%   Outputs:    Sd
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).
%

[B,C,K,n,m,p,Bnorm,nc] = comm_mat(W,ci);

win = sum(diag(B)); % number of intracluster edges
W=sum(sum(triu(B,1)));
Min = sum(nc.*(nc-1)/2); % number of intracluster pairs
M = p;

[win W Min M]
S=compute_surprise(M+W-1, Min+win-1, W, win);
