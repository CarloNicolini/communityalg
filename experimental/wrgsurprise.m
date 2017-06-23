function S = wrgsurprise(W,ci,base10)
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
if nargin==2
	base10=true;
end

[B,C,K,n,m,p,Bnorm,nc] = comm_mat(W,ci);

win = sum(diag(B)); % number of intracluster edges
W=sum(sum(triu(W,1)));
Min = sum(nc.*(nc-1)/2); % number of intracluster pairs
M = p;

[win W Min M]
% It's logbase10!!!
S=compute_surprise(M+W-1, Min+win-1, W, win,base10);
