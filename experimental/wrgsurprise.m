function S = wrgsurprise(W,ci,base10)
%weighted random graph surprise
%
%
%   Inputs      W,  undirected weighted network.
%               ci, membership vector
%
%   Outputs:    S the Bose Einstein Surprise
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).
%
if nargin==2
	base10=true;
end

[B,C,K,N,W,M,Bnorm,nc] = comm_mat(W,ci);

Win = sum(diag(B)); % number of intracluster edges

Min = sum(nc.*(nc-1)/2); % number of intracluster pairs

%[Win W Min M]
S=compute_surprise(M+W-1, Min+Win-1, W, Win,base10); %Checked with Mathematica
% compute_surprise(M+W, Min+Win-1, W, Win,base10) %this is the version which is a negative hypergeometric

% Compute the binomial approximation, too
% This is only for the dominant term in the hypergeometric sum
disp( -arrayfun(@(i)(logC(i,W)+i*log((i+Min-1)/(W+M-1))+(W-i)*(log(1-(i+Min-1)/(W+M-1)))),Win) )
disp( W*KL(Win/W, (Win+Min)/(W+M)) )