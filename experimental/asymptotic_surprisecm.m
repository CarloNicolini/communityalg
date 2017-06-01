function Scm = asymptotic_surprisecm(W, ci)
%ASYMPTOTIC_SURPRISECM      Compute asymptotical surprise of a vertex partition on a binary network.
%
%
%   Inputs:     A,  undirected weighted or unweighted network.
%               ci, membership vector
%
%   Outputs:    S,  value of Asymptotic Surprise with configuration model correction (base10 logs).
%
%   If you use this program, please cite:
%   V.A.Traag, R.Aldecoa. 2015. “Detecting Communities Using Asymptotical Surprise.” 
%   http://arxiv.org/abs/1503.0044.
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).

% Get the block matrix
[B,C,K,n,m,p,Bnorm,nc] = comm_mat(W,ci);
mc = diag(B);
lzeta = sum((K./(2*m)).^2);
% Use Kullback Leibler divergence to compute asymptotic surprise
Scm = m*KL(sum(mc)/m, lzeta);
