function S = consensus_entropy(p_ij)
%CONSENSUS_ENTROPY      Return the entropy of a probabilistic matrix, where
%                       p_ij is the probabilty that node i and node j are
%                       in the same community.
%
%   If you use this program, please cite:
%
%   Karrer, B., Levina E., and M. E. J. Newman. 
%   “Robustness of Community Structure in Networks.” 
%   Physical Review E 77 (4) (April): 046119.
%   doi:10.1103/PhysRevE.77.046119.
%   http://link.aps.org/doi/10.1103/PhysRevE.77.046119.
%
%   In this implementation we don't divide for the number of edges of the
%   network.
%
%   Input:
%           p_ij, probability matrix.
%   Output:
%           S, entropy as -(sum pij*log(pij)+((1-pij)*log(1-pij)))
%
% P_ij is the consensus matrix among |C| communities of a nxn adjacency
% matrix as obtained with the consensus_clustering function.

pp = p_ij(:);
% add eps to keep logarithm argument positive
S = -sum(xlogx(pp) + xlogx(1-pp) );

function y = xlogx(x)
y = x.*log(x);
y(x==0) = 0;
