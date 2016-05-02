function S = clustering_entropy(p_ij, m)
%CLUSTERING_ENTROPY Computes the clustering entropy of a normalized agreement matrix as described in Gfeller et al. 2005
%
%
% Inputs:   p_ij, the agreement matrix, (edges are probabilities)
%           m, the number of edges of the original graph
%
% If you use this program, please cite:
%
% Gfeller et al. "Finding instabilities in the community structure of
% complex networks". Phys.Rev.E 72.056135 (2005)
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).

p = nonzeros(p_ij);
if m == 0
    S = 0;
else
    S = -1.0/m.*sum( p.*log2(p+eps)+(1.0-p).*log2(1.0-p+eps));
end