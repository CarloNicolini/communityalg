function [ass, eij] = association_score(A,ci)
%ASSOCIATION_SCORE      Computes the association score of communities in the network A.
%
%   Inputs:     A, the adjacency matrix
%               ci, the membership vector
%   Outputs:    ass, the association score
%               eij, a square block matrix with the number of links from
%               community i to community j.
%
%   The p-value for observing at least eij edges between community i and community j can be computed by:
%   P = sum(binomial(aj,k)*binomial(M-aj,ai-k)/binomial(M,ai), {k,eij,min(ai,aj)})
%   Given the p-value of observing some number of edges between a pair of communities, we compute an association 
%   score between the two communities:
%   S(ci,cj) = -log10(P)
%
%   Note that is also defined for i=j, which can be used to define the statistical significance of a community.
%   We define two communities as associated if their association score between
%   two communities is greater than 2 (p<0.01) and affiliated if greater than 1 (p>0.1).
%
%   If you use this program, please cite:
%   Ruan, Jianhua, and Weixiong Zhang. 2008.
%   "Identifying Network Communities with a High Resolution."
%   Phys. Rev. E 77 (1) (January 14): 016104.
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).

ncomms = length(unique(ci));
if max(ci) ~= ncomms
    error('Non continuous membership vector. Must reindex membership.');
end

ass = zeros(ncomms,ncomms);
M = 2*number_of_edges(A);

eij = zeros(ncomms,ncomms);

for i=1:ncomms
    for j=1:ncomms
        B = A(i==ci,j==ci);
        eij(i,j) = sum(sum(B));
    end
end

% Divide diagonal by two, because communities i-i are considered twice.
eij(1:ncomms+1:ncomms*ncomms) = eij(1:ncomms+1:ncomms*ncomms)/2;

%  iterate on all communities
for i=1:ncomms
    ai = eij(i,i);
    for j=1:ncomms
        % comm j intradegree
        aj = eij(j,j);
        % links between community i and j
        ass(i,j) = compute_surprise(M, aj, ai, eij(i,j),true);
        if ass(i,j)<0
            ass(i,j)=0;
            %warning('Error computing HypergeometricAssociationScore');
            %fprintf(2,'i=%d j=%d p=%d pi=%d m=%d mi=%d S=%g\tM=%d aj=%d ai=%d eij=%d\n',i,j,M,aj,ai,eij(i,j),ass(i,j),M,aj,ai,eij(i,j));
        end
    end
end
