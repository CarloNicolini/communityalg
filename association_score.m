function [ass, eij] = association_score(A,membership)
%ASSOCIATION_SCORE      Computes the association score between communities
%                       in the network A.
%
%   If you use this program, please cite:
%
%   Ruan, Jianhua, and Weixiong Zhang. 2008.
%   “Identifying Network Communities with a High Resolution.”
%   Physical Review E 77 (1) (January 14): 016104.
%   doi:10.1103/PhysRevE.77.016104.
%   http://link.aps.org/doi/10.1103/PhysRevE.77.016104.
%
%
%   The theoretical p-value can be estimated analytically with a hypergeometric distribution with parameters M, ai, and aj.
%   This relationship can be best explained by considering  how  the  network  randomization  works.
%   Conceptually, to randomly rewire a network, we break each edge into two halves, or stubs, and then randomly  reconnect  these  stubs
%   into edges. Therefore, for a given network with m edges, we have a box of
%   M=2*m edge stubs.
%   The probability to observe exactly k edges between community ci and cj is equivalent to the probability of randomly drawing ai stubs from the box and observing k of them connected to the vertices in cj. This
%   probability is given by:
%
%   binomial(aj,k)*binomial(M-aj,ai-k)/binomial(M,ai)
%
%   the p-value for observing at least eij edges between ci and cj, therefor
%   can be computed by:
%
%   P = sum(binomial(aj,k)*binomial(M-aj,ai-k)/binomial(M,ai), {k,eij,min(ai,aj)})
%
%   Given the p-value of observing some number of edges between a pair of communities, we compute an association score between the two communities
%
%   S(ci,cj) = -log10(P)
%
%   Note that is also defined for i=j, which can be used to define the statistical significance of a community.
%   We define two communities as associated if their association score between
%   two communities is greater than 2 (p<0.01) and affiliated if greater than 1 (p>0.1).
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2015).
%

ncomms = length(unique(membership));
if max(membership) ~= ncomms
    error('Non continuous membership vector. Must reindex membership.');
end

ass = zeros(ncomms,ncomms);
M = 2*number_of_edges(A);

eij = zeros(ncomms,ncomms);

for i=1:ncomms
    for j=1:ncomms
        B = A(i==membership,j==membership);
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
