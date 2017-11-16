function [avgphi, allphi, allpvals, avgpval] = phicoeffavg(all_memberships)

% Implementation of the phi coefficient as in
% "Network community structure alterations in adult schizophrenia:
% "Identification and localization of alterations"
%
% Lerman-Sinkoff, Dov B. Barch, Deanna M.
% Neuroimage, Clinical (2016)

% all_memberships is a matrix of integers taking values from 1 to maxc
% where the number of rows is the number of independent repetitions of a
% community detection method, while the number of columns is the number of
% nodes in the graph

numnodes = size(all_memberships,2); % number of nodes
numsubjs = size(all_memberships,1); % number of subjects
delta = nan(numnodes,numnodes,numsubjs);
corr_test='right'; % tests that the correlation is greater than zero

for subj=1:numsubjs
    CijA = double(bsxfun(@eq,all_memberships(subj,:), unique(all_memberships(subj,:))'));
    delta(:,:,subj) = CijA'*CijA;
end

% Now computes the phi coefficient
allphi = nan(numnodes,numsubjs,numsubjs);
allpvals = nan(numnodes,numsubjs,numsubjs);
for node=1:numnodes
    for subjA=1:numsubjs
        for subjB=subjA+1:numsubjs %exclude diagonals
            [phival, pval] = corr(delta(node,:,subjA)',delta(node,:,subjB)','tail',corr_test);
            allphi(node,subjA,subjB) = phival;
            allpvals(node,subjA,subjB) = pval;
        end
    end
end

avgphi=nan(1,numnodes);

for node=1:numnodes
    % Select the upper triangular part of the phi ND array with the phi
    % coefficients of each node for every pair of subjects
    % put nan in the lower triangular and then ignores them
    phivals = triu(squeeze(allphi(node,:,:)),1) + tril(nan(numsubjs));
    phivals = phivals(~isnan(phivals));
    %phivals = triu(squeeze(phi(node,:,:)),1);
    % Compute the average over the pairs of subjects of the phi for each
    % node
    zphi = atanh(phivals-eps); % to avoid inf and nans around
    zphi = zphi(~isinf(zphi));
    avgphi(node) = tanh(mean(zphi));
end
avgpval = nanmean(nanmean(allpvals,3),2);
avgphi = nanmean(nanmean(allphi,3),2)';

end