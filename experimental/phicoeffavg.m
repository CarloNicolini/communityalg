function ret = phicoeffavg(all_memberships)

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

for subj=1:numsubjs
    CijA = double(bsxfun(@eq,all_memberships(subj,:), unique(all_memberships(subj,:))'));
    delta(:,:,subj) = CijA'*CijA;
end

% Now computes the phi coefficient
ret.phi = nan(numnodes,numsubjs,numsubjs);
ret.corr_pvalues = nan(numnodes,numsubjs,numsubjs);
for node=1:numnodes
    for subjA=1:numsubjs
        for subjB=1:numsubjs
            [phival, pval] = corr(delta(node,:,subjA)',delta(node,:,subjB)');
            [nvi, nmi] = partition_distance(delta(node,:,subjA)'+1,delta(node,:,subjB)'+1);
            ret.phi(node,subjA,subjB) = phival;
            ret.nvi(node,subjA,subjB) = nvi;
            ret.nmi(node,subjA,subjB) = nmi;
            ret.corr_pvalues(node,subjA,subjB) = pval;
        end
    end
end

ret.avgphi=nan(1,numnodes);
ret.avgnvi=nan(1,numnodes);
ret.avgnmi=nan(1,numnodes);

for node=1:numnodes
    % Select the upper triangular part of the phi ND array with the phi
    % coefficients of each node for every pair of subjects
    % put nan in the lower triangular and then ignores them
    phivals = triu(squeeze(ret.phi(node,:,:)),1) + tril(nan(numsubjs));
    phivals = phivals(~isnan(phivals));
    %phivals = triu(squeeze(phi(node,:,:)),1);
    % Compute the average over the pairs of subjects of the phi for each
    % node
    zphi = atanh(phivals-eps); % to avoid inf and nans around
    zphi = zphi(~isinf(zphi));
    ret.avgphi(node) = tanh(mean(zphi));
    
    nvis = triu(squeeze(ret.nvi(node,:,:)),1) + tril(nan(numsubjs));
    nmis = triu(squeeze(ret.nmi(node,:,:)),1) + tril(nan(numsubjs));
    ret.avgnvi(node) = nanmean(nvis(:));
    ret.avgnmi(node) = nanmean(nmis(:));
end
end