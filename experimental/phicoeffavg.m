function [avgphi, phi, p_values] = phicoeffavg(all_memberships)

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
phi = nan(numnodes,numsubjs,numsubjs);
p_values = nan(numnodes,numsubjs,numsubjs);
for node=1:numnodes
    for subjA=1:numsubjs
        for subjB=1:numsubjs
            [phival, pval] = corr(delta(node,:,subjA)',delta(node,:,subjB)');
            phi(node,subjA,subjB) = phival;
            p_values(node,subjA,subjB) = pval;
        end
    end
end

avgphi=nan(1,numnodes);

for node=1:numnodes
    % Select the upper triangular part of the phi ND array with the phi
    % coefficients of each node for every pair of subjects
    % put nan in the lower triangular and then ignores them
    phivals = triu(squeeze(phi(node,:,:)),1) + tril(nan(numsubjs));
    phivals = phivals(~isnan(phivals));
    %phivals = triu(squeeze(phi(node,:,:)),1);
    % Compute the average over the pairs of subjects of the phi for each
    % node
    zphi = atanh(phivals-eps); % to avoid inf and nans around
    zphi = zphi(~isinf(zphi));
    avgphi(node) = tanh(mean(zphi));
end
end