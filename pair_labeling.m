function [ S ] = pair_labeling( S1, S2 )
% S1 and S2 are two different partitions for the same set of nodes
% this function fix S1, relabelling the clusters for S2 to maximize the sum
% of the same labels
if size(S1,2) ~= 1
    S1 = S1';
end
if size(S2,2) ~= 1
    S2 = S2';
end
nS1 = max(S1);
nS2 = max(S2);
cluster_intersection = zeros(nS1, nS2);
for i = 1:nS1
    for j = 1:nS2
        cluster_intersection(i,j) = -sum(S1==i & S2==j);
    end
end
Smat = munkres(cluster_intersection);
S = S2;
unlabeledCluster = nS1+1;
for i = 1:nS2
    idx = find(Smat(:,i), 1);
    if isempty(idx)
        S(S2==i) = unlabeledCluster;
        unlabeledCluster = unlabeledCluster + 1;
    else
        S(S2==i) = idx;
    end
end
end

