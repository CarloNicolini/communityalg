function [all_levels_memb] = hierarchical_community_detection(A,cd_alg)
% HIERARCHICAL_COMMUNITY_DETECTION
% Apply recursively the community detection method "cd_alg" which takes as
% input an adjacency matrix.
% Returns an NNodes x NLevels array with the memberships of each node.
% The nodal memberships are indexed in a way such that is possible to
% reconstruct from which community at level l-1 the node i was at the level
% l.
% Each column has membership values that are strictly larger than the
% previous level, so to get the nodes at level 2 that where in the
% community 3 at level 1, you can do:
% N(N(:,1)==3,2)

all_levels_memb = cd_alg(A);
all_levels_memb=all_levels_memb(:);

condition=true;
while condition 
    cur_memb = hierarchical_helper(A, cd_alg, all_levels_memb(:,end)) + max(all_levels_memb(:,end));
    if length(unique(cur_memb)) == length(unique(all_levels_memb(:,end))) % No further splitting can find communities
        condition=false;
    else
        all_levels_memb = [all_levels_memb, cur_memb];
    end
end

function memb = hierarchical_helper(A, method, pre_memb)
memb = zeros(1,length(pre_memb));
for c=unique(pre_memb(:))'
    nodes = find(pre_memb==c);
    memb(nodes) = method(A(nodes,nodes)) + max(memb(:));
end
memb=memb(:);