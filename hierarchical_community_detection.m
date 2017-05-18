function [all_levels_memb] = hierarchical_community_detection(A,method)

N=length(A);
all_levels_memb = zeros(1,N);

level = 0;
condition=true;
while condition
	all_levels_memb = [all_levels_memb; hierarchical_helper(A,method,all_levels_memb(end,:)) + max(all_levels_memb(end,:))];
	level = level+1;
    if length(unique(all_levels_memb(end,:))) == length(unique(all_levels_memb(end-1,:))) % No further splitting can find communities
        condition=false;
    end
end

function memb = hierarchical_helper(A,method,pre_memb)
	memb = zeros(1,length(A));
	for c=unique(pre_memb(:))'
	    nodes = find(pre_memb==c);
	    memb(nodes) = method(A(nodes,nodes)) + max(memb(:));
	end
