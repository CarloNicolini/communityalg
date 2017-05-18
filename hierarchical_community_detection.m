function [dendr,quals] = hierarchical_community_detection(A,method)

N=length(A);
all_levels_memb = zeros(1,N);

while length(unique(all_levels_memb(end,:))) < N
	all_levels_memb = [all_levels_memb; hierarchical_helper(A,all_levels_memb(end,:)];
end

function memb = hierarchical_helper(A,pre_memb)
	memb = zeros(1,length(A));
	for c=unique(pre_memb(:))'
	    nodes = find(pre_memb==c);
	    memb(nodes) = method(A(nodes,nodes)) + max(memb);
	end
