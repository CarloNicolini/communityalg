function [mc, pc, m, p] = partition_params(A,memb)

C = length(unique(memb)); %ncomms
mc = zeros(C,1);
nc = zeros(C,1);
pc = zeros(C,1);
m = sum(sum(triu(A)));
p = length(A)*(length(A)-1)/2;

for i=1:C
    nodes = find(memb==i);
    subgraph = triu(A(nodes,nodes));
    mc(i) = sum(subgraph(:));
    nc(i) = length(subgraph);
    pc(i) = nc(i)*(nc(i)-1)/2;
end