function A = ring_of_custom_cliques(sizes)

n = sum(sizes);
A=zeros(n);
cumsizes = cumsum(sizes)
for i=1:length(sizes)
    nodeBeg = cumsizes(i)-sizes(i)+1;
    nodeEnd = cumsizes(i);
    A(nodeBeg:nodeEnd,nodeBeg:nodeEnd)=1;
end
