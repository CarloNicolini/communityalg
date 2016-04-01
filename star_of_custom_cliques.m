function [A,memb] = star_of_custom_cliques(sizes)

n = sum(sizes);
A=zeros(n);
cumsizes = cumsum(sizes);
memb = zeros(1,n);
for i=1:length(sizes)
    nodeBeg = cumsizes(i)-sizes(i)+1;
    nodeEnd = cumsizes(i);
    A(nodeBeg:nodeEnd,nodeBeg:nodeEnd)=1;
    
    % add the connecting edge
    A(cumsizes,cumsizes)=1;
    memb(nodeBeg:nodeEnd)=i;
end

A(1:n+1:n^2)=0;

