function [A,memb] = ring_of_custom_cliques(sizes)

n = sum(sizes);
A=zeros(n);
cumsizes = cumsum(sizes);
memb = zeros(1,n);
for i=1:length(sizes)
    nodeBeg = cumsizes(i)-sizes(i)+1;
    nodeEnd = cumsizes(i);
    A(nodeBeg:nodeEnd,nodeBeg:nodeEnd)=1;
    memb(nodeBeg:nodeEnd)=i;
end

% Add the interconnecting links
for i=1:length(sizes)-1
    A(cumsizes(i), cumsizes(i)+1)=1;
    A(cumsizes(i)+1, cumsizes(i))=1;
end
A(1,n)=1;
A(n,1)=1;
A(1:n+1:n^2)=0;
