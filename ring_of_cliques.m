function A = ring_of_cliques(n,r)

A=zeros(n*r);

for i=1:r
    nodeBeg = (i-1)*n+1;
    nodeEnd = i*n;
    A(nodeBeg:nodeEnd,nodeBeg:nodeEnd)=1;
end
