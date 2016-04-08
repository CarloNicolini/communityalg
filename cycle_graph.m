function A=cycle_graph(n)

A=zeros(n);
for i=1:n-1
    A(i+1,i)=1;
    A(i,i+1)=1;
end
A(n,1)=1;
A(1,n)=1;