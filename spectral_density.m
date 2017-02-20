clear all;
close all;

nsamp=50;
nnodes=500;
v=zeros(nnodes,1);

parfor k=1:nsamp
    v=v+eig(load(['~/Desktop/ergraph/er' num2str(k) '.txt']));
end

v=v/(nsamp*sqrt(nnodes));

v=eig(k_regular(nnodes,3))/(sqrt(nnodes));
ksf=zeros(nsamp,100);
xi=zeros(100,1);
[f,xi]=ksdensity(v);
plot(xi,f);

%compute spectral entropy but one needs to do that at different m/N and
%make a graph
-sum(f.*log(f))
