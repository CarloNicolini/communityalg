function memb=powerlaw_membership(nmin,nmax,tau,num_nodes)


c = powerlaw_comm_sizes(nmin,nmax,tau,num_nodes);

memb=zeros(1,num_nodes);
cumsizes = cumsum(c);

for i=1:length(c)
    nodeBeg = cumsizes(i)-c(i)+1;
    nodeEnd = cumsizes(i);
    memb(nodeBeg:nodeEnd)=i;
end
