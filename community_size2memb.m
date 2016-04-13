function memb = community_size2memb(c)
num_nodes=sum(c);
memb=zeros(1,num_nodes);
cumsizes = cumsum(c);

for i=1:length(c)
    nodeBeg = cumsizes(i)-c(i)+1;
    nodeEnd = cumsizes(i);
    memb(nodeBeg:nodeEnd)=i;
end
