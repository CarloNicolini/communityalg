function memb = community_size2memb(c)
%COMMUNITY_SIZE2MEMB Return a membership vector from a community size
%vector
%
%   Input:  c, vector containing the number of nodes of communities
%   Output: memb, vector of nodes membership, length(memb)=sum(c)
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).

num_nodes=sum(c);
memb=zeros(1,num_nodes);
cumsizes = cumsum(c);

for i=1:length(c)
    nodeBeg = cumsizes(i)-c(i)+1;
    nodeEnd = cumsizes(i);
    memb(nodeBeg:nodeEnd)=i;
end
