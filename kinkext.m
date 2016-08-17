function D = kinkext(W,Ci)


k=sum(W,2); % degrees of nodes

n=length(W);                        %number of vertices
Ko=sum(W,2);                        %(out)degree
L=(W~=0)*diag(Ci);                 %neighbor community affiliation
kin=zeros(n,1);                     %community-specific neighbors

for c=1:max(Ci)
   l = sum(W.*(L==c),2);
   nodesinc = find(Ci==c);
   kin(nodesinc) = l(nodesinc);
end

D.kin = kin;
D.kext = k-kin;
D.k = k;