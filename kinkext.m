function [kin,kext,kexttot,k,mut] = kinkext(W,Ci)

if (min(Ci)<1)
    warning('this function works only with communities from 1 to C');
end

n=length(W);                      % number of vertices
k=sum(W,2);                       % degrees of nodes
L=(W~=0)*diag(Ci);                % neighbor community affiliation
kin=zeros(n,1);                   % community-specific neighbors
kext=zeros(n,max(Ci));            % number neighbors of node i (row) to nodes in community c (column)

% qui possibile errore se si indicizzano le membership partendo da 0
for c=1:max(Ci)
   kic = sum(W.*(L==c),2);
   nodesinc = find(Ci==c);
   kext(:,c)=kic;
   kin(nodesinc) = kic(nodesinc);
end

kexttot = k-kin;
mut = 1-kin./k;
% to check one sees that sum(kext,2)-kexttot == kin for all nodes