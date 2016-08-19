function D = kinkext(W,Ci)

if (min(Ci)<1)
    warning('this function works only with communities from 1 to C');
end

n=length(W);                        % number of vertices
D.k=sum(W,2);                       % degrees of nodes
L=(W~=0)*diag(Ci);                  % neighbor community affiliation
D.kin=zeros(n,1);                   % community-specific neighbors
D.kext=zeros(n,max(Ci));            % number neighbors of node i (row) to nodes in community c (column)

% qui possibile errore se si indicizzano le membership partendo da 0
for c=1:max(Ci)
   kic = sum(W.*(L==c),2);
   nodesinc = find(Ci==c);
   D.kext(:,c)=kic;
   D.kin(nodesinc) = kic(nodesinc);
end

D.kexttot = D.k-D.kin;
D.mut = 1-D.kin./D.k;
% to check one sees that sum(D.kext,2)-D.kexttot == D.kin for all nodes