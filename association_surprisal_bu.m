function [p,KTot,Kin] = association_surprisal_bu(A,ci)
%bincoeff(kc,m_c) è il numero di combinazioni di kc stubs in insiemi di m_c/2 stubs
%il numero totale di stubs è 2m il numero totale di intraedges è sum(mc)
%bincoeff(kc,mc)/bincoeff(2m,m)

[B,C,KTot,~,m] = comm_mat(A,ci);
ncomm = size(B,1);
% number of stubs in comm c is twice number of edges in c 
Kin=2*diag(B)';
sKTot=sum(KTot);
sKin=sum(Kin);

fprintf('N=%d\n',sKTot);
fprintf('n=%d\n',sKin);

fprintf('x=\n');
disp(Kin(:)');
fprintf('M=\n');
disp(KTot(:)');

den = logbincoeff(sKTot,sKin);
logpis=bsxfun(@(i,j)(logbincoeff(i,j)-den),KTot(:),Kin(:));
pis=bsxfun(@(i,j)(bincoeff(i,j)/bincoeff(sKTot,sKin)),KTot(:),Kin(:));
sum(pis)
p=-sum(logpis);
