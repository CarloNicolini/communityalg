function [p,KTot,Kin] = surprisal(A,ci)

[B,C,KTot,n,m,p] = comm_mat(A,ci);
ncomm = size(B,1);
% number of stubs in comm c is twice number of edges in c 
Kin=2*diag(B)';
sKTot=sum(KTot);
sKin=sum(Kin);

% fprintf('N=%d\n',sKTot);
% fprintf('n=%d\n',sKin);

% fprintf('x=\n');
% disp(Kin(:)');
% fprintf('M=\n');
% disp(KTot(:)');

for c=1:ncomm
	p(c) = bincoeff(KTot(c),Kin(c))*bincoeff(sKTot-KTot(c),sKin-Kin(c))/bincoeff(sKTot,sKin);
end
p=-sum(log(p));