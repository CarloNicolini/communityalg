function [p,numLog] = association_surprisal_bu(A,ci)
%bincoeff(kc,m_c) è il numero di combinazioni di kc stubs in insiemi di m_c/2 stubs
%il numero totale di stubs è 2m il numero totale di intraedges è sum(mc)
%bincoeff(kc,mc)/bincoeff(2m,m)

[B,C,K,~,m] = comm_mat(A,ci);
ncomm = size(B,1);
KTot=2*sum(B);
Kin=2*diag(B)';
KTot
Kin
for i=1:ncomm
    p(i)=logbincoeff(KTot(i),Kin)-logbincoeff(sum(KTot),sum(Kin));
end
p
p=-sum(p);
%p=bsxfun(@(i,j)(logbincoeff(i,j)),K,2*M');
