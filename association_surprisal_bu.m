function p = association_surprisal_bu(A,ci)

[B,C,K,~,m] = comm_mat(A,ci);
ncomm = size(B,1);
K=sum(B);
%[i,j] = meshgrid(1:ncomm);
p=zeros(ncomm);

for i=1:ncomm
	for j=1:ncomm
		p(i,j) = log10(bincoeff(K(i)+K(j),B(i,j)));
	end
end