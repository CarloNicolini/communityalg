clear all;
[L,ml]=ring_of_custom_cliques([4 4]);
qmax=association_surprisal_bu(L,ml)

n=size(L,1);
q0=-Inf;
m0=1:n;
i=0;
m=m0;
while i<1000
	q = association_surprisal_bu(L,m);
	if q>q0
		q0=q;
		m0=m;
		q
	end
	m = [randi(n/2,1,n/2) randi(n/2,1,n/2) ];
	i=i+1;
end