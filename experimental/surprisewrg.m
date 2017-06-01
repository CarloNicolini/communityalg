function Scm = surprisewrg(A,memb)
error('still not implemented');
[B,~,Kc,~,m] = comm_mat(A,memb);
mc = diag(B);
Scm = compute_surprise( 4*m^2, sum(Kc.^2), m, sum(mc) );