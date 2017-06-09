function [exp_o, numerator, denominator] = expected_intracluster_ewrg(W,memb)

L = number_of_edges(W);
[Bw, C, K, N, Wtot, pairs, Bnorm, nc] = comm_mat(W,memb);
Bl = comm_mat(W~=0,memb);

pm = L^2/((nchoosek(N,2)-L)*(Wtot-L));
pw = (Wtot - L)/Wtot;

m = number_of_edges(W~=0);
pzeta = sum( nc.*(nc-1)/2  );
mzeta = sum(diag(Bl));
wzeta = sum(diag(Bw));

numerator = (1-pw)/(1-pw+pm*pw) ;
denominator = (N*(N-1)/2 - m + Wtot*P)*(1-P);

exp_o = numerator/denominator;