function [exp_o, numerator, denominator] = expected_intracluster_wrg(W,memb)

[Bw, C, K, N, Wtot, pairs, Bnorm, nc] = comm_mat(W,memb);
Bl = comm_mat(W~=0,memb);

m = number_of_edges(W~=0);
pzeta = sum( nc.*(nc-1)/2  );
mzeta = sum(diag(Bl));
wzeta = sum(diag(Bw));

P = 2*Wtot/(N*(N-1) + 2*Wtot);

numerator = (pzeta - mzeta + wzeta*P)*(1-P);
denominator = (N*(N-1)/2 - m + Wtot*P)*(1-P);

exp_o = numerator/denominator;