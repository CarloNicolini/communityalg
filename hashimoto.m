function H = hashimoto(A)
%% Compute the Hashimoto nonbacktracking matrix from the adjacency matrix A

I=full(incidence(digraph(A)));

H = double(I>0)'*double(I<0);