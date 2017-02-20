function logScm = cm_likelihood(A,Ci)
	[B,C,K,n,m,p,Bnorm,commsizes] = comm_mat(A,Ci);

	k=degrees_und(A);
	mc = diag(B);
	nc = commsizes;
	num_comms = length(commsizes);

	logScm = 0;
	for c=1:size(C,1)
		logScm = logScm + logmultinomialstirling(k(Ci==c));
	end

	logScm = logScm - logmultinomialstirling(k);

