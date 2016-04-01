function S = clustering_entropy(p_ij, m)
	p = nonzeros(p_ij);
	if m == 0
		S = 0;
	else
		S = -1.0/m.*sum( p.*log2(p+eps)+(1.0-p).*log2(1.0-p+eps));
	end