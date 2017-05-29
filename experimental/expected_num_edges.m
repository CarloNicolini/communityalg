function M = expected_num_edges(n,pm,pw)
M = n.*(n-1)./2.*pm.*pw./((1-pw+pm.*pw));