function W = expected_total_weight(n,pm,pw)

W = n.*(n-1)./2.*pm.*pw./((1-pw).*(1-pw+pm.*pw));