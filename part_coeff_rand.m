function pcoef_rand = part_coeff_rand(W, ci)

s = sum(W,2);

k = sum(W>0,2);

Gc = (W~=0)*diag(ci);

pij = k*k' / sum(k);

dev = 0
for c=1:max(ci)
	ki_obs = sum(W.*Gc==c,2);
	ki_rand = sum(pij.*Gc==c,2);
	dev += (ki_obs - ki_rand).^2;
end

pcoef_rand = 1 - sqrt(dev)./k;
