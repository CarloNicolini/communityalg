function logH = logHyperProbability(F, M, n, p, base10)
logH = logC(p,M) + logC(n-p,F-M) - logC(n,F);
if base10
	logH = logH/log(10);
end
