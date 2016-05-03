function logH = logHyperProbability(F, M, n, p)
logH = logC(p,M) + logC(n-p,F-M) - logC(n,F);
logH = logH/log(10);
