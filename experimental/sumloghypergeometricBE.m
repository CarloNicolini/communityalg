function logS = sumloghypergeometricBE(win,Min,w,M)

logS = 0;
logS = logC(win,win+Min-1) + logC(w-win,w-win+M-Min-1) - logC(w,w+M-1);
logS = -logS;
