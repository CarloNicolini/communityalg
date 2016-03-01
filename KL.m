function D = KL(q,p)
if (q==p)
    D=0;
    return;
end

D = 0.0;
if (q > 0.0 && p > 0.0)
    D = D + q*log(q/p);
end

if (q < 1.0 && p < 1.0)
    D = D + (1.0-q)*log((1.0-q)/(1.0-p));
end
