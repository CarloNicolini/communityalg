function logH = logHyperProbability(F, M, n, p)
logH = logC(p,M) + logC(n-p,F-M) - logC(n,F);
logH = logH/log(10);

function sum = sumFactorial(n)
sum = 0;
if(n > 1)
    for i = 2:n
        sum = sum + log(i);
    end
end
end

function sum = sumRange(min, max)
sum = 0;
for i = min:max
    sum = sum + log(i);
end
end

function x = logC(k, n)
if(k == n || ~k)
    x = 0;
else
    t = n - k;
    if(t < k)
        t = k;
    end
    x = sumRange(t+1, n) - sumFactorial(n - t);
end
end

function [stop, logP] = sumLogProbabilities(nextLogP, logP)
if(nextLogP == 0)
    stop = true;
else
    stop = false;
    if(nextLogP > logP)
        common = nextLogP;
        diffExponent = logP - common;
    else
        common = logP;
        diffExponent = nextLogP - common;
    end
    logP = common + ((log(1 + 10.^diffExponent)) / log(10));
    
    if(nextLogP - logP < -4)
        stop = true;
    end
end
end
end