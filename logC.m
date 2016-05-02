function x = logC(k, n)
%LOGC Returns the natural logarithm of binomial coefficient n over k (n k), 
% the number of k-elements subsets of an n-element set
%
% Input:    k, the size of the subset
%           n, the size of the set
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).

if (k == n || ~k)
    x = 0;
else
    t = n - k;
    if(t < k)
        t = k;
    end
    x = sumRange(t+1, n) - sumFactorial(n - t);
end
