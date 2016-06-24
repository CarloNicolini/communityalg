function D = KL(q,p)
%KL Binary Kullback-Leibler divergence
% Returns the binary Kullback-Leibler divergence with natural logarithms
%
%   Inputs:
%           q probability distribution q
%           p probability distribution p
% Note:
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).

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
