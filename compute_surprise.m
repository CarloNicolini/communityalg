function S = compute_surprise(F, M, n, p, ignore_check_args)
%COMPUTE_SURPRISE calculates the parameter Surprise for 4 given parameters
%
%   Input:  F, number of total vertex pairs in the graph
%           M, number of total intracluster vertex pairs
%           n, number of network edges
%           p, number of intracluster network edges.
%
% You should have received this program together with
%
% If you use this program, please cite:
%       Aldecoa R, Marin I (2011)
%       Deciphering network community structure by Surprise
%       PLoS ONE 6(9): e24195

% The program receives the four parameters needed for the computation
% of Surprise. It calculates the probability of the partition based on
% a cumulative hypergeometric distribution. All calculations are done
% by using logarithms and other optimizations to avoid underflow problems.
%
% Copyright (C) 2012 Rodrigo Aldecoa and Ignacio Marin
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   Contact info: Rodrigo Aldecoa <raldecoa@ibv.csic.es>
if nargin < 5 
    argsok = check_urn_model_validity(F,M,n,p);
end

min = M;
if(n < M)
    min = n;
end

logP = logHyperProbability(F,M,n,p);
stop = false;
while ~stop && (p < min)
    p = p+1;
    nextLogP = logHyperProbability(F,M,n,p);
    [stop, logP] = sumLogProbabilities(nextLogP, logP);
end

S = -logP;

%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%
    function logH = logHyperProbability(F, M, n, p)
        logH = logC(p,M) + logC(n-p,F-M) - logC(n,F);
        logH = logH/log(10);
    end

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
            
            if(nextLogP - logP > -4)
                stop = true;
            end
        end
    end
end

% This function check that the arguments provided to compute_surprise
% function are ok in the sense that they represent a realistic urn model:
% An urn with p balls inside, pi white and p-pi black, and m balls are
% drawm and one looks for the probability that at least mi of them are
% white.
function [argsok] = check_urn_model_validity(p,pi,m,mi)
argsok = true;
if ( mi<0 )
    error('Integer overflow: negative mi');
end

if ( m<0 )
    
    error('Integer overflow: negative m');
end

if ( pi<0 )
    
    error('Integer overflow: negative pi');
end

if ( p<0 )
    error('Integer overflow: negative p');
end

if (pi>p)
    error('Error: impossible urn model p<pi');
end

if (mi>m)
    error('Error: impossible urn model m<mi');
end

if ((m-mi) > (p-pi) )
    
    error('Error: impossible urn model (p-pi)<(m-mi)');
end

if (pi>p)
    error('Error: impossible urn model pi>p');
end

if (mi > pi)
    error('Error: impossible urn modelmi>pi');
end

if (mi==0)
    argsok=false;
    return;
end

if (pi==0)
    argsok=false;
    return;
end

if ((m-mi) == (p-pi) )
    argsok=false;
    return;
end

argsok=true;
end