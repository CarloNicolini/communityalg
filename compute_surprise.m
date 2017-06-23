function S = compute_surprise(F, M, n, p, base10)
%COMPUTE_SURPRISE calculates the parameter Surprise for 4 given parameters
%
%   Input:  F, number of total vertex pairs in the graph
%           M, number of total intracluster vertex pairs
%           n, number of network edges
%           p, number of intracluster network edges.
%			base10, if true use base 10 logarithms
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
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).


argsok = check_urn_model_validity(F,M,n,p);
if nargin<5
	base10=true;
end

min = M;
if(n < M)
    min = n;
end

logP = logHyperProbability(F,M,n,p,base10);
stop = false;
while ~stop && (p < min)
    p = p+1;
    nextLogP = logHyperProbability(F,M,n,p,base10);
    [stop, logP] = sumLogProbabilities(nextLogP, logP);
end

S = -logP;