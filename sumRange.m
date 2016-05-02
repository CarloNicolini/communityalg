function csum = sumRange(xmin, xmax)
%SUMRANGE Utility function for logC
%
%   Input:  xmin, lower bound
%           xmax, higher bound
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).

csum = 0;
for i = xmin:xmax
    csum = csum + log(i);
end
