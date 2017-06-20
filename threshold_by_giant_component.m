function [At, threshold] = threshold_by_giant_component(A)

% Do a series of bisection up and low to compute the threshold where the
% graph looses components


% INPUT: Function f, endpoint values a, b, tolerance TOL, maximum iterations NMAX
% CONDITIONS: a < b, either f(a) < 0 and f(b) > 0 or f(a) > 0 and f(b) < 0
% OUTPUT: value which differs from a root of f(x)=0 by less than TOL
%
% N ← 1
% While N ≤ NMAX # limit iterations to prevent infinite loop
%   c ← (a + b)/2 # new midpoint
%   If f(c) = 0 or (b – a)/2 < TOL then # solution found
%     Output(c)
%     Stop
%   EndIf
%   N ← N + 1 # increment step counter
%   If sign(f(c)) = sign(f(a)) then a ← c else b ← c # new interval
% EndWhile
% Output("Method failed.") # max number of steps exceeded
threshold=realmin('double');
if number_connected_components(A) > 1
    fprintf(2,'Already more than connected components\n');
    threshold = min(A(:));
    return;
end

a = max(A(:));
b = min(A(:));
TOL = 1E-12;
i = 1;
imax = 10000;

condition_iterations = true;
condition_tolerance = true;
while ( condition_iterations && condition_tolerance )
    condition_iterations = i<imax;
    condition_tolerance = abs(b-a)/2 > TOL;
    %fprintf(2,' %g < T < %g %g\n',b,a,b-a);
    c = (a + b)/2;
    if ( (number_connected_components(threshold_absolute(A,c)~=0) -1 == 0) || (b-a)/2 < TOL )
        threshold = c;
    end
    if sign(number_connected_components(threshold_absolute(A,c)~=0)-1) == sign(number_connected_components(threshold_absolute(A,a)~=0)-1)
        a = c;
    else
        b = c;
    end
    i = i + 1;
end

At = threshold_absolute(A,threshold);