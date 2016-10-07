function b = bincoeff (n, k)
%BINCOEFF
% Input: n, 
%        k,
%
% Author: KH <Kurt.Hornik@wu-wien.ac.at>
% Created: 8 October 1994
% Adapted-By: jwe
% 
% Adapted and simplified for Matlab by Carlo Nicolini
% Carlo Nicolini, Istituto Italiano di Tecnologia (2016).

if (nargin ~= 2)
    print_usage ();
end

b = zeros (size (n));
ok = (k >= 0) & (k == fix (k)) & (~isnan (n));
b(~ok) = nan;

n_int = (n == fix (n));
idx = n_int & (n < 0) & ok;
b(idx) = (-1) .^ k(idx) .* exp (gammaln (abs (n(idx)) + k(idx)) - gammaln (k(idx) + 1) - gammaln (abs (n(idx))));

idx = (n >= k) & ok;
b(idx) = exp (gammaln (n(idx) + 1) - gammaln (k(idx) + 1) - gammaln (n(idx) - k(idx) + 1));

idx = (~n_int) & (n < k) & ok;
b(idx) = (1/pi) * exp (gammaln (n(idx) + 1) - gammaln (k(idx) + 1) + gammaln (k(idx) - n(idx)) + log (sin (pi * (n(idx) - k(idx) + 1))));

% Clean up rounding errors.
b(n_int) = round (b(n_int));

idx = ~n_int;
b(idx) = real (b(idx));