function [X,iter] = nearcorr(A,tol,flag,maxits,n_pos_eig,w,prnt)
%NEARCORR    Nearest correlation matrix.
%   X = NEARCORR(A,TOL,FLAG,MAXITS,N_POS_EIG,W,PRNT)
%   finds the nearest correlation matrix to the symmetric matrix A.
%   TOL is a convergence tolerance, which defaults to 16*EPS.
%   If using FLAG == 1, TOL must be a 2-vector, with first component
%   the convergence tolerance and second component a tolerance
%   for defining "sufficiently positive" eigenvalues.
%   FLAG = 0: solve using full eigendecomposition (EIG).
%   FLAG = 1: treat as "highly non-positive definite A" and solve
%             using partial eigendecomposition (EIGS).
%   MAXITS is the maximum number of iterations (default 100, but may
%   need to be increased).
%   N_POS_EIG (optional) is the known number of positive eigenvalues of A.
%   W is a vector defining a diagonal weight matrix diag(W).
%   PRNT = 1 for display of intermediate output.
 
%   By N. J. Higham, 13/6/01, updated 30/1/13, 15/11/14, 07/06/15.
%   Reference:  N. J. Higham, Computing the nearest correlation
%   matrix---A problem from finance. IMA J. Numer. Anal.,
%   22(3):329-343, 2002.
 
if ~isequal(A,A'), error('A must by symmetric.'), end
if nargin < 2 || isempty(tol), tol = length(A)*eps*[1 1]; end
if nargin < 3 || isempty(flag), flag = 0; end
if nargin < 4 || isempty(maxits), maxits = 100; end
if nargin < 6 || isempty(w), w = ones(length(A),1); end
if nargin < 7, prnt = 0; end
 
n = length(A);
if flag >= 1
   if nargin < 5 || isempty(n_pos_eig)
      [V,D] = eig(A); d = diag(D);
      n_pos_eig = sum(d >= tol(2)*d(n));
   end
   if prnt, fprintf('n = %g, n_pos_eig = %g\n', n, n_pos_eig), end
end
 
X = A; Y = A;
iter = 0;
rel_diffX = inf; rel_diffY = inf; rel_diffXY = inf;
dS = zeros(size(A));
 
w = w(:); Whalf = sqrt(w*w');
 
while max([rel_diffX rel_diffY rel_diffXY]) > tol(1)
 
   Xold = X;
   R = Y - dS;
   R_wtd = Whalf.*R;
   if flag == 0
      X = proj_spd(R_wtd);
   elseif flag == 1
      [X,np] = proj_spd_eigs(R_wtd,n_pos_eig,tol(2));
   end
   X = X ./ Whalf;
   dS = X - R;
   Yold = Y;
   Y = proj_unitdiag(X);
   rel_diffX = norm(X-Xold,'fro')/norm(X,'fro');
   rel_diffY = norm(Y-Yold,'fro')/norm(Y,'fro');
   rel_diffXY = norm(Y-X,'fro')/norm(Y,'fro');
   iter = iter + 1;
   if prnt
      fprintf('%2.0f:  %9.2e  %9.2e  %9.2e', ...
               iter, rel_diffX, rel_diffY, rel_diffXY)
      if flag >= 1, fprintf('  np = %g\n',np), else fprintf('\n'), end
   end
   if iter > maxits
       error(['Stopped after ' num2str(maxits) ' its. Try increasing MAXITS.'])
   end
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%
function A = proj_spd(A)
%PROJ_SPD
 
if ~isequal(A,A'), error('Not symmetric!'), end
[V,D] = eig(A);
A = V*diag(max(diag(D),0))*V';
A = (A+A')/2; % Ensure symmetry.
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,n_pos_eig_found] = proj_spd_eigs(A,n_pos_eig,tol)
%PROJ_SPD_EIGS
 
if ~isequal(A,A'), error('Not symmetric!'), end
k = n_pos_eig + 10; % 10 is safety factor.
if k > length(A), k = n_pos_eig; end
opts.disp = 0;
[V,D] = eigs(A,k,'LA',opts); d = diag(D);
j = (d > tol*max(d));
n_pos_eig_found = sum(j);
A = V(:,j)*D(j,j)*V(:,j)';  % Build using only the selected eigenpairs.
A = (A+A')/2; % Ensure symmetry.
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = proj_unitdiag(A)
%PROJ_SPD
 
n = length(A);
A(1:n+1:n^2) = 1;