function res = rmtdecompose(C,T)
%RMTDECOMPOSE Decompose a correlation matrix into GlobalMode-RandomNoise-LocalModes correlation matrices
% Returns a struct with 3 variables:
% Cm is the correlation matrix of the global mode
% Cr is the correlation matrix of random noise
% Cg is the correlation matrix of remaining correlations (see them as local
% modes)
% References:
% "MacMahon, Garlaschelli", Community detection for correlation matrices,
% PhysRev X,5,021006.

N=length(C);
% Compute the predicted Marcenk-Pastur upper and lower bounds on eigen
% value distribution of random correlation
res.lambda_plus =  (1+sqrt(N/T))^2;
res.lambda_minus = (1-sqrt(N/T))^2;

% Predicted spectrum from random matrix theory
[V,D] = eig(C);
if ~issorted(diag(D))
    [V,D] = eig(A);
    [D,I] = sort(diag(D));
    V = V(:, I);
end
% Check that the eigendecomposition is fine
fprintf('Eigendecomposition abs error is %g\n', sum(sum(C*V-V*D)));
D=D.*(D>=0); % set the very small negative eigenvalues to zero
eigenvals=diag(D); % eigenvalues as 1D array sorted with the largest at the end

% Decompose the original correlation matrix into three components
% C = Cr + Cg + Cm
% Cm is the eigendecomposition of global mode Cm = lambda_max |v_max> <v_max|
% Cg is the eigendecomposition of remaining correlations Cg = sum_{res.lambda_plus<lambda_i<lambda_m} lambda_i |v_i> < v_i|
% Cr is the eigendecomposition of random-noise correlations Cr =
% sum_{lambda_i<=res.lambda_plus} lambda_i |v_i> < v_i|

% Since D is sorted with largest eigenvalues on the last diagonal element,
% the index of lambda_m is N
% Build the market mode correlation matrix
im = find(eigenvals==max(eigenvals));
res.eigenvals=eigenvals; % copy it to the output array
res.Cm = D(im,im).*V(:,im)*V(:,im)';


% Build the random-noise correlations
res.Cr = zeros(N);
ir = find(eigenvals<=res.lambda_plus)';
for k=ir
    res.Cr = res.Cr + D(k,k).*(V(:,k)*V(:,k)');
end

% Build the remaining correlations as C = Cm + Cr + Cg
res.Cg = zeros(N);
ig = find(eigenvals<max(eigenvals) & eigenvals>res.lambda_plus)';
for k=ig
    res.Cg = res.Cg + D(k,k).*(V(:,k)*V(:,k)');
end

res.Cs = res.Cg + res.Cm;

% Build the infinite time series correlation hypothesis (no correlation)
res.Cdelta = eye(N);

% Check that the partial eigendecomposition worked smoothly
fprintf('Max-min values of difference between C and its Marcenko-Pastur eigendecomposition [%g %g]\n', max(max(C-res.Cr-res.Cm-res.Cg)),min(min(C-res.Cr-res.Cm-res.Cg)));
