function res = rmtdecompose2(C,T)
%RMTDECOMPOSE Decompose a correlation matrix into GlobalMode-RandomNoise-LocalModes correlation matrices
% Returns a struct with 3 variables:
% Cm is the correlation matrix of the global mode
% Cr is the correlation matrix of random noise
% Cg is the correlation matrix of remaining correlations (see them as local
% modes)
% References:
% "MacMahon, Garlaschelli", Community detection for correlation matrices,
% PhysRev X,5,021006.
% Check modification of lambdaplus,minus from the original code of MacMahon

N=length(C);
% Decompose the correlation matrix into its eigenvalues and eigenvectors,
% store the indices of which columns the sorted eigenvalues come from
% and arrange the columns in this order
[V,D] = eig(C);
[eigvals, ind]=sort(diag(D),'ascend'); 
V = V(:,ind);
D=diag(sort(diag(D),'ascend')); 


% Find the index of the predicted lambda_max, ensuring to check boundary
% conditions
Q=T/N;
sigma = 1 - max(eigvals)/N;
RMTmaxEig = sigma*(1 + (1.0/Q) + 2*sqrt(1/Q));
RMTmaxIndex = find(eigvals > RMTmaxEig);
if isempty(RMTmaxIndex)
    RMTmaxIndex = N;
else
    RMTmaxIndex = RMTmaxIndex(1);
end

% Find the index of the predicted lambda_min, ensuring the check boundary
% conditions
RMTminEig = sigma*(1 + (1.0/Q) - 2*sqrt(1/Q));
RMTminIndex = find(eigvals < RMTminEig);
if isempty(RMTminIndex)
    RMTminIndex = 1;
else
    RMTminIndex = RMTminIndex(end);
end

% Determine the average Eigenvalue to rebalance the matrix after removing
% Any of the noise and/or market mode components
avgEigenValue = mean(eigvals(1:RMTmaxIndex));

% Build a new diagonal matrix consisting of the group eigenvalues
Dg = zeros(N,N);

% Replace the random component with average values.
Dg(1 : (N+1) : (RMTmaxIndex-1)*(N+1)) = avgEigenValue;

% Add the group component. The N+1 here is just used to increment to the 
% next diagonal element in the matrix
Dg(1+(N+1)*(RMTmaxIndex-1) : (N+1) : end-(N+1)) = D(1+(N+1)*(RMTmaxIndex-1) : (N+1) : end-(N+1));

% Build the component correlation matrix from the new diagonal eigenvalue
% matrix and eigenvector matrix. The eigenvectors corresponding to zero
% valued eigenvalue entries in Dg will not contribute to M

M = V * Dg * V.';

% Replace the diagonals with 1s
M = M - diag(diag(M)) + eye(N);

res.Cg = M;

