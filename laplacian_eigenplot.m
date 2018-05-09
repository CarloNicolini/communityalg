function y = laplacian_eigenplot(A)
n=length(A);
% 1. compute the diagonal degree matrix
D = diag(sum(A));
% 2. Compute the normalized Laplacian
L = eye(n) - D^(-1/2)*A*D^(-1/2);

% 3. Compute the eigenvalues lambda
lambda = eig(L);

% 4. Define the grid over which to compute the smoothing
nx = 500;
x = linspace(0,2,nx);

% 5. Draw histogram and the distribution
hold on;
nbins = 200;
histogram(eig(L),nbins,'Normalization','countdensity');
sigma=0.005;
y = kdesmooth(eig(L),x,sigma);
plot(x, y,'r','LineWidth',2);
