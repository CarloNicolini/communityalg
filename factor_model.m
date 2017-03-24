function [C,Y,P] = factor_model(ci,T,eta,mu)
% FACTOR_MODEL Computes a benchmark correlation matrix with
% controllable parameters for the local noise and market mode noise
% See MacMahon, Garlaschelli, Community detection for correlation matrices,
% PhysRev X,5,021006. Section IV, D "Benchmarking our methods"
% Parameters:
% membership: the membership vector to generate the correlation matrix
% T: the number of time samples on which the correlation are estimated
% eta: the local noise parameter
% mu: the market mode parameter
N=length(ci);

if T<N
    warning('This benchmark is not valid as curse of dimensionality is not respected, See MacMahon Garlaschelli, PhysRevX 2015');
end

Y=nan(T,N);

% It's important that the number of time samples is greater than the number
% of nodes, otherwise this method produces communities with external
% positive correlation as a by-product

for c=unique(ci(:))'
    Y(:,ci==c) = repmat(randn(T,1),[1,sum(ci==c)]);
end

% add local noise beta
Y = Y + eta*randn(T,N);
% add global signal alpha
Y = Y + mu*repmat(randn(T,1),[1,N]);

% to use with corrcoeff use the transpose as corrcoeff accepts a matrix
% where rows are variables and observations are columns
% C is the correlation matrix, P are the P-values of correlation
% coefficients
%Y = zscore(Y);
[C,P] = corrcoef(Y);
%C = nearcorr(C);


% function [C,y,P] = garlaschelli_benchmark_corr(membership,T,eta,mu)
% % GARLASCHELLI_BENCHMARK_CORR Computes a benchmark correlation matrix with
% % controllable parameters for the local noise and market mode noise
% % See MacMahon, Garlaschelli, Community detection for correlation matrices,
% % PhysRev X,5,021006. Section IV, D "Benchmarking our methods"
% % Parameters:
% % membership: the membership vector to generate the correlation matrix
% % T: the number of time samples on which the correlation are estimated
% % eta: the local noise parameter
% % mu: the market mode parameter

% c=length(unique(membership));
% N=length(membership);

% if T<N
%     error('This benchmark is not valid as curse of dimensionality is not respected, See MacMahon Garlaschelli, PhysRevX 2015');
% end
% comm_sizes = cellfun(@length,membership2groups(membership));

% gamma=randn(T,c); % community-time-series
% y=zeros(T,N); % final result

% % It's important that the number of time samples is greater than the number
% % of nodes, otherwise this method produces communities with external
% % positive correlation as a by-product

% for A=unique(membership)
%     y(:,membership==A) = repmat(gamma(:,A),1,comm_sizes(A));
% end

% % add local noise beta
% y = y + eta*randn(T,N);
% % add global signal alpha
% y = y + mu*randn(T,N);

% % to use with corrcoeff use the transpose as corrcoeff accepts a matrix
% % where rows are variables and observations are columns
% % C is the correlation matrix, P are the P-values of correlation
% % coefficients
% y = zscore(y);
% [C,P] = corrcoef(y);
% C = nearestSPD(C);