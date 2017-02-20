function M = FinRMT(priceTS)
%% FinRMT
% FinRMT uses Random Matrix Theory (RMT) to create a filtered correlation 
% matrix from a set of financial time series price data, for example the
% daily closing prices of the stocks in the S&P
%% Syntax
% M=FinRMT(priceTS)
%
%% Description
% This function eigendecomposes a correlation matrix of time series
% and splits it into three components, Crandom, Cgroup and Cmarket,
% according to techniques from literature (See, "Systematic Identification
% of Group Identification in Stock Markets, Kim & Jeong, (2008).") and
% returns a filtered correlation matrix containging only the Cgroup
% components.
% The function is intended to be used in conjunction with a community
% detection algorithm (such as the Louvain method) to allow for community 
% detecion on time series based networks.
%
%
%% Inputs arguments:
% priceTS : an mxn matrix containing timeseries' of stock prices. Each column
% should be a time series for one financial instrument and each row should 
% correspond to the value of each instrument at a point in time. For example
% 32.00   9.43   127.25   ...
% 32.07   9.48   126.98   ...
% 32.08   9.53   126.99   ...
%  ...    ...     ....    ...
% No header columns or timestamp columns should be included
%
%% Outputs:
% M : The filtered correlation matrix. This matrix can be passed directly to 
% a community detection algorithm in place of the modularity matrix 
%
%% Example:
% ModularityMatrix = FinRMT(myPriceData)
%  ...
% Communities = myCommunityDectionAlg(ModularityMatrix)
% 
%% Issues & Comments
% Note that the output of this function can serve as the Modularity
% Matrix (Not the Adjacency matrix) for a generalized Community Detection 
% Algorithm. Specifically, one which does not rely on properties of the 
% Adjaceny Matrix to create the Modularity Matrix. The Louvain Method 
% and methods based on spectral decompositon are examples of such.
%
%%

logs = log(priceTS);    % Convert price data
diffs = diff(logs);     % to log returns
    
N = size(diffs,2);      % N is the number of time series
T = size(diffs,1);      % T is the lenght of each series
    
    
%Create the initial correlation matrix and ensure it's symmetric
%It should be symmetric but sometimes matlab introduces small roundoff
%errors that would prevent an IsSymmetric call from returning true.
%This folding of the matrix will suffice to solve that.
    
C = corrcoef(diffs);    % Create a correlation matrix and ensure
C = .5 * (C+C');        % it's symmetric


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

end



