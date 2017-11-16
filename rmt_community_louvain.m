function [all_memb, filtmat] = rmt_community_louvain(X, varargin)
% This function returns a hierarchical decomposition of the timeseries
% contained in the multidimensional matrix X. The matrix has three indices
% in the form [NSubjects, NSamples, NRegions]. If your matrix has a
% different shape, you need to rearrange the dimensions the indices using
% the function "permute" in Matlab.
%
% This anonymous function is a wrapper for the network generator process
% that just starts from the time series and the selected ROIs
% This "generator" function helps to compute the average on the prefiltered subjects
% as follows:
% 1. Filter each subject with the included RMTFilter function
% 2. Convert the cells back to an array and Fisher transform it, by taking
% the atanh.
% 3. Compute the averages on the Fisher transformed matrix
% 4. Convert it back with inverse Fisher transformation to be it a proper
% correlation matrix.
% The inputs to the generator anonymouse function must be the TimeSeries
% and the ROIs one is interested to look. At the first run, all ROIs are
% analyzed. Subject level filtering
generator_pre_filt = @(TS,ROIs)tanh(mean(atanh(reshape(cell2mat(arrayfun(@(subj)RMTFilter(squeeze(TS(subj,:,ROIs))), 1:size(TS,1),'UniformOutput',false)), [length(ROIs),length(ROIs),size(TS,1)])),3));

% another possibility is to filter the inverse Fisher transformed of the average of the Fisher transformed
% correlation matrices of the subjects. Group level filtering
generator_post_filt = @(TS,ROIs)RMTFilter(tanh(mean(atanh(reshape(cell2mat(arrayfun(@(subj)corrcoef(squeeze(TS(subj,:,ROIs))),1:size(TS,1),'UniformOutput',false)),[length(ROIs),length(ROIs),size(TS,1)])),3)),length(ROIs),size(TS,2));

%% another possibility is to concatenate all the time series
generator_concatenate = @(TS,ROIs)RMTFilter(catalongfirstdim(squeeze(TS(:,:,ROIs))));

% Reference:
% Uncovering hidden functional brain organization by random matrix theory
% Assaf Almog, Ori Roethler, Renate Buijink, Stephan Michel, Johanna H Meijer, Jos H. T. Rohling, Diego Garlaschelli
% https://arxiv.org/abs/1708.07046
%
% Community Detection for Correlation Matrices
% Mel MacMahon and Diego Garlaschelli
% Phys. Rev. X 5, 021006 – Published 14 April 2015
% Code adapted from:
% https://it.mathworks.com/matlabcentral/fileexchange/49011-random-matrix-theory--rmt--filtering-of-financial-time-series-for-community-detection

if nargin==1
    generator = generator_post_filt;
end
if nargin==2
    switch varargin{1}
        case 'prefiltering'
            generator = generator_pre_filt;
        case 'postfiltering'
            generator = generator_post_filt;
        case 'concatenate'
            generator = generator_concatenate;
        case 'custom'
            generator = varargin{2};
        otherwise
            error('Unspecified network generator process');
    end
end
% This is the wrapper around the community detection method to use
% In this case we select the best result over 100 independent runs of the
% community_louvain method. Z will be every time a different matrix
% obtained from the subset of timeseries from a specific community
cd_wrapper = @(Z)method_best(Z, @(Z)community_louvain(Z,[],[], Z),10);

% Finally, this calls the hierarchical_community_detection on the input
% matrix X with the method defined previously and on the nodes selected
% The "method" wrapper is necessary to pass two variables to the  community detection method.
% It is responsible for the generation of the right input correlation filtered matrix, as selected from the nodes
% fed as input. ...
method_hierarchical = @(X,method, nodes) hierarchical_timeseries_community_detection( X, method, nodes);

% Apply the hierarchical decomposition on all ROIs defined by 1:size(X,3)
all_memb = method_hierarchical(X,...
    @(X,nodes)cd_wrapper(generator(X,nodes)),...
    1:size(X,3));

filtmat = generator(X,1:size(X,3));
end

function all_levels_memb = hierarchical_timeseries_community_detection(TS,method,nodes)
% HIERARCHICAL_TIMESERIES_COMMUNITY_DETECTION
% Apply recursively the community detection method "cd_alg" which takes as
% input an adjacency matrix. This function has a slight variation with
% respect to its parent, because in this variation we need to reconstruct
% an input matrix not by just slicing the adjacency matrix with indices of
% the nodes in the community, but by explicitly recomputing the correlation
% matrix from the timeseries itself.
%
%% Syntax
%
% Inputs:
% TS are the [NSubj x NSamples x NROIs] timeseries to analyze
% method : an anonymous function taking as input the timeseries, an
% anonymouse generator function and the subset of nodes.
% nodes: the subsect of nodes to analyze
%
% Outputs:
% Returns an NNodes x NLevels array with the memberships of each node.
% The nodal memberships are indexed in a way such that is possible to
% reconstruct from which community at level l-1 the node i was at the level
% l.
% Each column has membership values that are strictly larger than the
% previous level, so to get the nodes at level 2 that where in the
% community 3 at level 1, you can do:
% N(N(:,1)==3,2)
%
all_levels_memb = method(TS,nodes);
all_levels_memb=all_levels_memb(:);
condition=true;
% Hierarchically go into each community, by refiltering the timeseries
% belonging to that community and reapplying the community detection method
while condition
    cur_memb = hierarchical_helper(TS, @(X,nodes)method(X,nodes), all_levels_memb(:,end)) + max(all_levels_memb(:,end));
    if length(unique(cur_memb)) == length(unique(all_levels_memb(:,end))) % No further splitting can find communities
        condition=false;
    else %append
        all_levels_memb = [all_levels_memb, cur_memb];
    end
end
end

function memb = hierarchical_helper(A, method, pre_memb)
% This is an helper function to be called by the hierarchical community
% detection method
memb = zeros(1,length(pre_memb));
for c=unique(pre_memb(:))'
    nodes = find(pre_memb==c);
    memb(nodes) = method(A,nodes) + max(memb(:));
end
memb=memb(:); % always make the arrays column vectors
end

function M = RMTFilter(X,varargin)
%% RMTFilter
% RMTFilter uses Random Matrix Theory (RMT) to create a filtered correlation
% matrix from a set of time series data.
%% Syntax
% M=RMTFilter(X)
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
% X : an mxn matrix containing timeseries' of stock prices. Each column
% should be a time series for one financial instrument and each row should
% correspond to the value of each instrument at a point in time. For example
% 32.00   9.43   127.25   ...
% 32.07   9.48   126.98   ...
%
%% Outputs:
% M : The filtered correlation matrix. This matrix can be passed directly to
% a community detection algorithm in place of the modularity matrix
%
%% Example:
% ModularityMatrix = RMTFilter(myPriceData)
% Communities = myCommunityDectionAlg(ModularityMatrix)
%
%% Issues & Comments
% Note that the output of this function can serve as the Modularity
% Matrix (Not the Adjacency matrix) for a generalized Community Detection
% Algorithm.

if nargin==3
    C = X;
    N = varargin{1};
    T = varargin{2};
else
    % Check the input has the correct form
    % need to transform it like this because if the data comes from a squeezed dataset
    % N and T switch position and this function does not work!
    if size(X,1)==1 && size(X,2)~=1
        X=X(:);
    end
    N = size(X,2);      % N is the number of time series
    T = size(X,1);      % T is the lenght of each series
    
    %Create the initial correlation matrix and ensure it's symmetric
    %It should be symmetric but sometimes matlab introduces small roundoff
    %errors that would prevent an IsSymmetric call from returning true.
    %This folding of the matrix will suffice to solve that.
    C = corrcoef(X,'rows','complete');    % Create a correlation matrix and ensure it has no Nans
    C = .5 * (C+C');        % it's symmetric
end

if any(any(isnan(C))) || any(any(isinf(C)))
    error('There are Nan or Inf in the data, probably some values have zero variance or there are Inf and Nan in the source');
end
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

function TScat = catalongfirstdim(TS)
TScat = reshape(permute(TS,[3 2 1]),fliplr([size(TS,1)*size(TS,2),size(TS,3)]))';
end
