function [best_membership,best_qual] = parmethod_best(A, method, nreps)
%METHOD_BEST Returns the best membership of a stochastic community detection method repeated a given number of times but in parallel (it requires a little bit more memory instead)
% Inputs:        A is the binary or weighted adjacency matrix 
%                method is a function handle to a community detection method.  It works with functions in this form [membership, quality] = method(adjacency)
%                nreps: is the total number of times that `method` is run
%
% CommunityAlg toolbox: a toolbox for community detection utilities in Matlab/Octave
% Carlo Nicolini, Istituto Italiano Di Tecnologia (2016)

n = size(A,1);
qual = zeros(nreps,1);
membership = nan(nreps,n);
parfor i=1:nreps
    [membership(i,:),qual(i)] = method(A);
end
best_qual=max(qual);
ibest_qual=find(best_qual==qual);
ibest_qual=ibest_qual(1);
best_membership = membership(ibest_qual,:);
