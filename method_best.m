function [best_membership,best_qual] = method_best(A, method, nreps)
%METHOD_BEST Returns the best membership of a stochastic community detection method repeated a given number of times.
% Inputs:        A is the binary or weighted adjacency matrix 
%                method is a function handle to a community detection method.  It works with functions in this form [membership, quality] = method(adjacency)
%                nreps: is the total number of times that `method` is run
%
% CommunityAlg toolbox: a toolbox for community detection utilities in Matlab/Octave
% Carlo Nicolini, Istituto Italiano Di Tecnologia (2016)

best_qual = -inf;
best_membership = [];
for i=1:nreps
    [membership,qual] = method(A);
    if qual > best_qual
        best_membership = membership;
        best_qual = qual;
    end
end
