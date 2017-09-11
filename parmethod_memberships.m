function [memberships,quals] = parmethod_memberships(A, method, nreps)
%METHOD_MEMBERSHIPS Returns the all membership and quality values of a stochastic community detection method repeated a given number of times.
% Inputs:        A is the binary or weighted adjacency matrix 
%                method is a function handle to a community detection method.  It works with functions in this form [membership, quality] = method(adjacency)
%                nreps: is the total number of times that `method` is run
%
% CommunityAlg toolbox: a toolbox for community detection utilities in Matlab/Octave
% Carlo Nicolini, Istituto Italiano Di Tecnologia (2016)

n=length(A);
memberships=zeros(nreps,n);
quals = zeros(nreps,1);

parfor i=1:nreps
    [membership,qual] = method(A);
    membership=reindex_membership(membership);
    quals(i)=qual;
    memberships(i,:) = membership(:)';
end
