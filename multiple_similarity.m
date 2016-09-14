function [avgnmi,stdnmi] = multiple_similarity2(memberships)
%MULTIPLE_SIMILARITY Returns the average normalized mutual information and the standard deviation of a memberships matrix, where the nodes memberships of different repetition of a method are the rows 
% Inputs:        memberships a [nreps x n] 
%                method is a function handle to a community detection method.  It works with functions in this form [membership, quality] = method(adjacency)
%                nreps: is the total number of times that `method` is run
%
% CommunityAlg toolbox: a toolbox for community detection utilities in Matlab/Octave
% Carlo Nicolini, Istituto Italiano Di Tecnologia (2016)

nreps=size(memberships,1);
nmivals=[];

parfor i=1:nreps
    for j=i+1:nreps
        [~,nmival]=partition_distance(memberships(i,:),memberships(j,:));
        nmivals=[nmivals nmival];
    end
end

avgnmi=mean(nmivals);
stdnmi=std(nmivals);