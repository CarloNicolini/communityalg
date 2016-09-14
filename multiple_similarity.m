function [avgnmi,stdnmi] = multiple_similarity2(memberships)
%MULTIPLE_SIMILARITY Returns the average normalized mutual information and the standard deviation of a memberships matrix, where the nodes memberships of different repetition of a method are the rows 
% Inputs:        memberships a [nreps x nnodes] matrix of integers.
% Outputs:       avgnmi is the average nmi between every pair of distinct partitions
%                std is the standard deviation of nmi between every pair of distinct partitions
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