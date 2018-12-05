function [ SOptimal ] = multislice_pair_labeling(S )
% input:
%    S: pxn matrix, n is the numebr of partitions, p is the number of nodes
% output:
%    SOptimal: optimal labels re-assignament to optimize the persistence
SOptimal = zeros(size(S));
SOptimal(:,1) = S(:,1);
for i = 2:size(S,2)
    SOptimal(:,i) = pair_labeling(SOptimal(:,i-1), S(:,i));
end

