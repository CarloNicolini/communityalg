function g_ij = generate_connected_components(membership)
%GENERATE_CONNECTED_COMPONENTS      Get the connected components matrix from
%  membership vector
%
%   g_ij = generate_connected_components(membership) takes as input a set of 
%   vertex partitions of dimensions [vertex x 1] and returns a matrix g_ij 
%   containing ones if node i and node j are in the same community or zero
%   if they belong two different communities.
%
%   Inputs     membership,  affiliation (membership) vector.
%
%   Outputs:    g_ij
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2015).

if min(membership) == 0
    error('Membership vector starts from 1 to |C| included. Must reindex membership');
end
if max(membership) ~= length(unique(membership))
    error('Non continuous membership vector. Must reindex membership.');
end

groups = membership2groups(membership);
g_ij = zeros(length(membership));

for i=1:length(groups)
    g_ij(groups{i},groups{i})=1;
end
