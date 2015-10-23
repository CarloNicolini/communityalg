function groups = membership2groups(membership)
%MEMBERSHIP2GROUPS  Convert membership vector to grouping of vertices in a cell.
%
%   groups = membership2groups(membership)
%   Input:   membership, membership vector. Must be linearized.
%   Output:  groups, a cell containing for every single community, the set
%   of nodes belonging to it.
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2015).
%

if max(membership) ~= length(unique(membership))
    error('Non continuous membership vector. Must reindex membership.');
end

groups = {};
un = unique(membership);
for k = 1:length(un)
    l = un(k);
    groups{l} = find(membership == l);
end