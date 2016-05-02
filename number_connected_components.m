function [N, comps, comps_size] = number_connected_components(A)
%NUMBER_CONNECTED_COMPONENTS    Returns the number of connected components in the adjacency matrix A
%
% Inputs:   A, the adjacency matrix
% Ouputs:   N, the number of connected components
%           comps, a cell array with all the components as adjacency
%           matrices
%           comps_size, an array with the number of nodes of all components
[comps, comps_size] = get_components(A);
N = length(comps_size);
