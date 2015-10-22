function [N, comps, comps_size] = number_connected_components(A)
    [comps, comps_size] = get_components(A);
    N = length(comps_size);
