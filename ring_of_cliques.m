function [A,memb] = ring_of_cliques(n,r)

sizes = repmat(n,r,1);

[A,memb] = ring_of_custom_cliques(sizes);