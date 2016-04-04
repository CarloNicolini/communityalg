function rho = quantum_density(A)
% S. Braunstein, S. Ghosh, S. Severini, 
% The laplacian of a graph as a density matrix: 
% a basic combinatorial approach to separability of mixed states,
% Ann. of Combinatorics, 10, no 3 (2006), 291-317.
rho = graph_laplacian(A)/sum(sum(A));
