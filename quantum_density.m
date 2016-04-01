function rho = quantum_density(A)

rho = graph_laplacian(A)/sum(sum(A));
