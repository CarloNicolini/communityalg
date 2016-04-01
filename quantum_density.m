function sigma = quantum_density(A)
	%Annals of Combinatorics 10 (2006) 291-317 0218-0006/06/030291-27 DOI 10.1007/s00026-006-0289-3
	% The Laplacian of a Graph as a Density Matrix: A Basic Combinatorial Approach to Separability of Mixed States
	% Samuel L. Braunstein, Sibasish Ghosh, and Simone Severini
	sigma = graph_laplacian(A)/sum(strengths_und(A));