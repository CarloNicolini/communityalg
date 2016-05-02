function A=clique(n)
%CLIQUE Returns adjacency matrix of a clique graph with n nodes
%
%   Inputs: n, the number of nodes of the graph
%   Outputs: A, the adjacency matrix of the clique graph
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).

A=ones(n);
A(1:n+1:n*n)=0;
