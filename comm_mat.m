function [B,C,K,n,m,p] = comm_mat(W,ci)
%COMM_MAT 		Returns the block matrix of the community structure of a graph W with given membership vector ci
%
%   Inputs      W,  undirected weighted/binary network. 
%               ci, membership vector
%
%   Outputs:    B,  block matrix of community structure. Element B(i,j) contains the number of links 
%				from community to community j, half of them when i==j (no self-loops are considered).
%               C, community matrix, an ncommsX lenght(ci) matrix of ones and zeros. Each row holds a 1's where 
%               ci is equal to one of the unique value.
%               K, degree vector, K(i) is the number of stubs in community i
%               n, number of nodes in W
%               m, total weight in W
%               p, number of pairs in W, as bincoeff(n,2)
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).

if any(diag(W))
	warning('Input matrix has self-loops');
end

ci=ci(:)'; % force it to be a row vector

% The first matrix C is an ncommsX lenght(ci) matrix of ones and zeros. 
% Each row holds a 1's where ci is equal to one of the unique value.
C=bsxfun(@eq,ci,unique(ci)');
% Obtain the block matrix
B = C*W*C';
% Obtain the degrees matrix (useful for modularity and other configuration model-based methods)
K = sum(B,2);
% Divide the diagonal because self-loops are counted twice
B(logical(eye(size(B)))) = B(logical(eye(size(B))))./2;

% Other parameters
n = length(W); % number of nodes
m = sum(sum(triu(B))); % total sum of edge weight
p = n*(n-1)/2; % number of pairs