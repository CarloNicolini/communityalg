function Gpq = communicability_wu(Wij)
%COMMUNICABILITY
% Communicability centrality, also called subgraph centrality, of a node `n`
% is the sum of closed walks of all lengths starting and ending at node `n`.
%    References
% For weighted networks we refer to the implementation of Higham
%
% https://pure.strath.ac.uk/portal/files/381321/strathprints013675.pdf
%
% ----------
% [1] Ernesto Estrada, Juan A. Rodriguez-Velazquez,
%       "Subgraph centrality in complex networks",
%       Physical Review E 71, 056103 (2005).
%       http://arxiv.org/abs/cond-mat/0504730

%    .. [2] Ernesto Estrada, Naomichi Hatano,
%       "Communicability in complex networks",
%       Phys. Rev. E 77, 036111 (2008).
%       http://arxiv.org/abs/0707.0756
n = length(Wij);
D=eye(n);

% Renormalization step of weights as in
% https://pure.strath.ac.uk/portal/files/381321/strathprints013675.pdf
D(1:n+1:n*n) = 1.0./sqrt(D(1:n+1:n*n)); 
Wnormij = D*Wij*D;

% Compute communicability as matrix exponential.
% From MatlabR2015b expm is very good implemented, no need of manual eigendecomposition. For Octave the performance is not known.
Gpq = expm(Wnormij);
