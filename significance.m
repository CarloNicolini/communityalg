function S = significance(A, ci)
%SIGNIFICANCE      Compute significance of a vertex partition on a binary network.
%
%
%   Inputs      A,  undirected binary network. If weighted, weights are
%                   ignored.
%               ci, membership vector
%
%   Outputs:    S,  value of Significance.
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2015).
%

if length(unique(A(:))) ~= 2
    warning on;
    warning('Input matrix is not binary {0,1}. Ignoring edge weights to compute Surprise.');
end

groups = membership2groups(ci); % convert membership vector to groups
ncomms = length(groups); % number of communities

% Significance value
S = 0;
% Density of the graph
d = density_und(A);

for c=1:ncomms
    nodes = groups{c}; % nodes in the community c
    nc = length(nodes); % number of nodes in community c
    if nc>1
        %mc = sum(sum(triu(A(nodes,nodes))));
        dc = density_und(A(nodes,nodes)); % Density of community c
        pc = nc*(nc-1)/2; % number of nodes pairs in community c
        S = S + pc*KL(dc,d); % increment significance
    end
end
S = 2*S; % to be compliant with the Traag implementation

function D = KL(q,p)
if (q==p)
    D=0;
    return;
end
D = 0.0;
if (q > 0.0 && p > 0.0)
    D = D + q*log(q/p);
end
if (q < 1.0 && p < 1.0)
    D = D + (1.0-q)*log((1.0-q)/(1.0-p));
end