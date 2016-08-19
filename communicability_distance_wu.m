function Eta = communicability_distance_wu(A)
% COMMUNICABILITY_DISTANCE_BU
% Compute the communicability distance of undirected simple weighted graphs (no
% loops)
%
% Estrada, E. 
% The communicability distance in graphs.
% Linear Algebra Appl. 436, 4317–4328 (2012).
%
% Carlo Nicolini, Istituto Italiano di Tecnologia (2016).
%
n=length(A);
G = communicability_wu(A);
Eta=zeros(n);
for p=1:n
    for q=1:n
        Eta(p,q) = sqrt(G(p,p)+G(q,q)-2*G(p,q));
    end
end
