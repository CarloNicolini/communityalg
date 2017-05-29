function W = ewrg(n,pm,pw,wmax)
% Generate the weighted adjacency matrix from the Enhanced Weighted Random
% Graph model.

%L = n*(n-1)/2*pm*pw/(1-pw+pm*pw);
% This is the expected average total weight
%S = n*(n-1)/2*pm*pw/((1-pw)*(1-pw+pm*pw));

% First generate the topological adjacency matrix based on the topological
% edge picking probability
A = (rand(n)<pij(pm,pw));
% Generate the matrix of edge weights, based on sampling from the edge
% weights from the qij distribution
Q = sample_weights(n,pm,pw,wmax);
% The weighted adjacency matrix 
W = A.*Q;
W=triu(W);
W(1:n+1:n^2)=0;
W=W+W';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = sample_weights(n,pm,pw,wmax)
% Check the histogram with
% hold on; histogram(randp(Q,1,1000),'Normalization','probability'); plot(1:wmax+1,Q,'ro-'); hold off;
w = randp(qij(pm,pw,1:wmax),n);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y  = qij(pm,pw,w,varargin)
switch nargin
    case 3
        % The edge picking probability
        y = ((pm.^sign(w)).*(pw.^w).*(1-pw))./(1 - pw + pm*pw);
    case 4
        wmax = varargin{1};
        y = ((pm.^sign(w)).*(pw.^w))/(sum(arrayfun(@(v)(pm.^sign(v)).*(pw.^v),0:wmax)));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = pij(pm,pw)
p = pm.*pw./(1-pw+pm.*pw);
end