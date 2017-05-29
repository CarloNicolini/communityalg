function P = ewrg_prob(W)

n=size(W,1);
Wtriu=triu(W);

pairs = n*(n-1)/2;
L = sum(sum(Wtriu>0));
S = sum(sum(Wtriu));

pm = (L.^2)./((pairs - L).*(S-L));
pw = (S-L)/S;

logZ = pairs*(log((1-pw)/(1-pw+pm*pw)));

logP=0;
for i=1:n
    for j=(i+1):n
        w=W(i,j);
        logP = logP + log(pm^sign(w)*pw^w)/logZ;
    end
end

P=logP;
end


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