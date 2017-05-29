function Wij = wrg(n,p,wmax)

%p = (2*w)/(n*(n-1) + 2*w);

Wij = zeros(n);
for i=1:n
    for j=i+1:n
        % First put the edges
        w = randi([0 wmax]);
        if qij(p,w) < rand % randomly select if there is an edge or not
            Wij(i,j) = 1;
        end
    end
end

avgw = p/(1-p)

Wij(1:n+1:n^2)=0;
Wij = Wij+Wij';

function q = qij(p,w)
q = (p^w)*(1-p);