function Wij = wrg(n,p,wmax)

Wij = zeros(n);
for i=1:n
    for j=i+1:n
        % First put the edges
        w = randi([0 wmax]);
        if qij(p,w) > rand % randomly select if there is an edge or not
            Wij(i,j) = w;
        end
    end
end

Wij(1:n+1:n^2)=0;
Wij = Wij+Wij';

function q = qij(p,w)
q = (p^w)*(1-p);