function Wij = wrg(n,p,wmax)

Wij = zeros(n);
Wij = sample_weights(n,p,wmax);

% for i=1:n
%     for j=i+1:n
%         % First put the edges
%         w = randi([0 wmax]);
%         if qij(p,w) > rand % randomly select if there is an edge or not
%             Wij(i,j) = w;
%         end
%     end
% end

Wij(1:n+1:n^2)=0;
Wij = Wij+Wij';
end

function q = qij(p,w)
q = (p.^w).*(1-p);
end

function w = sample_weights(n,p,wmax)
% Check the histogram with
%hold on; histogram(randp(Q,1,1000),'Normalization','probability'); plot(1:wmax+1,Q,'ro-'); hold off;
w = randp(qij(p,0:wmax),n);
end