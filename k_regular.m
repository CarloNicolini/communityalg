function el = k_regular(n,k)
%KREGULAR Create a k-regular graph.
% Note: No solution for k and n both odd.
%
% @inputs n, # of nodes
% @input k, 1xN vector with degree of each vertex
% @output el, Kx3 edge list 
%
% Other routines used: symmetrizeEdgeL.m

% Updated: support for ismember 'rows' in MATLAB. 
% IB: last updated, 3/24/14

el=[0,0,0];

if k>n-1; fprintf('a simple graph with n nodes and k>n-1 does not exist\n'); return; end
if mod(k,2)==1 && mod(n,2)==1; fprintf('no solution for *n* and *k* both odd\n'); return; end


half_degree=floor(k/2);  % k/2 if k even, else (k-1)/2
    
for node=1:n
    for kk=1:half_degree
      
        node_f=mod(node+kk,n);
        if node_f==0; node_f=n; end
   
        if not(ismember([node,node_f,1],el,'rows'))
          el = [el; node node_f 1]; %#ok<AGROW>
        end
          
        node_b=mod(node-kk,n);
        if node_b==0; node_b=n; end
        
        if not(ismember([node,node_b,1],el,'rows'))
          el = [el; node node_b 1]; %#ok<AGROW>
        end

        
    end
end

if mod(k,2)==1 && mod(n,2)==0
    % connect mirror nodes
    for node=1:n/2
        
        node_m=mod(node+n/2,n);
        if node_m==0; node_m=n; end
        
        if not(ismember([node,node_m,1],el,'rows'))
          el = [el; node node_m 1]; %#ok<AGROW>
        end
        
    end 
end

el = el(2:end, :);
el=symmetrizeEdgeL(el);
el=edgeL2adj(el);

function el=symmetrizeEdgeL(el)

el2=[el(:,1), el(:,2)];

for e=1:size(el,1)
    ind=ismember(el2,[el2(e,2),el2(e,1)],'rows');
    if sum(ind)==0; el=[el; el(e,2), el(e,1), el(e,3)]; end
end


function adj=edgeL2adj(el)

nodes=sort(unique([el(:,1) el(:,2)])); % get all nodes, sorted
adj=zeros(numel(nodes));         % initialize adjacency matrix

% across all edges
for i=1:size(el,1); adj(find(nodes==el(i,1)),find(nodes==el(i,2)))=el(i,3); end