function [effs,nodes_eff,eff] = effcommplot(graph, membership)
% effcomplot
% Makes a plot with the local efficiency of vertices grouped by community
% (different colors)
% Carlo Nicolini, 2016
% 22 March 2016, first version
% 23 March 2016, switching 12 colors with colorbrewer qualitative maps
n = size(graph,1);
[~,membsortind] = sort(membership);

nc = length(unique(membership));
colors=brewermap(12,'Set3');

% Compute community weighted efficiency for every community
effs = [];
for c=unique(membership)
    nodes = find(membership==c);
    commsubgraph = graph(nodes,nodes);
    effs(c)=efficiency_wei(commsubgraph);
end

% Compute local efficiency for all nodes
nodes_eff = efficiency_wei(graph,true);

% Compute global efficiency for the entire graph
eff = efficiency_wei(graph);

% Now sort the nodes local efficiency by the index of their membership
nodes_eff = nodes_eff(membsortind);
figure;
hold on;
for c=unique(membership)
    curcolor = colors(mod(c,11)+1,:);
    nodes = find(membership(membsortind)==c);
    bar(nodes,nodes_eff(nodes),'FaceColor',curcolor,'LineStyle','none');
    set(gca,'XTick',1:size(graph,1));
    set(gca,'XTickLabel',membsortind);
    x = linspace(min(nodes)-0.5,max(nodes)+0.5,10);
    y = ones(10,1)*effs(c);
    plot(x,y,'Color',colors(mod(c,11)+1,:)*0.8,'LineWidth',4);
end
plot(1:n,ones(n,1)*eff,'Color',[1,0,0],'LineWidth',2);

set(gca,'XTickLabelRotation',90);
set(gca,'YLim',[0,1.05]);
title('Community efficiency');
xlabel('Vertex');
ylabel('Efficiency');
hold off;