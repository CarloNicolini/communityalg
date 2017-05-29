function c = size_giant_component(G)
[~,compsize]=get_components(G);
c = max(compsize);