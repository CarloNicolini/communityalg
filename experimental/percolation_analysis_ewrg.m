function percolation_analysis_ewrg()
n=100;
pm=0.5;
wmax=25;
g=[];
pws=linspace(0.05,0.95,25);
nreps=25;
errg=[];
for pw=pws
    G = ewrg(n,pm,pw,wmax);
    comps = arrayfun(@(i)percolation_analysis(G), 1:nreps);
    g=[g mean(comps)];
    errg=[errg; std(comps)];
end

hold on;
plot(pws,g);
errorbar(pws,g,errg);
hold off;

function giant = percolation_analysis(T)

w=unique(nonzeros(T(:)));
%w=0.01:0.005:0.3355;
all_comps_w=[];
all_comps_size=[];
giant=0;
for t=w(:)'
    L=threshold_absolute(T,t);
    [~, comps_size]=get_components(L);
    giant = max([giant, comps_size]);
end
giant = max(all_comps_size);
