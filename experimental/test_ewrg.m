clc;
close all;
clear all;

n=250;
reps=50;
wmax=10;
pm=0.4;
pw=0.2;
% Expected edges
L = (n*(n-1)/2) * pm*pw/(1-pw+pm*pw);
S = (n*(n-1)/2) * pm*pw/((1-pw)*(1-pw+pm*pw));

% Sample an element from this ensemble
W = ewrg(n,pm,pw,wmax);

% Variance on the expected total weight
varS = nchoosek(n,2).*pm*pw*(1+pw^2*(pm-1))./((1-pw)^2 * (1-pw+pm*pw)^2);

all_edges = arrayfun(@(i)number_of_edges(ewrg(n,pm,pw,wmax)),1:reps);
all_weight = arrayfun(@(i)total_weight(ewrg(n,pm,pw,wmax)),1:reps);
all_degrees = cellfun(@mean,arrayfun(@(i)degrees_und(ewrg(n,pm,pw,wmax)>0),1:reps,'UniformOutput',false));
all_strengths = cellfun(@mean,arrayfun(@(i)strengths_und(ewrg(n,pm,pw,wmax)>0),1:reps,'UniformOutput',false));

% Check that the observed and sampled distributions are the same
[L S]./[mean(all_edges) mean(all_weight)]
% Check that the variance of the expected total weight is correctly
% predicted
var(all_weight)./varS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check the average degree
avgK = (n-1)*pm*pw/(1-pw+pm*pw);
mean(all_degrees)./avgK

%% Check the average strength
avgS = (n-1)*pm*pw/((1-pw+pm*pw)*(1-pw));
mean(all_strengths)./avgS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% PERCOLATION ANALYSIS %%%%%%%%%
npm=10;
npw=10;
nreps=100;
[pm,pw]=meshgrid(linspace(0.05,0.2,npm),linspace(0.05,0.2,npw));

Comps = arrayfun(@(x,y)number_connected_components(ewrg(n,x,y,wmax)),pm,pw);
Giant = arrayfun(@(x,y)...
    mean(arrayfun(@(i)size_giant_component(ewrg(n,x,y,wmax)),1:nreps))...
    ,pm,pw);
surfc(pm,pw,Comps);
surfc(pm,pw,Giant);