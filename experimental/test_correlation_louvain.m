clear all;
close all;
addpath('~/workspace/communityalg/');
addpath('~/workspace/BCT/2017_01_15_BCT');
clc;


%% TRY WITH THE SAME DATA OF GARLASCHELLI
T=50000; % number of time points
localnoise=0.4; % parameter of local noise
globalmode=0.4; % parameter of global mode
planted_memb = reindex_membership(community_size2memb([35,60,85,110,140,165,190,215]));
[C,X]=garlaschelli_benchmark_corr(planted_memb,T,localnoise,globalmode);

M = FinRMT(zscore(X));
imagesc(M)
%[rmt_best_memb, rmt_best_qual] = correlation_louvain(C,T,[],'FTWM');
[rmt_best_memb, rmt_best_qual] = community_louvain(C,[],[],C-M);
[vi,nmi]=partition_distance(rmt_best_memb(:),planted_memb(:))

%figure;
%[Xd,Yd,indices_detected] = grid_communities(rmt_best_memb); 
%subplot(1,2,1);
%imagesc(C(indices_detected,indices_detected));
%title('RMT CORRECTION');
%hold on;
%plot(Xd,Yd,'r','linewidth',1);
%hold off;

%subplot(1,2,2);
%[Xp,Yp,indices_planted] = grid_communities(planted_memb); 
%imagesc(C(indices_planted,indices_planted));
%title('Planted membership');
%hold on;
%plot(Xp,Yp,'r','linewidth',1);
%hold off;
