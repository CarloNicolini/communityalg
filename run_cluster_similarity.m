function [clujacc, cluhyp] = run_cluster_similarity(W, pert_thresh, pert_ampl, prob_thresh, n_reps, method)

edges = (W(W~=0));

% First compute the consensus membership on unperturbed network
[consensus_unperturbed] = consensus_robustness(W, 0, n_reps_agreement, method, consensus_reps, consensus_thresh)
[clujacc,cluhyp] = cluster_similarity(c1,c2);

%
scrsz = get(groot,'ScreenSize');
fig = figure('Position',[1 1 scrsz(3) scrsz(4)]);

colormap(parula);

%% Start plotting things
subplot(1,3,1);
hist(edges,length(unique(edges)));

hold on;
line([pert_thresh pert_thresh],ylim,'Color',[0 1 0],'LineWidth',2);
line([pert_thresh-pert_ampl pert_thresh-pert_ampl], ylim, 'Color', [1 0.2 0], 'LineWidth',1);
line([pert_thresh+pert_ampl pert_thresh+pert_ampl], ylim, 'Color', [1 0.2 0], 'LineWidth',1);
hold off;
title('Weights perturbation');
subplot(1,3,2);

imagesc(clujacc);

col = colorbar;
col.Label.String = 'Jaccard';

xlabel('Orig. comm.');
ylabel('Pert. comm.');
title('Jaccard Similarity');

%%
subplot(1,3,3);

maxlogp = max(cluhyp(:));
imagesc(cluhyp);

col = colorbar;
col.Label.String = strcat('-log_{10}(p)');
xlabel('Orig. comm.');
ylabel('Pert. comm.');
title('p-score similarity (normalized)');
