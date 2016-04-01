function [orig_membership, pert_membership] = robustness_edge_weight_und(a_ij, perturb_thresh, perturb_ampl, n_reps, method, rethreshold, prob_threshold)
if nargin==0
    error('Invalid call to robustness_edge_weight_und. Correct usage is:\n [c1,c2]=community_robustness(a_ij, sigma, theta, n_reps, method)');
end
% Evaluate the robustness of community structure by perturbing the network
% number of nodes
n = length(a_ij);
% number of edges
m = number_of_edges(a_ij);
% probability matrix
p_ij = zeros(size(a_ij));

fprintf(2,'Original |E|=%d\n', m);
% Compute the original community structure, membership and groups
orig_membership = reindex_membership(method(a_ij));
orig_groups = membership2groups(orig_membership);

parfor i=1:n_reps
    if isoctave()
        text_waitbar(i/n_reps,'Network perturbation...');
    end
    % Create the noise to perturb the network
    noise = rand(size(a_ij))*2*perturb_ampl - perturb_ampl; % generate random noise in [-perturb_ampl,perturb_ampl]
    noise = tril(noise); noise = noise+noise';  % to keep distribution uniform and noise symmetric
    % Clear diagonal noise (no noise on self-loops)
    noise(1:n+1:n*n) = 0;
    % Noise acts only on edges
    noise = noise .* (a_ij~=0);
    % Create the perturbed matrix, by summing noise
    pert_a_ij = a_ij + noise;

    % h1 = histogram(nonzeros(pert_a_ij));  set(h1,'FaceColor','r','EdgeColor','b','facealpha',0.1);
    % h2 = histogram(nonzeros(a_ij));  set(h2,'FaceColor','r','EdgeColor','r','facealpha',0.1);
    % hold off;
    % pause;
    
    % Get the minimum weight after perturbation
    minw = min(nonzeros(pert_a_ij));
    % Get the maximum weight after perturbation
    maxw = max(nonzeros(pert_a_ij));
    
    if rethreshold
        % Threshold the network to binary, given threshold theta
        pert_a_ij  = double(threshold_absolute(pert_a_ij,perturb_thresh)~=0);
    end
    
    % Count statistics on number of edges and components
    npert_edges = number_of_edges(pert_a_ij);
    npert_conn_comps = length(unique(get_components(pert_a_ij)));
    
    % Get the community membership on the thresholded perturbed network
    memb = method(pert_a_ij); memb = reindex_membership(memb);
    % Compute the comembership matrix from membership vector
    d_ij = generate_connected_components(memb);
    % Sum the current d_ij, 1 if two nodes belong the same community, 0 otherwise.
    p_ij = p_ij + d_ij;
    % Print some info
    fprintf(2,'Iter: %2.1f\tw∈[%f,%f]\tPerturbed |E|=%d, |nComp|=%d, |ncomms|=%d\n',i/n_reps*100,minw,maxw,npert_edges,npert_conn_comps,length(unique(memb)));
end
% Get the average probability that two nodes are in the same community.
p_ij = p_ij/n_reps;

perc_thresh_pij = threshold_by_giant_component(p_ij);
fprintf(2,'Used threshold on consensus matrix=%g, suggested minimum threshold to keep consensus matrix connected=%g\n',prob_threshold, perc_thresh_pij);
% Get the subcomponents
p_ij_thr = p_ij;
p_ij_thr(p_ij_thr<prob_threshold) = 0.0;
p_ij_thr(p_ij_thr>=prob_threshold) = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pert_membership = reindex_membership(get_components(p_ij_thr));
pert_groups = membership2groups(pert_membership);
