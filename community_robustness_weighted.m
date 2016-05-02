function [orig_groups, pert_groups, p_ij] = community_robustness_weighted(a_ij, perturb_ampl, n_reps, method, prob_threshold)
%COMMUNITY_ROBUSTNESS_WEIGHTED Implementation of the weighted consensus
%method
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).

if nargin==0
        error('Invalid call to community_robustness. Correct usage is:\n [c1,c2]=community_robustness(a_ij, sigma, theta, n_reps, method)');
    end
    
    % If using louvain method, which is quite noisy in its results,
    % Lancichinetti et al (Consensus clustering in complex networks)
    % suggest to set prob_threshold ~0.3
    % If using Newman spectral method on modularity, higher prob_threshold
    % has to be used, ~0.7
    
    % Evaluate the robustness of community structure by perturbing the network
    % number of nodes
    n = length(a_ij);
    % number of edges
    m = number_of_edges(a_ij);
    
    % Agreement matrix as in Lancichinetti
    p_ij = double(zeros(size(a_ij)));

    fprintf(2,'Original |E|=%d\tThresholded |E|=%d\n', m, number_of_edges(a_ij));
    % Compute the original community structure, membership and groups
    orig_membership = reindex_membership(method(a_ij));
    orig_groups = membership2groups(orig_membership);

    parfor i=1:n_reps
        if isoctave()
            text_waitbar(i/n_reps,'Network perturbation...');
        end
        % Create the noise to perturb the network
        noise = rand(size(a_ij))*2*perturb_ampl - perturb_ampl; % generate [-sigma,sigma]
        % make the noise symmetric, important to keep the network undirected
        noise = (noise + noise')/2;
        % Clear diagonal noise (no noise on self-loops)
        noise(1:n+1:n*n) = 0; 
        % Create the perturbed matrix
        pert_a_ij = a_ij + noise.*(a_ij.*(a_ij~=0)); % only add noise on edges

        % Get the minimum weight after perturbation
        minw = min(pert_a_ij(pert_a_ij~=0));
        % Get the maximum weight after perturbation
        maxw = max(pert_a_ij(:));
                
        % Count statistics on number of edges and components
        npert_edges = number_of_edges(pert_a_ij);
        npert_conn_comps = length(unique(get_components(pert_a_ij)));
        
        % Get the community membership on the thresholded perturbed network
        memb = reindex_membership(method(pert_a_ij));
        % Compute the comembership matrix from membership vector
        d_ij = generate_connected_components(memb);
        % Sum the current d_ij, 1 if two nodes belong the same community, 0 otherwise.
        p_ij = p_ij + d_ij;
        % Print some info
        fprintf(2,'Iter: %2.1f\tw∈[%f,%f]\tPerturbed |E|=%d, |nComp|=%d, |ncomms|=%d\n',i/n_reps*100,minw,maxw,npert_edges,npert_conn_comps,length(unique(memb)));
    end
    % Get the average probability that two nodes are in the same community.
    p_ij = p_ij/n_reps;

    % Run the consensuns clustering using the same method on the agreement
    % matrix p_ij
    % Get the subcomponents
    p_ij_thr = p_ij;
    p_ij_thr(p_ij_thr<prob_threshold) = 0.0;
    p_ij_thr(p_ij_thr>=prob_threshold) = 1.0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pert_membership = reindex_membership(get_components(p_ij_thr));
    pert_groups = membership2groups(pert_membership);
