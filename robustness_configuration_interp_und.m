function [orig_groups, pert_groups, p_ij] = robustness_configuration_interp_und(a_ij, alpha, n_reps, method, theta)
    if nargin==0
        error('Invalid call to community_robustness. Correct usage is:\n [c1,c2]=community_robustness(a_ij, sigma, theta, n_reps, method)');
    end
    % Evaluate the robustness of community structure by perturbing the network
    % number of nodes
    n = length(a_ij);
    % number of edges
    m = number_of_edges(a_ij);
    % probability matrix
    p_ij = double(zeros(size(a_ij)));

    fprintf('Original |E|=%d\tThresholded |E|=%d\n', m, number_of_edges(a_ij));

    % Compute the original community structure, membership and groups
    orig_membership = reindex_membership(method(a_ij));
    orig_groups = sort_group_by_size(membership2groups(orig_membership));

    for i=1:n_reps
        % Create the perturbed matrix (very slow)
        pert_a_ij = randomizer_bin_und(a_ij,alpha);
        
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
        fprintf('Iter: %2.1f\tPerturbed |E|=%d, |nComp|=%d, |ncomms|=%d\n',i/n_reps*100,npert_edges,npert_conn_comps,length(unique(memb)));
    end
    % Get the average probability that two nodes are in the same community.
    p_ij = p_ij/n_reps;

    % Get the subcomponents
    p_ij_thr = p_ij;
    p_ij_thr(p_ij_thr<theta) = 0.0;
    p_ij_thr(p_ij_thr>=theta) = 1.0;

    % Compute the components from the probability matrix
    pert_membership = reindex_membership(get_components(p_ij_thr));
    pert_groups = sort_group_by_size(membership2groups(reindex_membership(pert_membership)));
