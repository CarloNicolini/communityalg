function [consensus_memb_pert, p_ij] = consensus_robustness(a_ij, pert_ampl, n_reps_agreement, method, consensus_reps, consensus_thresh)
    % Detect the consensus partition on a network with its edge weights
    % perturbed by a random amount pert_ampl
    n = length(a_ij);
    % Now perturb the network edge weights
    % Create the agreement matrix for the perturbed network
    p_ij = zeros(size(a_ij));
    parfor i=1:n_reps_agreement
        % Create the noise to perturb the network
        noise = rand(size(a_ij))*2*pert_ampl - pert_ampl; % generate [-sigma,sigma]
        % make the noise symmetric, important to keep the network undirected
        noise = (noise + noise')/2;
        % Clear diagonal noise (no noise on self-loops)
        noise(1:n+1:n*n) = 0; 
        % Create the perturbed matrix
        pert_a_ij = a_ij + noise.*(a_ij.*(a_ij~=0)); % only add noise on edges
        % Get the community membership on the thresholded perturbed network
        memb = reindex_membership(method(pert_a_ij));
        % Compute the comembership matrix from membership vector
        d_ij = generate_connected_components(memb);
        % Sum the current d_ij, 1 if two nodes belong the same community, 0 otherwise.
        p_ij = p_ij + d_ij;
    end
    % Get the perturbed agreement matrix
    p_ij = p_ij/n_reps_agreement;

    % Now run the consensus clustering on the agreement matrix
    consensus_memb_pert = consensus_clustering(p_ij, method, consensus_reps, consensus_thresh);

    % Finally return the perturbed consensus membership 
    consensus_memb_pert = reindex_membership(consensus_memb_pert);