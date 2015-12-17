function D = generate_agreement_weighted(W, n_reps_agreement, method)
%GENERATE_AGREEMENT

%GENERATE_AGREEMENT      Compute the agreement matrix for a given method
%
%   Inputs      W,  weighted or binary undirected network.
%               n_reps_agreement, number of samples
%               method, handle to a community detection method that returns
%               a row vector of community membership
%
%   Outputs:    P, the agreement matrix, average number of times that
%   vertex i and vertex j are clustered in the same community over all
%   n_reps_agreeement
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2015).

    all_memberships = [];
    qual_weights=zeros(n_reps_agreement,1);
    parfor i=1:n_reps_agreement
        [memb,qual] = method(W);
        % make memb ALWAYS a ROW vector
        memb = memb(:)';
        all_memberships = [all_memberships; memb];
        qual_weights(i)=qual;
    end
    D = agreement_weighted(all_memberships'+1,qual_weights);
    % no need of division to n_reps_agreement
