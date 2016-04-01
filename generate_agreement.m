function D = generate_agreement(W, n_reps_agreement, method)
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

    ci = [];
    parfor i=1:n_reps_agreement
        memb = method(W);
        ci = [ci; memb];
    end
    D = agreement(ci'+1)/n_reps_agreement;
