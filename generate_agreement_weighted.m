function D = generate_agreement_weighted(W, n_reps_agreement, method, percentQuant)
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

if nargin<4
    percentQuant=0;
end

all_memberships = [];
qual_weights=zeros(n_reps_agreement,1);
for i=1:n_reps_agreement
    [memb,qual] = method(W);
    % make memb ALWAYS a ROW vector
    memb = memb(:)';
    all_memberships = [all_memberships; memb];
    qual_weights(i)=qual;
end

quant_weights = quantile(qual_weights,percentQuant);
ind_quant_weights = find(qual_weights > quant_weights);
qual_weights = qual_weights(ind_quant_weights);
disp(['Agreement selected nb of repetition ' num2str(length((qual_weights)))]);
D = agreement_weighted(all_memberships(ind_quant_weights,:)'+1,qual_weights);
