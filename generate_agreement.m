function D = generate_agreement(W, n_reps_agreement, method)
    ci = [];
    parfor i=1:n_reps_agreement
        memb = method(W);
        ci = [ci; memb];
    end
    D = agreement(ci'+1)/n_reps_agreement;
