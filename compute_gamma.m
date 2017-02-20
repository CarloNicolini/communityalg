function g = compute_gamma(W)
% implements optimal gamma computation as in Newman 2016
% Community detection in networks:
% Modularity optimization and maximum likelihood are equivalent
% arXiv:1606.02319v

% First guess of gamma
g = 1;
for i=1:5
    Ci = community_louvain(W,g);
    %length(unique(Ci))
    [~,K,~,~,~,~,Bnorm,~] = comm_mat(W,Ci);
    omega_in = mean(diag(Bnorm));
    omega_out = mean(nonzeros(Bnorm-eye(length(Bnorm)).*Bnorm));
    g = (omega_in - omega_out)/(log(omega_in)-log(omega_out))
end
