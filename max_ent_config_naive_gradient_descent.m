function [C_con, alpha, beta, it] = max_ent_config_naive_gradient_descent(C, tolerance, r, transform_to_corr_mat)
% Credits to 
% https://journals.aps.org/pre/pdf/10.1103/PhysRevE.98.012312
% Code from Naoki Masuda
% https://github.com/naokimas/config_corr
% 
%
% Input
%   C: covariance matrix
%   tolerance: tolerance in relative error
%   r: learning rate
%   transform_to_corr_mat: whether to trasnform the input and output covariance matrix
%   to the correlation matrix.
%      0: no transform
%      1: transform (default)%
% Output
%   C_con: estimated covariance matrix
%   alpha, beta: parameters of the estimated covariance/precision matrix
%   it: iteration number
%
% estimated precision matrix K_est is given by
% K_est(i,i) = alpha(i)
% K_est(i,j) = beta(i) + beta(j)

N = size(C,1); % number of nodes

if transform_to_corr_mat == 1 % transform the input matrix to the correlation matrix
    C = diag(diag(C))^(-1/2) * C * diag(diag(C))^(-1/2); % corr matrix; overwritten on the original variable C
end
s = sum(C,2); % node strength including the self-loop
% 's' is a column vector

K = inv(C); % precision matrix from data
alpha = diag(K); % initialization
beta = zeros(N,1); % initialization
error = 1e+3; % initialization
it=0; % number of iteration

while (error > tolerance) 

    % construct the precision matrix from the current estimate
    K_est = repmat(beta,1,N) + repmat(beta,1,N)';
    for i=1:N
        K_est(i,i) = K_est(i,i) + alpha(i);
    end
    
    % covariance matrix from the current estimate
    C_con = inv(K_est);

    % gradient descent on the log likelihood
    alpha = alpha + r * (diag(C_con) - diag(C));
    beta = beta + r * (1/N) * (sum(C_con,2) - s); % gradient descent on the covariance matrix

    if (mod(it,1000) == 0) % Then measure the relative error
        error = (sum(abs((diag(C) - diag(C_con)) ./ diag(C))) + ...
            sum(abs(((s - sum(C_con,2)) ./ s)))) / N;
        fprintf('%f %f %f\n', error, mean(abs(alpha)), mean(abs(beta)));
    end
    it = it + 1;
end % end of while loop
end
