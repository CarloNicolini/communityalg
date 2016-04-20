function [best_membership,best_qual] = method_best(A, method, nreps)

best_qual = -inf;
best_membership = [];
for i=1:nreps
    [membership,qual] = method(A);
    if qual > best_qual
        best_membership = membership;
        best_qual = qual;
    end
end
