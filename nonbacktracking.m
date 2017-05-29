function B = nonbacktracking(A)
I = eye(size(A));
D=diag(sum(A));

B = [ zeros(size(A)), D-I; -I, A];