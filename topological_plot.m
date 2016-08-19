function [retval] = topological_plot (W)
	n=length(W);
	D=sum(W).*eye(n);
	D(1:n+1:n*n) = 1.0./sqrt(D(1:n+1:n*n)); %make the renormalization step
	Wn = D*W*D;

