function y = logmultinomialstirling(k)
	n = sum(k);
	ki = k./n;
	y = -n*sum(xlogx(ki));