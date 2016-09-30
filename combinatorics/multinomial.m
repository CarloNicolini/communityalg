function y = multinomial(k)
	x=sum(k);
	y = factorial(x)/(prod(factorial(k)));