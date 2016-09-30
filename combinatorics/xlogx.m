function y=xlogx(x)
	if any(x<=eps)
		y=zeros(size(x));
	end
	y = x.*log(x);
