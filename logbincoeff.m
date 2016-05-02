function c = logbincoeff(n,k)

m1 = floor(n/2);
m2 = ceil(n/2);

o = [0 cumsum(log(n-k+1)-log(k))];
o(n+1:-1:m1+2) = o(1:m2);

c=sum(o);
