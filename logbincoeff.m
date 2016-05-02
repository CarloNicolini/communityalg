function c = logbincoeff(n,k)
c=log(bincoeff(n,k)); % to replace with faster implementation, 
					   % not use nchoosek that works only with scalars
%c=logC(k,n);