function y = kdesmooth(lambda,x,sigma)

y=zeros(1,length(x));
for i=1:length(lambda)
    y = y + 1/(sqrt(2*pi*sigma.^2)).*exp(-((x-lambda(i))/(sqrt(2)*sigma)).^2);
end
%y  = y/sum(y);