function [pmean,ppred]= lfrpartcoeff(N,k,maxk,mut,muw,minc,maxc)

[W,Ci]=lfrw_mx('N',N,'k',k,'maxk',maxk,'mut',mut,'muw',muw,'minc',minc,'maxc',maxc);
W=W~=0;
ki=sum(W,2);
nc=length(unique(Ci));
P=participation_coef(W,Ci);
pmean=mean(P)
%ppred=1-((1-mut)^2 + (nc-1)*(mut^2))*(1/(N^2));
ppred=1-((1-mut).^2) - (nc-1)*((mut*k)^2)./(ki.^2);
%ppred=ones(N,1)*ppred;
ppredmean=mean(ppred)
plot(1:N,P,'bo',1:N,ppred,'r.');

