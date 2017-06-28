function [sol,M,exitflag] = maxent_null_models(method,M)
% Credits to Rossana Mastrandrea, taken from the MAXandSAM package
% http://nl.mathworks.com/matlabcentral/fileexchange/46912-max-sam-package-zip?focused=3825900&tab=function
% Readapted and fixed some typos here and there
% 2017 Carlo Nicolini
% Lorentz Institute, Leiden, The Netherlands
%
% M is the input adjacency matrix
% Procedure for solving the constrained maximization entropy problem for
% all the models when the input data are in matrix-form

options = optimset('Display','off','Algorithm','interior-point','MaxFunEvals',(10^4),'MaxIter',5*10^4,'TolX',10^(-32),'TolFun',10^(-32),'TolCon',10^(-32));

switch upper(method)
    case 'UBCM' 
        n = length(M);
        x0 = 0.9*ones(n,1);
        lb = zeros(n,1);
        ub=Inf*ones(n,1);        
        [sol,fval,exitflag] = fmincon(@(x) Likelihood(x,M,method),x0,[],[],[],[],lb,ub,[],options);

    case 'UWCM' 
        n = length(M);
        x0 = 0.9*ones(n,1);
        lb = zeros(n,1);
        ub=ones(n,1);        
        [sol,fval,exitflag] = fmincon(@(x) Likelihood(x,M,method),x0,[],[],[],[],lb,ub,[],options);
    
    case 'DBCM' 
        n = 2*length(M);
        x0 = 0.9*ones(n,1);
        lb = zeros(n,1);
        ub=Inf*ones(n,1);        
        [sol,fval,exitflag] = fmincon(@(x) Likelihood(x,M,method),x0,[],[],[],[],lb,ub,[],options);
    
    case 'DWCM' 
        n = 2*length(M);
        x0 = 0.9*ones(n,1);
        lb = zeros(n,1);
        ub=[ones(n,1)];        
        [sol,fval,exitflag] = fmincon(@(x) Likelihood(x,M,method),x0,[],[],[],[],lb,ub,[],options);  
        
    case 'UECM'
        n = 2*length(M);
        x0 = 0.9*ones(n,1);
        lb = zeros(n,1);
        ub=[Inf*ones(n/2,1),ones(n/2,1)];        
        [sol,fval,exitflag] = fmincon(@(x) Likelihood(x,M,method),x0,[],[],[],[],lb,ub,[],options); 
    
    case 'RBCM'
        n = 3*length(M);
        x0 = 0.9*ones(n,1);
        lb = zeros(n,1);
        ub=[Inf*ones(n,1)];        
        [sol,fval,exitflag] = fmincon(@(x) Likelihood(x,M,method),x0,[],[],[],[],lb,ub,[],options); 
        
    case 'RWCM'
        n = 3*length(M);
        x0 = 0.9*ones(n,1);
        lb = zeros(n,1);
        ub=[ones(n,1)];
        [sol,fval,exitflag] = fmincon(@(x) Likelihood(x,M,method),x0,[],[],[],[],lb,ub,[],options); 
end

function f = Likelihood(h,M,method)
%Procedure for computing the likelihood function for each model when the
%input are in matrix-form

n = length(h);
f = 0;
switch upper(method)
    case 'UBCM' % Undirected Binary Configuration Model
        %ftest = f + sum(sum(triu(M.*log(h*h') - log(1+h*h'),1)));
        for i = 1:n 
            for j = (i+1):n
                f = f + (M(i,j)*log(h(i)*h(j)) - log(1 + h(i)*h(j)));
            end
        end
       %disp('f-ftest=');
       %[f ftest]
    case 'UWCM'   % Undirected Weighted Configuration Model
        x=h;
        for i=1:n
            for j=(i+1):n
                f=f+(M(i,j)*log(h(i)*h(j))+log(1-h(i)*h(j)));
            end
        end
        
        
    case 'DBCM' % Directed Binary Configuration Model
        m=n/2;
        x=h(1:m);
        y=h((m+1):end);
        
        for i = 1:m 
            for j = 1:m
                if i ~= j
                    f = f + (M(i,j)*log(x(i)*y(j)) - log(1 + x(i)*y(j)));
                end
            end
        end
        
    case 'DWCM'  % Directed Weighted Configuration Model
        m=n/2;
        x=h(1:m);
        y=h((m+1):end);
        
        for i = 1:m 
            for j = 1:m
                if i ~= j
                    f = f + (M(i,j)*log(x(i)*y(j)) + log(1 - x(i)*y(j)));
                end
            end
        end
        
    case 'UECM'  % Undirected Enhanced Configuration Model
        m=n/2;
        x=h(1:m);
        y=h((m+1):end);
        
        A=fix(M>0);
        
        f=0;
        for i=1:m
            for j=(i+1):m
                f=f+A(i,j)*log(x(i)*x(j))+M(i,j)*log(y(i)*y(j))+log(1-y(i)*y(j))-log(1-(y(i)*y(j))+x(i)*x(j)*y(i)*y(j));
            end
        end
        
    case 'RBCM'  % Reciprocal Bianry Configuration Model
        q=n/3;
        x=h(1:q);
        y=h((q+1):2*q);
        z=h((2*q+1):end);
        
       
        [Aout,Ain,Arec]=rec(M);
  
        f=0;
        for i=1:q
            for j=(i+1):q
                f=f+(Aout(i,j)*log(x(i)*y(j))+Ain(i,j)*log(x(j)*y(i))+Arec(i,j)*log(z(i)*z(j))-log(1+x(i)*y(j)+x(j)*y(i)+z(i)*z(j)));
            end
        end
        
    case 'RWCM'  %Reciprocal Weighted Configuration Model        
        q=n/3;
        x=h(1:q);
        y=h((q+1):2*q);
        z=h((2*q+1):end);
        
        [Wout,Win,Wrec]=recW(M);
        
        for i=1:q
            for j=(i+1):q
                Z(i,j) = (1-(x(i)*x(j)*y(i)*y(j)))/((1-(x(i)*y(j)))*(1-(x(j)*y(i)))*(1-(z(i)*z(j))));
                f=f+((Wout(i,j)*log(x(i)*y(j))) + (Win(i,j)*log(x(j)*y(i))) + (Wrec(i,j)*log(z(i)*z(j))) - (log(Z(i,j))));
            end
        end
        
end

f=-f;
