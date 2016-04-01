% Script to compute the Standardized Mutual Information (SMI) between 
% two clusterings. 
% --------------------------------------------------------------------------
% INPUT: A contingency table T
%        OR
%        Cluster labels of the two clusterings in two vectors
%        eg: true_mem=[1 2 4 1 3 5]
%                 mem=[2 1 3 1 4 5]
%        Cluster labels are coded using positive integers. 
% OUTPUT: SMI 

function [SMI_]=smi(true_mem,mem)
  if nargin==1
    T=true_mem; %contingency table pre-supplied
  elseif nargin==2
    %build the contingency table from membership arrays
    r=max(true_mem);
    c=max(mem);

    %identify & removing the missing labels
    list_t=ismember(1:r,true_mem);
    list_m=ismember(1:c,mem);
    T=Contingency(true_mem,mem);
    T=T(list_t,list_m);
  end

  [r c]=size(T);
  if (c == 1 || r == 1)
   error('Clusterings should have at least 2 clusters')
   return
  end
  
  N = sum(sum(T)); % total number of records
 
  % update the true dimensions
  a=sum(T,2)';
  b=sum(T);

  % compute useful things
  maxNij = min(max(a),max(b));
  NijLogNij=(1:maxNij).*log2(1:maxNij);
  NijLogNij = [0 NijLogNij]; % 0log0 added
  x = -(a(a ~= 0))*log2(a(a ~= 0))' - (b(b ~= 0))*log2(b(b ~= 0))' + N*log2(N);

  % calculate nLogn
  sum_nLogn=0;
  for i=1:r
    for j=1:c
      if T(i,j)>0 
          sum_nLogn = sum_nLogn + NijLogNij(T(i,j)+1);
      end;
    end
  end    
  N_MI = x + sum_nLogn;
  
  % check carefully this ____________________________
  if (N/r/c > 200)
    fprintf('High number of records, N/(rc) > 200, SMI computed using chi-square');
    SMI_ = (2*log(2)*N_MI - (r-1)*(c-1)) / sqrt( 2*(r-1)*(c-1) );
    return
  end
  
  %____________________________
  

  % calculate E[N MI] 
  EP=zeros(r,c);
  for i=1:r
    for j=1:c
      EP(i,j) = E_nLogn(a(i),b(j),N,NijLogNij);
    end
  end

  E_sum_nLogn = sum(sum(EP));
  E_N_MI = x + E_sum_nLogn;

  % calculate E[(N MI)^2] 

  % transpose the contingecy table because of 
  % consideration in Section 3.3 of the paper.
  if (c > r)
    T = T';
    [r c]=size(T);
    a=sum(T,2)';
    b=sum(T);
  end  
  
  % will il take awhile to compute?
  if (r * c * N^3 > 1E7)
    fprintf('Computing MI variance (it might take awhile).');
  else
    fprintf('Computing MI variance.');
  end
  
  EP = zeros(r,c);
  for i=1:r    
    for j=1:c
      fprintf('.'); % just to show it is computing
    
      p = getP(a(i),b(j),N);        
      for nij=max(0,a(i)+b(j)-N):min(a(i), b(j))
        sumP = 0;
        
        % i=i' j=~j' and i=~i' j=~j' 
        % (Lines (3,4) of E[sum_nLogn^2] formula in the Read Me)
        N_ = N - b(j);
        a_ = a(i) - nij;
        for jp=(j+1):c  
          b_ = b(jp);                
          p_= getP(a_,b_,N_);

          for nijp=max(0,a_+b_-N_):min(a_, b_);
            sumP_ = 0;
            for ip=[1:i-1 i+1:r]
              sumP_ = sumP_ + E_nLogn(a(ip), b(jp)-nijp, N-a(i), NijLogNij);
            end
            
            sumP_ = sumP_ + NijLogNij(nijp+1); 
            
            sumP = sumP + sumP_*p_;                  
            p_=incrP(p_,a_,b_,nijp,N_);
          end
        end

        % i=~i' j=j' (Line (2) of E[sum_nLogn^2] formula in the Read Me)
        N_ = N - a(i);
        b_ = b(j) - nij;
        for ip=(i+1):r  
          a_ = a(ip);                
          p_= getP(a_,b_,N_);
          for nipj=max(0,a_+b_-N_):min(a_, b_);
            sumP = sumP + NijLogNij(nipj+1)*p_;
            
            p_=incrP(p_,a_,b_,nipj,N_);
          end
        end            

        % i=i' j=j' (Line (2) of E[sum_nLogn^2] formula in the Read Me)
        sumP = 2*sumP + NijLogNij(nij+1);
        
        Lpnij = NijLogNij(nij+1)*p;
        EP(i,j) = EP(i,j) + Lpnij*sumP;            

        p=incrP(p,a(i),b(j),nij,N);                   
      end
    end
  end

  E_sum_nLogn_2 = sum(sum(EP));    
  E_N_MI_2 = x^2 + 2*x*E_sum_nLogn + E_sum_nLogn_2;

  % Just compute the final value

  SMI_ = (N_MI - E_N_MI)/sqrt(E_N_MI_2 - E_N_MI^2);
   
end
%---------------------auxiliary functions---------------------

% create a contingecy table

function Cont=Contingency(Mem1,Mem2)
  if nargin < 2 || min(size(Mem1)) > 1 || min(size(Mem2)) > 1
     error('Contingency: Requires two vector arguments')
     return
  end

  Cont=zeros(max(Mem1),max(Mem2));

  for i = 1:length(Mem1);
     Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
  end
end

% gets the the probability of the smallest number
% of success for a r.v. Hyp(a,b,N)

function p = getP(a,b,N)
    nij=max(0,a+b-N);
    X=sort([nij N-a-b+nij]);
    if N-b>X(2)
      nom=[[a-nij+1:a] [b-nij+1:b] [X(2)+1:N-b]];
      dem=[[N-a+1:N] [1:X(1)]];
    else
      nom=[[a-nij+1:a] [b-nij+1:b]];       
      dem=[[N-a+1:N] [N-b+1:X(2)] [1:X(1)]];
    end
    p = prod(nom./dem);        
end

% given the probability of n successes for a Hyp(a,b,N)
% computes the probability of n+1 successes

function p = incrP( p,a,b,n,N)
  p = p*(a-n)*(b-n)/(n+1)/(N-a-b+n+1);  
end

% computes E[nLogn] when n is 
% distributed as Hyp(a,b,N)

function [e]= E_nLogn(a,b,N,nLogn)
  p = getP(a,b,N);
  e = 0;
  for n=max(0,a+b-N):min(a,b);
     e = e + nLogn(n+1)*p;
     p = incrP(p,a,b,n,N);
  end
end




            