function [S,mcd,K] = degree_corrected_surprise(W,ci)
%DEGREE_CORRECTED_SURPRISE

%DEGREE_CORRECTED_SURPRISE      Compute degree corrected surprise of a vertex partition on a binary network.
%
%
%   Inputs      W,  undirected weighted or unweighted network.
%               ci, membership vector
%
%   Outputs:    Q,  value of Degree Corrected Surprise (natural logs).
%               pars, partition parameters. [intraclusteredges,
%               intracluster_pairs, number of edges, number of pairs]
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).
%
n = length(W);
m = sum(nonzeros(triu(W)));
ncomms = length(unique(ci)); % number of communities

mcd = zeros(ncomms);
k=sum(W,1);
K = zeros(ncomms,1);

for c=1:ncomms
    nodesc = find(ci==c);
    K(c) = sum(k(nodesc));
    for d=1:ncomms
        nodesd = find(ci==d);
        M = W(nodesc,nodesd);
        mcd(c,d) = sum(sum(M));
    end
end

% Divide diagonal half because counted twice
mcd(1:ncomms+1:ncomms*ncomms)=mcd(1:ncomms+1:ncomms*ncomms)/2;
mcd

% This is a super fast implementation
C = bsxfun(@eq, ci,unique(ci)');
mcd = C*W*C';
mcd(logical(eye(size(mcd)))) = mcd(logical(eye(size(mcd))))./2;


% The probabilistic definition
denom = bincoeff(4*m,m);
S=1;
for c=1:ncomms
    for d=c:ncomms
        num = K(c) + K(d);
        %fprintf('c=%d d=%d K=%d Kd=%d mcd=%d logC=%g\n',c,d,K(c),K(d),mcd(c,d),logC(mcd(c,d),num));
        S = S + logC(2*mcd(c,d),num) - logC(m,4*m);
    end
end
S=-S;

% % The log definition logS = -log(S)
% logS = 0;
% for c=1:ncomms
%     for d=1:ncomms
%         %logS = logS + logC(K(c)+K(d),mcd(c,d)) - log(denom);
%     end
% end
% logS = -logS;





%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%
    function logH = logHyperProbability(F, M, n, p)
        logH = logC(p,M) + logC(n-p,F-M) - logC(n,F);
        logH = logH/log(10);
    end

    function sum = sumFactorial(n)
        sum = 0;
        if(n > 1)
            for i = 2:n
                sum = sum + log(i);
            end
        end
    end

    function sum = sumRange(min, max)
        sum = 0;
        for i = min:max
            sum = sum + log(i);
        end
    end

    function x = logC(k, n)
        if(k == n || ~k)
            x = 0;
        else
            t = n - k;
            if(t < k)
                t = k;
            end
            x = sumRange(t+1, n) - sumFactorial(n - t);
        end
    end

    function [stop, logP] = sumLogProbabilities(nextLogP, logP)
        if(nextLogP == 0)
            stop = true;
        else
            stop = false;
            if(nextLogP > logP)
                common = nextLogP;
                diffExponent = logP - common;
            else
                common = logP;
                diffExponent = nextLogP - common;
            end
            logP = common + ((log(1 + 10.^diffExponent)) / log(10));
            
            if(nextLogP - logP > -4)
                stop = true;
            end
        end
    end
end

% This function check that the arguments provided to compute_surprise
% function are ok in the sense that they represent a realistic urn model:
% An urn with p balls inside, pi white and p-pi black, and m balls are
% drawm and one looks for the probability that at least mi of them are
% white.
function [argsok] = check_urn_model_validity(p,pi,m,mi)
argsok = true;
if ( mi<0 )
    error('Integer overflow: negative mi');
end

if ( m<0 )
    
    error('Integer overflow: negative m');
end

if ( pi<0 )
    
    error('Integer overflow: negative pi');
end

if ( p<0 )
    error('Integer overflow: negative p');
end

if (pi>p)
    error('Error: impossible urn model p<pi');
end

if (mi>m)
    error('Error: impossible urn model m<mi');
end

if ((m-mi) > (p-pi) )
    
    error('Error: impossible urn model (p-pi)<(m-mi)');
end

if (pi>p)
    error('Error: impossible urn model pi>p');
end

if (mi > pi)
    error('Error: impossible urn modelmi>pi');
end

if (mi==0)
    argsok=false;
    return;
end

if (pi==0)
    argsok=false;
    return;
end

if ((m-mi) == (p-pi) )
    argsok=false;
    return;
end

argsok=true;
end