function [M, Q, res]=correlation_louvain(CIJ,T,M0,ECIJ)
%CORRELATION_LOUVAIN     Optimal community structure from correlation matrix
%   Following the article of MacMahon, Garlaschelli, PhysRev X (2015)
%   Community detection for correlation matrices
%
%   M = correlation_louvain(CIJ,T);
%   [M,Q] = correlation_louvain(CIJ,T);
%   [M,Q] = correlation_louvain(CIJ,T,M0);
%   [M,Q] = correlation_louvain(CIJ,T,M0);
%
%   The optimal community structure is a subdivision of the network into
%   nonoverlapping groups of nodes which maximizes the number of within-
%	group edges, and minimizes the number of between-group edges.
%
%   This function is a fast and accurate multi-iterative generalization of the
%	Louvain community detection algorithm.
%
%   Input:      C       a positive semidefinite correlation matrix
%               T       number of time samples used to build the
%               correlation matrix, (compulsory)
%               M0,    	initial community affiliation vector (optional)
%               B,      objective-function type or custom objective-function matrix (optional)
%                          	'ITNM' Infinite time series with no global mode
%                           'FTNM' Finite time series without global mode
%                           'FTWM' Finite time series with global mode
%
%   Outputs:    M,     	community structure
%               Q,    	optimized quality value
%
%   References: Blondel et al. (2008)  J. Stat. Mech. P10008.
%               Reichardt and Bornholdt (2006) Phys. Rev. E 74, 016110.
%               Ronhovde and Nussinov (2008) Phys. Rev. E 80, 016109
%               Sun et al. (2008) Europhysics Lett 86, 28004.
%               MacMahon et al. (2015) Phys Rev X.
%
%   Basic community_louvain version by "Mika Rubinov, U Cambridge 2015"
%   RMT modification to handle positive/negative correlation matrices by 
%   "Carlo Nicolini, Istituto Italiano di Tecnologia, 2016-2017"
%
% Check positive semidefiniteness by checking the second argument of chol
% (see documentation of chol for details)
[~,p]=chol(CIJ);

if ~isreal(CIJ)
    error('Correlation matrix must be real matrix.');
end

if p~=0
    error('Correlation matrix must be real semipositive definite matrix');
end

n=length(CIJ); % number of nodes
Cnorm=2*sum(sum(triu(CIJ))); % sum of edges (each undirected edge is counted twice)

% Compute the random matrix theory RMT spectrum of CIJ
res = rmtdecompose(CIJ,T);
if ~exist('ECIJ','var') || isempty(ECIJ)
    ECIJ = eye(n); % assume the null model for infinitely long time series without global mode
end

% Check if membership vector already exists
if ~exist('M0','var') || isempty(M0)
    M0=1:n;
elseif numel(M0)~=n
    error('M0 must contain n elements.');
end

% reindex the membership vector
[~,~,Mb] = unique(M0);
M = Mb; % start with a membership vector already reindexed to have values in 1:c where c is the number of communities

if ischar(ECIJ)
    switch ECIJ
        case 'ITNM';
            ECIJ = (CIJ - eye(n)); % Infinite time series without global mode
        case 'FTNM';
            ECIJ = (CIJ-res.Cr); % Finite time series without global mode
        case 'FTWM';
            ECIJ = CIJ - res.Cr - res.Cg; % Finite time series with global mode
        otherwise;
            error('Unknown null model.');
    end
else
    ECIJ = double(ECIJ);
    if ~isequal(size(CIJ),size(ECIJ))
        error('CIJ and ECIJ must have the same size.');
    end
    if max(max(abs(ECIJ-ECIJ.')))>1e-10
        warning('ECIJ is not symmetric, enforcing symmetry.');
    end
end

%ECIJ = (ECIJ+ECIJ.')/2;                                 % symmetrize null model
Hnm=zeros(n,n);                                         % node-to-module degree
for m=1:max(Mb)                                         % loop over modules
    Hnm(:,m)=sum(ECIJ(:,Mb==m),2);
end
H=sum(Hnm,2);                                           % node degree
Hm=sum(Hnm,1);                                          % module degree

Q0 = -inf;
Q = sum(ECIJ(bsxfun(@eq,M0,M0.')));%;/Cnorm             % compute RMT modularity
first_iteration = true;
while Q-Q0>1e-16
    flag = true;                                        % flag for within-hierarchy search
    while flag;
        flag = false;
        for u=randperm(n)                               % loop over all nodes in random order
            ma = Mb(u);                                 % current module of u
            dQ = Hnm(u,:) - Hnm(u,ma) + ECIJ(u,u);
            dQ(ma) = 0;                                 % (line above) algorithm condition
            
            [max_dQ, mb] = max(dQ);                     % maximal increase in modularity and corresponding module
            if max_dQ>1e-10;                            % if maximal increase is positive
                flag = true;
                Mb(u) = mb;                             % reassign module
                Hnm(:,mb) = Hnm(:,mb)+ECIJ(:,u);        % change node-to-module strengths
                Hnm(:,ma) = Hnm(:,ma)-ECIJ(:,u);
                Hm(mb) = Hm(mb)+H(u);                   % change module strengths
                Hm(ma) = Hm(ma)-H(u);
            end
        end
    end
    [~,~,Mb] = unique(Mb);                              % new module assignments
    
    M0 = M;
    if first_iteration
        M=Mb;
        first_iteration=false;
    else
        for u=1:n                                       % loop through initial module assignments
            M(M0==u)=Mb(u);                             % assign new modules
        end
    end
    
    n=max(Mb);                                          % new number of modules
    B1=zeros(n);                                        % new weighted matrix
    for u=1:n
        for v=u:n
            bm=sum(sum(ECIJ(Mb==u,Mb==v)));             % pool weights of nodes in same module
            B1(u,v)=bm;
            B1(v,u)=bm;
        end
    end
    ECIJ=B1/sum(sum(triu(ECIJ)));
    
    Mb=1:n;                                             % initial module assignments
    Hnm=ECIJ;                                           % node-to-module strength
    H=sum(ECIJ);                                        % node strength
    Hm=H;                                               % module strength
    
    Q0=Q;
    Q=trace(ECIJ)/Cnorm;                                % compute modularity
end
Q
