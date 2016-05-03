function [Sd] = degree_corrected_sbm(W,ci)
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

[B,C,KTot,n,m,p] = comm_mat(W,ci);
ncomm = size(B,1);

%Sd = 0;
L = 0;
%den = logbincoeff(p,m);
nc = sum(C,2)
for c=1:ncomm
    for d=c:ncomm
        %Sd = Sd + logbincoeff(nc(c)*nc(d),B(c,d))-den;
        ncd = nc(c)*nc(d);
        mcd = B(c,d);
        L = L + mcd*log(mcd/ncd)+(ncd-mcd)*log((ncd-mcd)/ncd);
    end
end
%Sd = Sd*-1;
Sd = L;