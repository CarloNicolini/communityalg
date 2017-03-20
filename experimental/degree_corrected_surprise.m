function [Sd] = degree_corrected_surprise(W,ci)
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

[B,C,KTot,~,m] = comm_mat(W,ci);
ncomm = size(B,1);
% number of stubs in comm c is twice number of edges in c 
Kin=2*diag(B)';
sKTot=sum(KTot);
sKin=sum(Kin);

Sd = 0;
den = logbincoeff(sKTot,sKin);
for c=1:ncomm
    Sd = Sd + logbincoeff(KTot(c),Kin(c))-den;
end
Sd = Sd*-1;

SdAsympt = 0;
for c=1:ncomm
	SdAsympt = SdAsympt + Kin(c)/sKin*log((Kin(c)/sKin)/(KTot/sKTot));
end
SdAsympt = SdAsympt*sKin;