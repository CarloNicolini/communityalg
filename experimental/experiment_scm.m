clear all;
close all;
clc;
addpath(genpath('/data1/workspace/algonet'));
addpath(genpath('/data1/workspace/communityalg'));
addpath(genpath('/data1/workspace/BCT'));

N=100;
C=2;
pout=0.1;
nP = 10;

%AS = [];
%S = [];

for pin=linspace(0,1,nP)
	[A,groups] = randomModularGraphPinPout(N,C,pin,pout);
	memb = group2membership(groups);
	%AS = [S asymptotic_surprisecm(A,memb) ];
	%S = [S asymptotic_surprise(A,memb) ];
	%[~,~,K,~,m,pairs,Bnorm,nc] = comm_mat(A,memb);
	% Questi sono 3 possibili null models che producono diverse funzioni modularità
	%[ sum(nc.*(nc-1)/2)/pairs sum(nc.^2/(N^2)) sum((K./(2*m)).^2)]
	% ma tutti invece hanno la stessa "Surprise" generalizzata
	Wtot = sum(sum(triu(A,1)));
	Pstar = 2*Wtot/(N*(N-1)+2*Wtot);
	%Pstar = density_und(A);
	[predicted_exp_o, predicted_numerator, predicted_denominator] = expected_intracluster_wrg(A,memb);

	obs_numerator = sum(sum(triu(Pstar.^A(1:N/2,1:N/2),1))) + sum(sum(triu(Pstar.^A(N/2+1:N,N/2+1:N),1)));
	obs_denominator = sum(sum(triu(Pstar.^A,1)));
	
	%[ sum(sum(triu(Pstar.^A,1))) 	predicted_denominator ] % GIUSTO
	[ obs_numerator/obs_denominator, obs_numerator*(1-Pstar), obs_denominator*(1-Pstar) ; ...
	  predicted_exp_o, predicted_numerator, predicted_denominator]
end
 

% % 1. Tries to recreate the null model
% k = sum(A);
% kikj2m = k'*k/(2*m);

% deltaCiCj = bsxfun(@eq,memb(:)',unique(memb(:)')')'*bsxfun(@eq,memb(:)',unique(memb(:)')');
% sum(sum(triu(kikj2m.*deltaCiCj)))
% sum(K.^2/(4*m))

% % This is the exact value of modularity
% sum(sum( (A - kikj2m).*deltaCiCj ))/(2*m)
% modularity(A,memb)

% sum(sum( k'*k ))