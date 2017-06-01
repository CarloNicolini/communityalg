clear all;
close all;
clc;
addpath(genpath('/data1/workspace/algonet'));
addpath(genpath('/data1/workspace/communityalg'));
addpath(genpath('/data1/workspace/BCT'));

% Calcolo valore atteso ER
%% 1. Fisso la partizione, due blocchi di uguale dimensioen
m = 1;
N = 8;
p = 4/nchoosek(N,2);
memb = [ones(1,N/2) 2*ones(1,N/2)];
% 2. Ora genero tantissimi grafi random e calcolo la percentuale di lati intracluster osservata ed attesa
ntrials=nchoosek(nchoosek(N,2),m)*2000
obs_mzeta 	 = zeros(1,ntrials);
exp_mzeta_cm = zeros(1,ntrials);
exp_mzeta_er  = zeros(1,ntrials);

for trial=1:ntrials
	[B, C, K, ~, ~, pairs, Bnorm, nc] = comm_mat(randomGraph(N,[],m),memb);
	obs_mzeta(trial) = sum(diag(B));
	exp_mzeta_cm(trial) = sum( K.^2/(4*m) );
	exp_mzeta_er(trial) = sum(nc.*(nc-1)/2)*m/pairs;
end

%% 2. Ora questi sono grafi ER con due blocchi di uguali dimensioni
% allora il numero di lati ER expected si può semplicemente scrivere nel modello nullo ER come
[ 2*nchoosek(N/2,2)*p , mean(exp_mzeta_er), std(exp_mzeta_er) ]
% e vedo che questa corrisponde a quella empirica.

% Invece per il valore atteso di mz sotto il CM si ottiene:
[m*(1/(2*m) * N/2 *(N-1)*p) , mean(exp_mzeta_cm), std(exp_mzeta_cm) ]

%% 3. Ora valuto le probabilità che fra tutte queste densità 
% risulti una densità empirica maggiore o uguale di quella expected
mean_er = mean(obs_mzeta >= exp_mzeta_er);
% La deviazione standard inoltre sarà
sqrt(mean_er*(1-mean_er)/(N-1))
% per questo motivo mi aspetto che la Surprise di questa partizione sia 
% molto bassa, faccio ntrials e vedo la Surprise media
allS = arrayfun( @(i)surprise(randomGraph(N,[],m),memb), 1:ntrials );
fprintf('Mean Surprise = 1E-%2.3f\n',mean(allS),'+/- 1E-%2.3f\n,',std(allS))
% questa dovrebbe essere molto vicina a quella calcolata con la probabilità corretta




%[~,sortidx ] = sort(obs_mzeta);
%obs_mzeta = obs_mzeta(sortidx);
%exp_mzeta_cm = exp_mzeta_cm(sortidx);
%exp_mzeta_er = exp_mzeta_er(sortidx);

% plot(1:ntrials,obs_mzeta,'r','LineWidth',2,...
% 	1:ntrials,exp_mzeta_cm,'b','LineWidth',1,...
% 	1:ntrials,exp_mzeta_er,'g','LineWidth',1);
% legend({'Observed','Expected CM','Expected ER'});
% grid;
% % nc=5:25;

% close all;
% %% Here tries to see if Surprise always increases for ring of cliques
% nc=5:50;
% hold on;
% plot(nc, arrayfun( @(i)surprisecm( nthargout([1,2],@ring_of_cliques,i,10){:}), nc ),'r' );
% plot(nc, arrayfun( @(i)asymptotic_surprisecm( nthargout([1,2],@ring_of_cliques,i,10){:})/2, nc ),'b' );
% hold off;
