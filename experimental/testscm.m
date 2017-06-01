close all;
clear all;
addpath('/data1/workspace/BCT');
addpath(genpath('/data1/workspace/communityalg'));


% Load the data

A=load('/data1/workspace/communityalg/data/karate.adj');
n = length(A);
m = number_of_edges(A);
memb = community_louvain(A);

qcm=[];
scm=[];
s=[];
qer=[];

reps=250;
for t=1:reps
  memb1 = memb(randperm(n));
  qcm = [qcm modularity(A,memb1)];
  qer = [qer modularity_er(A,memb1)];
  scm = [scm surprisecm(A,memb1)];
  s = [s surprise(A,memb1)];
end

figure;
hold on;
%hist(scm,100);
%hist(qcm.^2*(2*m),100);
plot(1:reps, s,'r',1:reps,2*m*qer.^2);
legend({'Surprise','2m Qer^2'});
title('Surprise VS Modularity ER');
hold off

figure;
hold on;
%hist(scm,100);
%hist(qcm.^2*(2*m),100);
plot(1:reps, scm,'r',1:reps,2*m*qcm.^2);
legend({'SurpriseCM','2m Qcm^2'});
title('SurpriseCM VS Modularity CM');
hold off