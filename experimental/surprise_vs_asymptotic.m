close all;
clear all;
clc;
addpath('/data1/workspace/BCT');
addpath(genpath('/data1/workspace/communityalg'));

addpath('/data1/workspace/algonet/network-toolbox/');

% Load the data
[A,memb]=randomModularGraphPinPout(2000,4,0.5,0.05);
memb=group2membership(memb);
disp('generated graph');

n = length(A);
m = number_of_edges(A);

qcm=[];
scm=[];
s=[];
qer=[];
as = [];
ascm=[];


reps=10;
for t=1:reps
  t
  memb1 = memb(randperm(n));
  qcm = [qcm modularity(A,memb1)];
  qer = [qer modularity_er(A,memb1)];
  scm = [scm surprisecm(A,memb1)];
  s = [s surprise(A,memb1)];
  as = [as asymptotic_surprise(A,memb1)];
  ascm = [ascm asymptotic_surprisecm(A,memb1)];
end

close all;
subplot(1,2,1);
plot(1:reps, sort(s),'r','LineWidth',2,1:reps, sort(as),'b','LineWidth',2);
title('ER');
subplot(1,2,2);
plot(1:reps, sort(scm),'r','LineWidth',2,1:reps, sort(ascm),'b','LineWidth',2);
title('CM');