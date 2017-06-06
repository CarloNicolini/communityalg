addpath('/data1/workspace/BCT');
addpath('/data1/workspace/communityalg');
A=load('/data1/workspace/communityalg/data/karate.adj');
n = length(A);
m = number_of_edges(A);
memb = community_louvain(A);
zq=[];
scm=[];
s=[];
for t=1:2500
  zq = [zq zmodularity(A,memb(randperm(n)))];
  scm = [scm; surprisecm(A,memb(randperm(n)))];
  s = [s surprise(A,memb)];
end
figure;
hold on;
hist(scm,100);
hist(s,100);
hist(zq,100);
hold off