clc;
close all;
clear all;
base10=false;
zqsingle=[];
zqpair = [];
qpair = [];
qsingle = [];
asym_surp_single=[];
asym_surp_pair=[];
surp_single=[];
surp_pair=[];
rcliques=4:2:50;
nodes_per_clique=6;
for r=rcliques
	[A,membsingle] = ring_of_cliques(nodes_per_clique,r);
	membpair = membsingle-1.*(mod(membsingle,2)==0);
	% ZScore modularity
	zqsingle = [zqsingle; zscoremodularity(A,membsingle)*(2*number_of_edges(A))];
	zqpair = [zqpair; zscoremodularity(A,membpair)*(2*number_of_edges(A))];
	% Newman Modularity
	qpair = [qpair; modularity(A,membpair)*(2*number_of_edges(A))];
	qsingle = [qsingle; modularity(A,membsingle)*(2*number_of_edges(A))];
	% Asymptotic Surprise
	asym_surp_single = [asym_surp_single; asymptotic_surprise(A,membsingle)];
	asym_surp_pair = [asym_surp_pair; asymptotic_surprise(A,membpair)];
	% Surprise
	surp_single = [surp_single; surprise(A,membsingle,base10)];
	surp_pair = [surp_pair; surprise(A,membpair,base10)];
end

h = figure;
subplot(2,2,1);
plot(rcliques,zqsingle,'r',rcliques,zqpair,'b');
grid;
xlabel('r');
ylabel('Quality');
title('2L*Qzscore');

subplot(2,2,2);
plot(rcliques,qsingle,'r',rcliques,qpair,'b');
grid;
xlabel('r');
ylabel('Quality');
title('2L*QNewman');

subplot(2,2,3);
plot(rcliques,asym_surp_single,'r',rcliques,asym_surp_pair,'b');
grid;
xlabel('r');
ylabel('Quality');
title('Asymptotic Surprise')

subplot(2,2,4);
plot(rcliques,surp_single,'r',rcliques,surp_pair,'b');
grid;
xlabel('r');
ylabel('Quality');
title('Surprise');

print(h,'test.tex','-dpdflatexstandalone');

% %%%%% Compute the zscore modularity for the ring of clique (single) %%%%%%
% clear all;
% clc;
% r=10;
% c=8;
% L = r*c*(c-1)/2+r;
% [A,membsingle] = ring_of_cliques(c,r);
% membpairs = membsingle-1.*(mod(membsingle,2)==0);
% qsingle = (1-2/(c*(c-1)+2))-1/r;
% qpairs = 1-1/(c*(c-1)+2)-2/r;
% qzsingle = qsingle/(sqrt( (c*(c-1)*(c*(c-1)+1)*(r-c*(c-1)-2))/r	))
% qnewmansingle = (r*c*(c-1)/2)/L - ((c*(c-1)+2)/r)/(2*L);

% qzsingle = zscoremodularity(A,membsingle);
% qzpairs = zscoremodularity(A,membpairs);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Make an histogram of how Qz is distributed
% addpath(genpath('/data1/workspace/algonet/'));
% p=0.1;
% memb=[ones(1,50) 2*ones(1,50)];
% zz=[];
% qq=[];
% for g=1:1000
% 	G=randomGraph(100,p);
% 	zz=[zz; zscoremodularity(G,memb)];
% 	qq=[qq; modularity(G,memb)];
% end
% hist(zz)
% hist(qq)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% qz=[];
% for rep=1:5000
% 	qz = [qz; zscoremodularity(A,membsingle(randperm(length(A))))];
% end
% hist(qz)