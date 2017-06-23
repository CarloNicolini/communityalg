function [stdB,meanB] = distributionzscoremodularity(Aij,ci)

m = number_of_edges(Aij);
n=size(Aij,1);
ki = sum(Aij);
kikj = (ki(:)*ki(:)');
pij = kikj/(2*m-1);
C = bsxfun(@eq,ci,unique(ci)');
Dij = C'*C;
Bij = Aij-pij; % this is the modularity matrix
Min = sum(sum(triu(Dij)))-sum(diag(Dij)); % this is the sum of all possible node pairs within communities
Q = sum(sum((Bij.*Dij)))/(2*m); % standard modularity
stdBij = sqrt(pij.*(1-pij)*Min); % This is the correction on the denominator to standardize edge weights
Qz = sum(sum(Bij./stdBij .*Dij))/(2*m);
Qnewman = sum(sum(((Aij-pij).*Dij)))/(2*m);

stdB=std((Bij(:)./stdBij(:).*Dij(:)));
meanB=mean(( Bij(:)./stdBij(:).*Dij(:) ));