function UBCM_Bmat = ubcmmodularity(Aij,ci)
addpath('~/workspace/maxandsam');
xi = matrix_case('UBCM',Aij);
xixj = xi*xi(:)';
Pij = xixj./(1+xixj);
C = bsxfun(@eq,ci, unique(ci)');
Dij = double(C)'*double(C);
Min = sum(sum(triu(Dij)))-sum(diag(Dij)); % this is the sum of all possible node pairs within communities

UBCM_Bmat = (Aij-Pij)./sqrt(Pij.*(1-Pij).*Min);

