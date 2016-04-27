function blockmat = comm_mat(W,ci)
%COMM_MAT comm_mat returns the community matrix, with the links intra and inter communities
% See http://stackoverflow.com/questions/36694384/accelerate-matlab-nested-for-loop-with-bsxfun
% The first matrix C is an ncommsX lenght(ci) matrix of ones and zeros. where each row holds a 1's where ci is equal to one of the unique value
C=bsxfun(@eq,ci,unique(ci)');
%  Now lets say you want to calculate dsm(1,2): you can take the first you can mulitpy W by the first row from the left and the second row from the right. This basically adds (multiply by 1 or zero and then add) all the places in W where both the but only the rows where C(1,:) is equal to one and the colums where C(2,:) is equal to 1. By doing the matrix multiplication, you can achive all at once
blockmat = C*W*C';
blockmat(logical(eye(size(blockmat)))) = blockmat(logical(eye(size(blockmat))))./2;
% The diagonal of comm_links contains the sum of edge weights of edges inside the i-th community, the numbers outside are edges outside community, this is similar to the block matrix