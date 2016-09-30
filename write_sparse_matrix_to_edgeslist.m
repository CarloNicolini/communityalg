function write_sparse_matrix_to_edgeslist(A,filename)

N=length(A);
[ii,jj,ss] = find(A);

dlmwrite(filename, [ii jj ss], 'delimiter', '\t');
