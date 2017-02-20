function write_sparse_matrix_to_pajek(A,filename)

N=length(A);
[ii,jj,ss] = find(A);

fid = fopen(filename, 'w');
fprintf(fid, '*vertices %6i\n', N);
for i = 1:N
    fprintf(fid, '%d\n',i);
end

fprintf(fid, '*edges\n');
for k=1:length(ii)
    if ii(k)<jj(k)
        fprintf(fid, '%d %d %g\n', ii(k), jj(k), ss(k));
    end
end
fprintf(fid,'\n');
fclose(fid);