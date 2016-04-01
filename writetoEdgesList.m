function writetoEdgesList(CIJ, fname, arcs)
%WRITETOPAJ         Write to Pajek
%
%   writetoEdgesList(CIJ, fname, arcs);
%
%   This function writes an edges list .net file from a MATLAB matrix
%
%   Inputs:     CIJ,        adjacency matrix
%               fname,      filename minus .net extension
%               arcs,       1 for directed network
%                           0 for an undirected network
%
%   Chris Honey, Indiana University, 2007


N = size(CIJ,1);
fid = fopen(cat(2,fname,'.net'), 'w');

for i = 1:N
    for j = 1:N
        if CIJ(i,j) ~= 0
            fprintf(fid, '%d %d %g \n', i, j, CIJ(i,j));
        end
    end
end

fclose(fid);