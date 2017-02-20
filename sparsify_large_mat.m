function S = sparsify_large_mat(matobj, fieldname, threshold_percentage)

[nrows, ncols] = size(matobj, fieldname);

if nrows ~= ncols
    error('non square adjacency matrix');
end

S = logical(sparse(0,0));

% Find maximum weight
zmax=-inf;
handle_waitbar = waitbar(0,'Computing maximum edge weight');
for row=1:nrows
    zrow = matobj.(fieldname)(row, row+1:nrows);
    %histogram(zrow);
    zrow(isinf(zrow))=0;
    %error(['contains inf or nan at row ' num2str(row) ]);
    zrowmax = max(zrow);
    if zrowmax > zmax
        zmax = zrowmax;
    end
    waitbar(row/nrows,handle_waitbar,sprintf('Computing maximum edge weight %.1f percent done',row/nrows*100));
end
close(handle_waitbar);

handle_waitbar = waitbar(0,'Generating the sparse graph');
for row=1:nrows
    z = matobj.(fieldname)(row, row+1:nrows);
    [~,s,~] = find(z./zrowmax >= threshold);
    S(row,s)=true;
    S(s,row)=true;
    waitbar(row/nrows,handle_waitbar,sprintf('Generating sparsified graph %.1f percent done',row/nrows*100));
end