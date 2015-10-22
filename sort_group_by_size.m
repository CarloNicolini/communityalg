function gs = sort_group_by_size(g)
    [dummy, index] = sort(cellfun('length', g), 'descend');
    gs = g(index);
