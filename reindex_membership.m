function c = reindex_membership(m)
    m = reshape(m,[1,length(m)]);
    u = unique(m);
    nu = length(u);
    % counts the relative frequencies of the species
    % histcounts(m,length(unique(m)))
    mapsu = [sort(u); 1:nu];
    c = zeros(size(m));
    for i=1:length(m)
        c(i) = mapsu(2,find(m(i)==mapsu(1,:)));
    end
    c = group2membership(sort_group_by_size(membership2groups(c)))';