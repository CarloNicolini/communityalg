function m = group2membership(group)
	n = 0;
	for i=1:length(group)
		n = n+length(group{i});
	end
	
	m = zeros([n,1]);

	for i=1:length(group)
		for j=1:length(group{i})
			m(group{i}(j)) = i;
		end
	end
