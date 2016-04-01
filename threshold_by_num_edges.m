function B = threshold_by_num_edges(A, desidered_num_edges)

n = size(A,1);
d = 2*desidered_num_edges/((n-1)*n);

B = threshold_proportional(A,d);