[L,ml] = ring_of_custom_cliques([3 3]);
L(1,6)=0;
L(6,1)=0;

degree_corrected_sbm(L,ml)