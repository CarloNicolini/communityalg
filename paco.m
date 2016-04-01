function [m,qual]=paco(W)
[m,qual]=paco_mx(W,'quality',2,'nrep',1);
m=reindex_membership(m);
