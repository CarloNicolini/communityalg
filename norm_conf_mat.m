function norm_conf_mat(memb1, memb2)

C = confusionmat(memb1,memb2);
C = normr(C);