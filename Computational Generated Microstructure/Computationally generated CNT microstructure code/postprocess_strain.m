load('CNT_Strain_Vf0-5p_L300n_c_2n_li100n_D1n.mat')
strain = -0.006:0.001:0.006;
plot(100*strain,100*mean(result1))
plot(100*strain,100*mean(Rvity))
plot(100*strain,mean(Rvity))