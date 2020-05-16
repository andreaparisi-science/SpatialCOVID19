#if not available

library(deSolve)



model <- function( t, x, params )  {

	with (as.list(x, params),  {

		ff_0 <- beta/eigen*(KK_0_0*ss_0*(ii_0+aa_0)/NN_0+KK_0_1*ss_0*(ii_1+aa_1)/NN_1+KK_0_2*ss_0*(ii_2+aa_2)/NN_2+KK_0_3*ss_0*(ii_3+aa_3)/NN_3+KK_0_4*ss_0*(ii_4+aa_4)/NN_4+KK_0_5*ss_0*(ii_5+aa_5)/NN_5+KK_0_6*ss_0*(ii_6+aa_6)/NN_6+KK_0_7*ss_0*(ii_7+aa_7)/NN_7+KK_0_8*ss_0*(ii_8+aa_8)/NN_8+KK_0_9*ss_0*(ii_9+aa_9)/NN_9+KK_0_10*ss_0*(ii_10+aa_10)/NN_10+KK_0_11*ss_0*(ii_11+aa_11)/NN_11+KK_0_12*ss_0*(ii_12+aa_12)/NN_12+KK_0_13*ss_0*(ii_13+aa_13)/NN_13+KK_0_14*ss_0*(ii_14+aa_14)/NN_14+KK_0_15*ss_0*(ii_15+aa_15)/NN_15+KK_0_16*ss_0*(ii_16+aa_16)/NN_16)
		ds_0 <- -ff_0
		de_0 <- +ff_0 - sigma*ee_0
		di_0 <- +zz_0*sigma*ee_0 - gamma*ii_0
		da_0 <- +(1.0-zz_0)*sigma*ee_0 - gamma*aa_0
		ff_1 <- beta/eigen*(KK_1_0*ss_1*(ii_0+aa_0)/NN_0+KK_1_1*ss_1*(ii_1+aa_1)/NN_1+KK_1_2*ss_1*(ii_2+aa_2)/NN_2+KK_1_3*ss_1*(ii_3+aa_3)/NN_3+KK_1_4*ss_1*(ii_4+aa_4)/NN_4+KK_1_5*ss_1*(ii_5+aa_5)/NN_5+KK_1_6*ss_1*(ii_6+aa_6)/NN_6+KK_1_7*ss_1*(ii_7+aa_7)/NN_7+KK_1_8*ss_1*(ii_8+aa_8)/NN_8+KK_1_9*ss_1*(ii_9+aa_9)/NN_9+KK_1_10*ss_1*(ii_10+aa_10)/NN_10+KK_1_11*ss_1*(ii_11+aa_11)/NN_11+KK_1_12*ss_1*(ii_12+aa_12)/NN_12+KK_1_13*ss_1*(ii_13+aa_13)/NN_13+KK_1_14*ss_1*(ii_14+aa_14)/NN_14+KK_1_15*ss_1*(ii_15+aa_15)/NN_15+KK_1_16*ss_1*(ii_16+aa_16)/NN_16)
		ds_1 <- -ff_1
		de_1 <- +ff_1 - sigma*ee_1
		di_1 <- +zz_1*sigma*ee_1 - gamma*ii_1
		da_1 <- +(1.0-zz_1)*sigma*ee_1 - gamma*aa_1
		ff_2 <- beta/eigen*(KK_2_0*ss_2*(ii_0+aa_0)/NN_0+KK_2_1*ss_2*(ii_1+aa_1)/NN_1+KK_2_2*ss_2*(ii_2+aa_2)/NN_2+KK_2_3*ss_2*(ii_3+aa_3)/NN_3+KK_2_4*ss_2*(ii_4+aa_4)/NN_4+KK_2_5*ss_2*(ii_5+aa_5)/NN_5+KK_2_6*ss_2*(ii_6+aa_6)/NN_6+KK_2_7*ss_2*(ii_7+aa_7)/NN_7+KK_2_8*ss_2*(ii_8+aa_8)/NN_8+KK_2_9*ss_2*(ii_9+aa_9)/NN_9+KK_2_10*ss_2*(ii_10+aa_10)/NN_10+KK_2_11*ss_2*(ii_11+aa_11)/NN_11+KK_2_12*ss_2*(ii_12+aa_12)/NN_12+KK_2_13*ss_2*(ii_13+aa_13)/NN_13+KK_2_14*ss_2*(ii_14+aa_14)/NN_14+KK_2_15*ss_2*(ii_15+aa_15)/NN_15+KK_2_16*ss_2*(ii_16+aa_16)/NN_16)
		ds_2 <- -ff_2
		de_2 <- +ff_2 - sigma*ee_2
		di_2 <- +zz_2*sigma*ee_2 - gamma*ii_2
		da_2 <- +(1.0-zz_2)*sigma*ee_2 - gamma*aa_2
		ff_3 <- beta/eigen*(KK_3_0*ss_3*(ii_0+aa_0)/NN_0+KK_3_1*ss_3*(ii_1+aa_1)/NN_1+KK_3_2*ss_3*(ii_2+aa_2)/NN_2+KK_3_3*ss_3*(ii_3+aa_3)/NN_3+KK_3_4*ss_3*(ii_4+aa_4)/NN_4+KK_3_5*ss_3*(ii_5+aa_5)/NN_5+KK_3_6*ss_3*(ii_6+aa_6)/NN_6+KK_3_7*ss_3*(ii_7+aa_7)/NN_7+KK_3_8*ss_3*(ii_8+aa_8)/NN_8+KK_3_9*ss_3*(ii_9+aa_9)/NN_9+KK_3_10*ss_3*(ii_10+aa_10)/NN_10+KK_3_11*ss_3*(ii_11+aa_11)/NN_11+KK_3_12*ss_3*(ii_12+aa_12)/NN_12+KK_3_13*ss_3*(ii_13+aa_13)/NN_13+KK_3_14*ss_3*(ii_14+aa_14)/NN_14+KK_3_15*ss_3*(ii_15+aa_15)/NN_15+KK_3_16*ss_3*(ii_16+aa_16)/NN_16)
		ds_3 <- -ff_3
		de_3 <- +ff_3 - sigma*ee_3
		di_3 <- +zz_3*sigma*ee_3 - gamma*ii_3
		da_3 <- +(1.0-zz_3)*sigma*ee_3 - gamma*aa_3
		ff_4 <- beta/eigen*(KK_4_0*ss_4*(ii_0+aa_0)/NN_0+KK_4_1*ss_4*(ii_1+aa_1)/NN_1+KK_4_2*ss_4*(ii_2+aa_2)/NN_2+KK_4_3*ss_4*(ii_3+aa_3)/NN_3+KK_4_4*ss_4*(ii_4+aa_4)/NN_4+KK_4_5*ss_4*(ii_5+aa_5)/NN_5+KK_4_6*ss_4*(ii_6+aa_6)/NN_6+KK_4_7*ss_4*(ii_7+aa_7)/NN_7+KK_4_8*ss_4*(ii_8+aa_8)/NN_8+KK_4_9*ss_4*(ii_9+aa_9)/NN_9+KK_4_10*ss_4*(ii_10+aa_10)/NN_10+KK_4_11*ss_4*(ii_11+aa_11)/NN_11+KK_4_12*ss_4*(ii_12+aa_12)/NN_12+KK_4_13*ss_4*(ii_13+aa_13)/NN_13+KK_4_14*ss_4*(ii_14+aa_14)/NN_14+KK_4_15*ss_4*(ii_15+aa_15)/NN_15+KK_4_16*ss_4*(ii_16+aa_16)/NN_16)
		ds_4 <- -ff_4
		de_4 <- +ff_4 - sigma*ee_4
		di_4 <- +zz_4*sigma*ee_4 - gamma*ii_4
		da_4 <- +(1.0-zz_4)*sigma*ee_4 - gamma*aa_4
		ff_5 <- beta/eigen*(KK_5_0*ss_5*(ii_0+aa_0)/NN_0+KK_5_1*ss_5*(ii_1+aa_1)/NN_1+KK_5_2*ss_5*(ii_2+aa_2)/NN_2+KK_5_3*ss_5*(ii_3+aa_3)/NN_3+KK_5_4*ss_5*(ii_4+aa_4)/NN_4+KK_5_5*ss_5*(ii_5+aa_5)/NN_5+KK_5_6*ss_5*(ii_6+aa_6)/NN_6+KK_5_7*ss_5*(ii_7+aa_7)/NN_7+KK_5_8*ss_5*(ii_8+aa_8)/NN_8+KK_5_9*ss_5*(ii_9+aa_9)/NN_9+KK_5_10*ss_5*(ii_10+aa_10)/NN_10+KK_5_11*ss_5*(ii_11+aa_11)/NN_11+KK_5_12*ss_5*(ii_12+aa_12)/NN_12+KK_5_13*ss_5*(ii_13+aa_13)/NN_13+KK_5_14*ss_5*(ii_14+aa_14)/NN_14+KK_5_15*ss_5*(ii_15+aa_15)/NN_15+KK_5_16*ss_5*(ii_16+aa_16)/NN_16)
		ds_5 <- -ff_5
		de_5 <- +ff_5 - sigma*ee_5
		di_5 <- +zz_5*sigma*ee_5 - gamma*ii_5
		da_5 <- +(1.0-zz_5)*sigma*ee_5 - gamma*aa_5
		ff_6 <- beta/eigen*(KK_6_0*ss_6*(ii_0+aa_0)/NN_0+KK_6_1*ss_6*(ii_1+aa_1)/NN_1+KK_6_2*ss_6*(ii_2+aa_2)/NN_2+KK_6_3*ss_6*(ii_3+aa_3)/NN_3+KK_6_4*ss_6*(ii_4+aa_4)/NN_4+KK_6_5*ss_6*(ii_5+aa_5)/NN_5+KK_6_6*ss_6*(ii_6+aa_6)/NN_6+KK_6_7*ss_6*(ii_7+aa_7)/NN_7+KK_6_8*ss_6*(ii_8+aa_8)/NN_8+KK_6_9*ss_6*(ii_9+aa_9)/NN_9+KK_6_10*ss_6*(ii_10+aa_10)/NN_10+KK_6_11*ss_6*(ii_11+aa_11)/NN_11+KK_6_12*ss_6*(ii_12+aa_12)/NN_12+KK_6_13*ss_6*(ii_13+aa_13)/NN_13+KK_6_14*ss_6*(ii_14+aa_14)/NN_14+KK_6_15*ss_6*(ii_15+aa_15)/NN_15+KK_6_16*ss_6*(ii_16+aa_16)/NN_16)
		ds_6 <- -ff_6
		de_6 <- +ff_6 - sigma*ee_6
		di_6 <- +zz_6*sigma*ee_6 - gamma*ii_6
		da_6 <- +(1.0-zz_6)*sigma*ee_6 - gamma*aa_6
		ff_7 <- beta/eigen*(KK_7_0*ss_7*(ii_0+aa_0)/NN_0+KK_7_1*ss_7*(ii_1+aa_1)/NN_1+KK_7_2*ss_7*(ii_2+aa_2)/NN_2+KK_7_3*ss_7*(ii_3+aa_3)/NN_3+KK_7_4*ss_7*(ii_4+aa_4)/NN_4+KK_7_5*ss_7*(ii_5+aa_5)/NN_5+KK_7_6*ss_7*(ii_6+aa_6)/NN_6+KK_7_7*ss_7*(ii_7+aa_7)/NN_7+KK_7_8*ss_7*(ii_8+aa_8)/NN_8+KK_7_9*ss_7*(ii_9+aa_9)/NN_9+KK_7_10*ss_7*(ii_10+aa_10)/NN_10+KK_7_11*ss_7*(ii_11+aa_11)/NN_11+KK_7_12*ss_7*(ii_12+aa_12)/NN_12+KK_7_13*ss_7*(ii_13+aa_13)/NN_13+KK_7_14*ss_7*(ii_14+aa_14)/NN_14+KK_7_15*ss_7*(ii_15+aa_15)/NN_15+KK_7_16*ss_7*(ii_16+aa_16)/NN_16)
		ds_7 <- -ff_7
		de_7 <- +ff_7 - sigma*ee_7
		di_7 <- +zz_7*sigma*ee_7 - gamma*ii_7
		da_7 <- +(1.0-zz_7)*sigma*ee_7 - gamma*aa_7
		ff_8 <- beta/eigen*(KK_8_0*ss_8*(ii_0+aa_0)/NN_0+KK_8_1*ss_8*(ii_1+aa_1)/NN_1+KK_8_2*ss_8*(ii_2+aa_2)/NN_2+KK_8_3*ss_8*(ii_3+aa_3)/NN_3+KK_8_4*ss_8*(ii_4+aa_4)/NN_4+KK_8_5*ss_8*(ii_5+aa_5)/NN_5+KK_8_6*ss_8*(ii_6+aa_6)/NN_6+KK_8_7*ss_8*(ii_7+aa_7)/NN_7+KK_8_8*ss_8*(ii_8+aa_8)/NN_8+KK_8_9*ss_8*(ii_9+aa_9)/NN_9+KK_8_10*ss_8*(ii_10+aa_10)/NN_10+KK_8_11*ss_8*(ii_11+aa_11)/NN_11+KK_8_12*ss_8*(ii_12+aa_12)/NN_12+KK_8_13*ss_8*(ii_13+aa_13)/NN_13+KK_8_14*ss_8*(ii_14+aa_14)/NN_14+KK_8_15*ss_8*(ii_15+aa_15)/NN_15+KK_8_16*ss_8*(ii_16+aa_16)/NN_16)
		ds_8 <- -ff_8
		de_8 <- +ff_8 - sigma*ee_8
		di_8 <- +zz_8*sigma*ee_8 - gamma*ii_8
		da_8 <- +(1.0-zz_8)*sigma*ee_8 - gamma*aa_8
		ff_9 <- beta/eigen*(KK_9_0*ss_9*(ii_0+aa_0)/NN_0+KK_9_1*ss_9*(ii_1+aa_1)/NN_1+KK_9_2*ss_9*(ii_2+aa_2)/NN_2+KK_9_3*ss_9*(ii_3+aa_3)/NN_3+KK_9_4*ss_9*(ii_4+aa_4)/NN_4+KK_9_5*ss_9*(ii_5+aa_5)/NN_5+KK_9_6*ss_9*(ii_6+aa_6)/NN_6+KK_9_7*ss_9*(ii_7+aa_7)/NN_7+KK_9_8*ss_9*(ii_8+aa_8)/NN_8+KK_9_9*ss_9*(ii_9+aa_9)/NN_9+KK_9_10*ss_9*(ii_10+aa_10)/NN_10+KK_9_11*ss_9*(ii_11+aa_11)/NN_11+KK_9_12*ss_9*(ii_12+aa_12)/NN_12+KK_9_13*ss_9*(ii_13+aa_13)/NN_13+KK_9_14*ss_9*(ii_14+aa_14)/NN_14+KK_9_15*ss_9*(ii_15+aa_15)/NN_15+KK_9_16*ss_9*(ii_16+aa_16)/NN_16)
		ds_9 <- -ff_9
		de_9 <- +ff_9 - sigma*ee_9
		di_9 <- +zz_9*sigma*ee_9 - gamma*ii_9
		da_9 <- +(1.0-zz_9)*sigma*ee_9 - gamma*aa_9
		ff_10 <- beta/eigen*(KK_10_0*ss_10*(ii_0+aa_0)/NN_0+KK_10_1*ss_10*(ii_1+aa_1)/NN_1+KK_10_2*ss_10*(ii_2+aa_2)/NN_2+KK_10_3*ss_10*(ii_3+aa_3)/NN_3+KK_10_4*ss_10*(ii_4+aa_4)/NN_4+KK_10_5*ss_10*(ii_5+aa_5)/NN_5+KK_10_6*ss_10*(ii_6+aa_6)/NN_6+KK_10_7*ss_10*(ii_7+aa_7)/NN_7+KK_10_8*ss_10*(ii_8+aa_8)/NN_8+KK_10_9*ss_10*(ii_9+aa_9)/NN_9+KK_10_10*ss_10*(ii_10+aa_10)/NN_10+KK_10_11*ss_10*(ii_11+aa_11)/NN_11+KK_10_12*ss_10*(ii_12+aa_12)/NN_12+KK_10_13*ss_10*(ii_13+aa_13)/NN_13+KK_10_14*ss_10*(ii_14+aa_14)/NN_14+KK_10_15*ss_10*(ii_15+aa_15)/NN_15+KK_10_16*ss_10*(ii_16+aa_16)/NN_16)
		ds_10 <- -ff_10
		de_10 <- +ff_10 - sigma*ee_10
		di_10 <- +zz_10*sigma*ee_10 - gamma*ii_10
		da_10 <- +(1.0-zz_10)*sigma*ee_10 - gamma*aa_10
		ff_11 <- beta/eigen*(KK_11_0*ss_11*(ii_0+aa_0)/NN_0+KK_11_1*ss_11*(ii_1+aa_1)/NN_1+KK_11_2*ss_11*(ii_2+aa_2)/NN_2+KK_11_3*ss_11*(ii_3+aa_3)/NN_3+KK_11_4*ss_11*(ii_4+aa_4)/NN_4+KK_11_5*ss_11*(ii_5+aa_5)/NN_5+KK_11_6*ss_11*(ii_6+aa_6)/NN_6+KK_11_7*ss_11*(ii_7+aa_7)/NN_7+KK_11_8*ss_11*(ii_8+aa_8)/NN_8+KK_11_9*ss_11*(ii_9+aa_9)/NN_9+KK_11_10*ss_11*(ii_10+aa_10)/NN_10+KK_11_11*ss_11*(ii_11+aa_11)/NN_11+KK_11_12*ss_11*(ii_12+aa_12)/NN_12+KK_11_13*ss_11*(ii_13+aa_13)/NN_13+KK_11_14*ss_11*(ii_14+aa_14)/NN_14+KK_11_15*ss_11*(ii_15+aa_15)/NN_15+KK_11_16*ss_11*(ii_16+aa_16)/NN_16)
		ds_11 <- -ff_11
		de_11 <- +ff_11 - sigma*ee_11
		di_11 <- +zz_11*sigma*ee_11 - gamma*ii_11
		da_11 <- +(1.0-zz_11)*sigma*ee_11 - gamma*aa_11
		ff_12 <- beta/eigen*(KK_12_0*ss_12*(ii_0+aa_0)/NN_0+KK_12_1*ss_12*(ii_1+aa_1)/NN_1+KK_12_2*ss_12*(ii_2+aa_2)/NN_2+KK_12_3*ss_12*(ii_3+aa_3)/NN_3+KK_12_4*ss_12*(ii_4+aa_4)/NN_4+KK_12_5*ss_12*(ii_5+aa_5)/NN_5+KK_12_6*ss_12*(ii_6+aa_6)/NN_6+KK_12_7*ss_12*(ii_7+aa_7)/NN_7+KK_12_8*ss_12*(ii_8+aa_8)/NN_8+KK_12_9*ss_12*(ii_9+aa_9)/NN_9+KK_12_10*ss_12*(ii_10+aa_10)/NN_10+KK_12_11*ss_12*(ii_11+aa_11)/NN_11+KK_12_12*ss_12*(ii_12+aa_12)/NN_12+KK_12_13*ss_12*(ii_13+aa_13)/NN_13+KK_12_14*ss_12*(ii_14+aa_14)/NN_14+KK_12_15*ss_12*(ii_15+aa_15)/NN_15+KK_12_16*ss_12*(ii_16+aa_16)/NN_16)
		ds_12 <- -ff_12
		de_12 <- +ff_12 - sigma*ee_12
		di_12 <- +zz_12*sigma*ee_12 - gamma*ii_12
		da_12 <- +(1.0-zz_12)*sigma*ee_12 - gamma*aa_12
		ff_13 <- beta/eigen*(KK_13_0*ss_13*(ii_0+aa_0)/NN_0+KK_13_1*ss_13*(ii_1+aa_1)/NN_1+KK_13_2*ss_13*(ii_2+aa_2)/NN_2+KK_13_3*ss_13*(ii_3+aa_3)/NN_3+KK_13_4*ss_13*(ii_4+aa_4)/NN_4+KK_13_5*ss_13*(ii_5+aa_5)/NN_5+KK_13_6*ss_13*(ii_6+aa_6)/NN_6+KK_13_7*ss_13*(ii_7+aa_7)/NN_7+KK_13_8*ss_13*(ii_8+aa_8)/NN_8+KK_13_9*ss_13*(ii_9+aa_9)/NN_9+KK_13_10*ss_13*(ii_10+aa_10)/NN_10+KK_13_11*ss_13*(ii_11+aa_11)/NN_11+KK_13_12*ss_13*(ii_12+aa_12)/NN_12+KK_13_13*ss_13*(ii_13+aa_13)/NN_13+KK_13_14*ss_13*(ii_14+aa_14)/NN_14+KK_13_15*ss_13*(ii_15+aa_15)/NN_15+KK_13_16*ss_13*(ii_16+aa_16)/NN_16)
		ds_13 <- -ff_13
		de_13 <- +ff_13 - sigma*ee_13
		di_13 <- +zz_13*sigma*ee_13 - gamma*ii_13
		da_13 <- +(1.0-zz_13)*sigma*ee_13 - gamma*aa_13
		ff_14 <- beta/eigen*(KK_14_0*ss_14*(ii_0+aa_0)/NN_0+KK_14_1*ss_14*(ii_1+aa_1)/NN_1+KK_14_2*ss_14*(ii_2+aa_2)/NN_2+KK_14_3*ss_14*(ii_3+aa_3)/NN_3+KK_14_4*ss_14*(ii_4+aa_4)/NN_4+KK_14_5*ss_14*(ii_5+aa_5)/NN_5+KK_14_6*ss_14*(ii_6+aa_6)/NN_6+KK_14_7*ss_14*(ii_7+aa_7)/NN_7+KK_14_8*ss_14*(ii_8+aa_8)/NN_8+KK_14_9*ss_14*(ii_9+aa_9)/NN_9+KK_14_10*ss_14*(ii_10+aa_10)/NN_10+KK_14_11*ss_14*(ii_11+aa_11)/NN_11+KK_14_12*ss_14*(ii_12+aa_12)/NN_12+KK_14_13*ss_14*(ii_13+aa_13)/NN_13+KK_14_14*ss_14*(ii_14+aa_14)/NN_14+KK_14_15*ss_14*(ii_15+aa_15)/NN_15+KK_14_16*ss_14*(ii_16+aa_16)/NN_16)
		ds_14 <- -ff_14
		de_14 <- +ff_14 - sigma*ee_14
		di_14 <- +zz_14*sigma*ee_14 - gamma*ii_14
		da_14 <- +(1.0-zz_14)*sigma*ee_14 - gamma*aa_14
		ff_15 <- beta/eigen*(KK_15_0*ss_15*(ii_0+aa_0)/NN_0+KK_15_1*ss_15*(ii_1+aa_1)/NN_1+KK_15_2*ss_15*(ii_2+aa_2)/NN_2+KK_15_3*ss_15*(ii_3+aa_3)/NN_3+KK_15_4*ss_15*(ii_4+aa_4)/NN_4+KK_15_5*ss_15*(ii_5+aa_5)/NN_5+KK_15_6*ss_15*(ii_6+aa_6)/NN_6+KK_15_7*ss_15*(ii_7+aa_7)/NN_7+KK_15_8*ss_15*(ii_8+aa_8)/NN_8+KK_15_9*ss_15*(ii_9+aa_9)/NN_9+KK_15_10*ss_15*(ii_10+aa_10)/NN_10+KK_15_11*ss_15*(ii_11+aa_11)/NN_11+KK_15_12*ss_15*(ii_12+aa_12)/NN_12+KK_15_13*ss_15*(ii_13+aa_13)/NN_13+KK_15_14*ss_15*(ii_14+aa_14)/NN_14+KK_15_15*ss_15*(ii_15+aa_15)/NN_15+KK_15_16*ss_15*(ii_16+aa_16)/NN_16)
		ds_15 <- -ff_15
		de_15 <- +ff_15 - sigma*ee_15
		di_15 <- +zz_15*sigma*ee_15 - gamma*ii_15
		da_15 <- +(1.0-zz_15)*sigma*ee_15 - gamma*aa_15
		ff_16 <- beta/eigen*(KK_16_0*ss_16*(ii_0+aa_0)/NN_0+KK_16_1*ss_16*(ii_1+aa_1)/NN_1+KK_16_2*ss_16*(ii_2+aa_2)/NN_2+KK_16_3*ss_16*(ii_3+aa_3)/NN_3+KK_16_4*ss_16*(ii_4+aa_4)/NN_4+KK_16_5*ss_16*(ii_5+aa_5)/NN_5+KK_16_6*ss_16*(ii_6+aa_6)/NN_6+KK_16_7*ss_16*(ii_7+aa_7)/NN_7+KK_16_8*ss_16*(ii_8+aa_8)/NN_8+KK_16_9*ss_16*(ii_9+aa_9)/NN_9+KK_16_10*ss_16*(ii_10+aa_10)/NN_10+KK_16_11*ss_16*(ii_11+aa_11)/NN_11+KK_16_12*ss_16*(ii_12+aa_12)/NN_12+KK_16_13*ss_16*(ii_13+aa_13)/NN_13+KK_16_14*ss_16*(ii_14+aa_14)/NN_14+KK_16_15*ss_16*(ii_15+aa_15)/NN_15+KK_16_16*ss_16*(ii_16+aa_16)/NN_16)
		ds_16 <- -ff_16
		de_16 <- +ff_16 - sigma*ee_16
		di_16 <- +zz_16*sigma*ee_16 - gamma*ii_16
		da_16 <- +(1.0-zz_16)*sigma*ee_16 - gamma*aa_16
		list( c(ds_0, de_0, di_0, da_0, ds_1, de_1, di_1, da_1, ds_2, de_2, di_2, da_2, ds_3, de_3, di_3, da_3, ds_4, de_4, di_4, da_4, ds_5, de_5, di_5, da_5, ds_6, de_6, di_6, da_6, ds_7, de_7, di_7, da_7, ds_8, de_8, di_8, da_8, ds_9, de_9, di_9, da_9, ds_10, de_10, di_10, da_10, ds_11, de_11, di_11, da_11, ds_12, de_12, di_12, da_12, ds_13, de_13, di_13, da_13, ds_14, de_14, di_14, da_14, ds_15, de_15, di_15, da_15, ds_16, de_16, di_16, da_16) )
	})
}
dat <- read.table( '../../../Test/Contacts/KenyaContactMatrix.csv', header=FALSE )
KK <- data.matrix(dat)
KK <- cbind(KK, KK[,16])
KK <- rbind(KK, KK[16,])
zz = c( 0.0014432192790456596, 0.0009357685327922587, 0.0010761342759786922, 0.0013458086273749326, 
     0.013695268824959014, 0.016889531485125703, 0.03167807081661833, 0.024290277337257693, 0.028326981798033054, 
     0.040739451216136126, 0.07430476961279056, 0.13743204603117717, 0.17923146887955166, 0.3315015659109177, 
     0.2923192122329743, 0.375935907058367, 1.00)
zz_0 = zz[1]
zz_1 = zz[2]
zz_2 = zz[3]
zz_3 = zz[4]
zz_4 = zz[5]
zz_5 = zz[6]
zz_6 = zz[7]
zz_7 = zz[8]
zz_8 = zz[9]
zz_9 = zz[10]
zz_10 = zz[11]
zz_11 = zz[12]
zz_12 = zz[13]
zz_13 = zz[14]
zz_14 = zz[15]
zz_15 = zz[16]
zz_16 = zz[17]
dat <- read.table( '../../../Test/Setup/Test_5000km_stats.dat', header=FALSE)
NN <- 47983657
NN_0 <- dat$V3[1]
NN_1 <- dat$V3[2]
NN_2 <- dat$V3[3]
NN_3 <- dat$V3[4]
NN_4 <- dat$V3[5]
NN_5 <- dat$V3[6]
NN_6 <- dat$V3[7]
NN_7 <- dat$V3[8]
NN_8 <- dat$V3[9]
NN_9 <- dat$V3[10]
NN_10 <- dat$V3[11]
NN_11 <- dat$V3[12]
NN_12 <- dat$V3[13]
NN_13 <- dat$V3[14]
NN_14 <- dat$V3[15]
NN_15 <- dat$V3[16]
NN_16 <- dat$V3[17]
R0 <- 2.5
sigma <- 1/3.0
gamma <- 1/4.0
eigen <- 16.7448

beta <- R0*gamma
KK_0_0=KK[1, 1]
KK_0_1=KK[1, 2]
KK_0_2=KK[1, 3]
KK_0_3=KK[1, 4]
KK_0_4=KK[1, 5]
KK_0_5=KK[1, 6]
KK_0_6=KK[1, 7]
KK_0_7=KK[1, 8]
KK_0_8=KK[1, 9]
KK_0_9=KK[1, 10]
KK_0_10=KK[1, 11]
KK_0_11=KK[1, 12]
KK_0_12=KK[1, 13]
KK_0_13=KK[1, 14]
KK_0_14=KK[1, 15]
KK_0_15=KK[1, 16]
KK_0_16=KK[1, 17]
KK_1_0=KK[2, 1]
KK_1_1=KK[2, 2]
KK_1_2=KK[2, 3]
KK_1_3=KK[2, 4]
KK_1_4=KK[2, 5]
KK_1_5=KK[2, 6]
KK_1_6=KK[2, 7]
KK_1_7=KK[2, 8]
KK_1_8=KK[2, 9]
KK_1_9=KK[2, 10]
KK_1_10=KK[2, 11]
KK_1_11=KK[2, 12]
KK_1_12=KK[2, 13]
KK_1_13=KK[2, 14]
KK_1_14=KK[2, 15]
KK_1_15=KK[2, 16]
KK_1_16=KK[2, 17]
KK_2_0=KK[3, 1]
KK_2_1=KK[3, 2]
KK_2_2=KK[3, 3]
KK_2_3=KK[3, 4]
KK_2_4=KK[3, 5]
KK_2_5=KK[3, 6]
KK_2_6=KK[3, 7]
KK_2_7=KK[3, 8]
KK_2_8=KK[3, 9]
KK_2_9=KK[3, 10]
KK_2_10=KK[3, 11]
KK_2_11=KK[3, 12]
KK_2_12=KK[3, 13]
KK_2_13=KK[3, 14]
KK_2_14=KK[3, 15]
KK_2_15=KK[3, 16]
KK_2_16=KK[3, 17]
KK_3_0=KK[4, 1]
KK_3_1=KK[4, 2]
KK_3_2=KK[4, 3]
KK_3_3=KK[4, 4]
KK_3_4=KK[4, 5]
KK_3_5=KK[4, 6]
KK_3_6=KK[4, 7]
KK_3_7=KK[4, 8]
KK_3_8=KK[4, 9]
KK_3_9=KK[4, 10]
KK_3_10=KK[4, 11]
KK_3_11=KK[4, 12]
KK_3_12=KK[4, 13]
KK_3_13=KK[4, 14]
KK_3_14=KK[4, 15]
KK_3_15=KK[4, 16]
KK_3_16=KK[4, 17]
KK_4_0=KK[5, 1]
KK_4_1=KK[5, 2]
KK_4_2=KK[5, 3]
KK_4_3=KK[5, 4]
KK_4_4=KK[5, 5]
KK_4_5=KK[5, 6]
KK_4_6=KK[5, 7]
KK_4_7=KK[5, 8]
KK_4_8=KK[5, 9]
KK_4_9=KK[5, 10]
KK_4_10=KK[5, 11]
KK_4_11=KK[5, 12]
KK_4_12=KK[5, 13]
KK_4_13=KK[5, 14]
KK_4_14=KK[5, 15]
KK_4_15=KK[5, 16]
KK_4_16=KK[5, 17]
KK_5_0=KK[6, 1]
KK_5_1=KK[6, 2]
KK_5_2=KK[6, 3]
KK_5_3=KK[6, 4]
KK_5_4=KK[6, 5]
KK_5_5=KK[6, 6]
KK_5_6=KK[6, 7]
KK_5_7=KK[6, 8]
KK_5_8=KK[6, 9]
KK_5_9=KK[6, 10]
KK_5_10=KK[6, 11]
KK_5_11=KK[6, 12]
KK_5_12=KK[6, 13]
KK_5_13=KK[6, 14]
KK_5_14=KK[6, 15]
KK_5_15=KK[6, 16]
KK_5_16=KK[6, 17]
KK_6_0=KK[7, 1]
KK_6_1=KK[7, 2]
KK_6_2=KK[7, 3]
KK_6_3=KK[7, 4]
KK_6_4=KK[7, 5]
KK_6_5=KK[7, 6]
KK_6_6=KK[7, 7]
KK_6_7=KK[7, 8]
KK_6_8=KK[7, 9]
KK_6_9=KK[7, 10]
KK_6_10=KK[7, 11]
KK_6_11=KK[7, 12]
KK_6_12=KK[7, 13]
KK_6_13=KK[7, 14]
KK_6_14=KK[7, 15]
KK_6_15=KK[7, 16]
KK_6_16=KK[7, 17]
KK_7_0=KK[8, 1]
KK_7_1=KK[8, 2]
KK_7_2=KK[8, 3]
KK_7_3=KK[8, 4]
KK_7_4=KK[8, 5]
KK_7_5=KK[8, 6]
KK_7_6=KK[8, 7]
KK_7_7=KK[8, 8]
KK_7_8=KK[8, 9]
KK_7_9=KK[8, 10]
KK_7_10=KK[8, 11]
KK_7_11=KK[8, 12]
KK_7_12=KK[8, 13]
KK_7_13=KK[8, 14]
KK_7_14=KK[8, 15]
KK_7_15=KK[8, 16]
KK_7_16=KK[8, 17]
KK_8_0=KK[9, 1]
KK_8_1=KK[9, 2]
KK_8_2=KK[9, 3]
KK_8_3=KK[9, 4]
KK_8_4=KK[9, 5]
KK_8_5=KK[9, 6]
KK_8_6=KK[9, 7]
KK_8_7=KK[9, 8]
KK_8_8=KK[9, 9]
KK_8_9=KK[9, 10]
KK_8_10=KK[9, 11]
KK_8_11=KK[9, 12]
KK_8_12=KK[9, 13]
KK_8_13=KK[9, 14]
KK_8_14=KK[9, 15]
KK_8_15=KK[9, 16]
KK_8_16=KK[9, 17]
KK_9_0=KK[10, 1]
KK_9_1=KK[10, 2]
KK_9_2=KK[10, 3]
KK_9_3=KK[10, 4]
KK_9_4=KK[10, 5]
KK_9_5=KK[10, 6]
KK_9_6=KK[10, 7]
KK_9_7=KK[10, 8]
KK_9_8=KK[10, 9]
KK_9_9=KK[10, 10]
KK_9_10=KK[10, 11]
KK_9_11=KK[10, 12]
KK_9_12=KK[10, 13]
KK_9_13=KK[10, 14]
KK_9_14=KK[10, 15]
KK_9_15=KK[10, 16]
KK_9_16=KK[10, 17]
KK_10_0=KK[11, 1]
KK_10_1=KK[11, 2]
KK_10_2=KK[11, 3]
KK_10_3=KK[11, 4]
KK_10_4=KK[11, 5]
KK_10_5=KK[11, 6]
KK_10_6=KK[11, 7]
KK_10_7=KK[11, 8]
KK_10_8=KK[11, 9]
KK_10_9=KK[11, 10]
KK_10_10=KK[11, 11]
KK_10_11=KK[11, 12]
KK_10_12=KK[11, 13]
KK_10_13=KK[11, 14]
KK_10_14=KK[11, 15]
KK_10_15=KK[11, 16]
KK_10_16=KK[11, 17]
KK_11_0=KK[12, 1]
KK_11_1=KK[12, 2]
KK_11_2=KK[12, 3]
KK_11_3=KK[12, 4]
KK_11_4=KK[12, 5]
KK_11_5=KK[12, 6]
KK_11_6=KK[12, 7]
KK_11_7=KK[12, 8]
KK_11_8=KK[12, 9]
KK_11_9=KK[12, 10]
KK_11_10=KK[12, 11]
KK_11_11=KK[12, 12]
KK_11_12=KK[12, 13]
KK_11_13=KK[12, 14]
KK_11_14=KK[12, 15]
KK_11_15=KK[12, 16]
KK_11_16=KK[12, 17]
KK_12_0=KK[13, 1]
KK_12_1=KK[13, 2]
KK_12_2=KK[13, 3]
KK_12_3=KK[13, 4]
KK_12_4=KK[13, 5]
KK_12_5=KK[13, 6]
KK_12_6=KK[13, 7]
KK_12_7=KK[13, 8]
KK_12_8=KK[13, 9]
KK_12_9=KK[13, 10]
KK_12_10=KK[13, 11]
KK_12_11=KK[13, 12]
KK_12_12=KK[13, 13]
KK_12_13=KK[13, 14]
KK_12_14=KK[13, 15]
KK_12_15=KK[13, 16]
KK_12_16=KK[13, 17]
KK_13_0=KK[14, 1]
KK_13_1=KK[14, 2]
KK_13_2=KK[14, 3]
KK_13_3=KK[14, 4]
KK_13_4=KK[14, 5]
KK_13_5=KK[14, 6]
KK_13_6=KK[14, 7]
KK_13_7=KK[14, 8]
KK_13_8=KK[14, 9]
KK_13_9=KK[14, 10]
KK_13_10=KK[14, 11]
KK_13_11=KK[14, 12]
KK_13_12=KK[14, 13]
KK_13_13=KK[14, 14]
KK_13_14=KK[14, 15]
KK_13_15=KK[14, 16]
KK_13_16=KK[14, 17]
KK_14_0=KK[15, 1]
KK_14_1=KK[15, 2]
KK_14_2=KK[15, 3]
KK_14_3=KK[15, 4]
KK_14_4=KK[15, 5]
KK_14_5=KK[15, 6]
KK_14_6=KK[15, 7]
KK_14_7=KK[15, 8]
KK_14_8=KK[15, 9]
KK_14_9=KK[15, 10]
KK_14_10=KK[15, 11]
KK_14_11=KK[15, 12]
KK_14_12=KK[15, 13]
KK_14_13=KK[15, 14]
KK_14_14=KK[15, 15]
KK_14_15=KK[15, 16]
KK_14_16=KK[15, 17]
KK_15_0=KK[16, 1]
KK_15_1=KK[16, 2]
KK_15_2=KK[16, 3]
KK_15_3=KK[16, 4]
KK_15_4=KK[16, 5]
KK_15_5=KK[16, 6]
KK_15_6=KK[16, 7]
KK_15_7=KK[16, 8]
KK_15_8=KK[16, 9]
KK_15_9=KK[16, 10]
KK_15_10=KK[16, 11]
KK_15_11=KK[16, 12]
KK_15_12=KK[16, 13]
KK_15_13=KK[16, 14]
KK_15_14=KK[16, 15]
KK_15_15=KK[16, 16]
KK_15_16=KK[16, 17]
KK_16_0=KK[17, 1]
KK_16_1=KK[17, 2]
KK_16_2=KK[17, 3]
KK_16_3=KK[17, 4]
KK_16_4=KK[17, 5]
KK_16_5=KK[17, 6]
KK_16_6=KK[17, 7]
KK_16_7=KK[17, 8]
KK_16_8=KK[17, 9]
KK_16_9=KK[17, 10]
KK_16_10=KK[17, 11]
KK_16_11=KK[17, 12]
KK_16_12=KK[17, 13]
KK_16_13=KK[17, 14]
KK_16_14=KK[17, 15]
KK_16_15=KK[17, 16]
KK_16_16=KK[17, 17]
times <- seq( from = 0, to = 366, by = 0.1)
S0 <- NN_0-2
E0 <- 0
I0 <- 2
A0 <- 0
S1 <- NN_1-2
E1 <- 0
I1 <- 2
A1 <- 0
S2 <- NN_2-2
E2 <- 0
I2 <- 2
A2 <- 0
S3 <- NN_3-2
E3 <- 0
I3 <- 2
A3 <- 0
S4 <- NN_4-2
E4 <- 0
I4 <- 2
A4 <- 0
S5 <- NN_5-2
E5 <- 0
I5 <- 2
A5 <- 0
S6 <- NN_6-2
E6 <- 0
I6 <- 2
A6 <- 0
S7 <- NN_7-2
E7 <- 0
I7 <- 2
A7 <- 0
S8 <- NN_8-2
E8 <- 0
I8 <- 2
A8 <- 0
S9 <- NN_9-2
E9 <- 0
I9 <- 2
A9 <- 0
S10 <- NN_10-2
E10 <- 0
I10 <- 2
A10 <- 0
S11 <- NN_11-2
E11 <- 0
I11 <- 2
A11 <- 0
S12 <- NN_12-2
E12 <- 0
I12 <- 2
A12 <- 0
S13 <- NN_13-2
E13 <- 0
I13 <- 2
A13 <- 0
S14 <- NN_14-2
E14 <- 0
I14 <- 2
A14 <- 0
S15 <- NN_15-2
E15 <- 0
I15 <- 2
A15 <- 0
S16 <- NN_16-2
E16 <- 0
I16 <- 2
A16 <- 0
params <- c( NN_0=NN_0, NN_1=NN_1, NN_2=NN_2, NN_3=NN_3, NN_4=NN_4, NN_5=NN_5, NN_6=NN_6, NN_7=NN_7, NN_8=NN_8, NN_9=NN_9, NN_10=NN_10, NN_11=NN_11, NN_12=NN_12, NN_13=NN_13, NN_14=NN_14, NN_15=NN_15, NN_16=NN_16, beta=beta, sigma=sigma, gamma=gamma, eigen=eigen, KK_0_0=KK[1, 1], KK_0_1=KK[1, 2], KK_0_2=KK[1, 3], KK_0_3=KK[1, 4], KK_0_4=KK[1, 5], KK_0_5=KK[1, 6], KK_0_6=KK[1, 7], KK_0_7=KK[1, 8], KK_0_8=KK[1, 9], KK_0_9=KK[1, 10], KK_0_10=KK[1, 11], KK_0_11=KK[1, 12], KK_0_12=KK[1, 13], KK_0_13=KK[1, 14], KK_0_14=KK[1, 15], KK_0_15=KK[1, 16], KK_0_16=KK[1, 17], KK_1_0=KK[2, 1], KK_1_1=KK[2, 2], KK_1_2=KK[2, 3], KK_1_3=KK[2, 4], KK_1_4=KK[2, 5], KK_1_5=KK[2, 6], KK_1_6=KK[2, 7], KK_1_7=KK[2, 8], KK_1_8=KK[2, 9], KK_1_9=KK[2, 10], KK_1_10=KK[2, 11], KK_1_11=KK[2, 12], KK_1_12=KK[2, 13], KK_1_13=KK[2, 14], KK_1_14=KK[2, 15], KK_1_15=KK[2, 16], KK_1_16=KK[2, 17], KK_2_0=KK[3, 1], KK_2_1=KK[3, 2], KK_2_2=KK[3, 3], KK_2_3=KK[3, 4], KK_2_4=KK[3, 5], KK_2_5=KK[3, 6], KK_2_6=KK[3, 7], KK_2_7=KK[3, 8], KK_2_8=KK[3, 9], KK_2_9=KK[3, 10], KK_2_10=KK[3, 11], KK_2_11=KK[3, 12], KK_2_12=KK[3, 13], KK_2_13=KK[3, 14], KK_2_14=KK[3, 15], KK_2_15=KK[3, 16], KK_2_16=KK[3, 17], KK_3_0=KK[4, 1], KK_3_1=KK[4, 2], KK_3_2=KK[4, 3], KK_3_3=KK[4, 4], KK_3_4=KK[4, 5], KK_3_5=KK[4, 6], KK_3_6=KK[4, 7], KK_3_7=KK[4, 8], KK_3_8=KK[4, 9], KK_3_9=KK[4, 10], KK_3_10=KK[4, 11], KK_3_11=KK[4, 12], KK_3_12=KK[4, 13], KK_3_13=KK[4, 14], KK_3_14=KK[4, 15], KK_3_15=KK[4, 16], KK_3_16=KK[4, 17], KK_4_0=KK[5, 1], KK_4_1=KK[5, 2], KK_4_2=KK[5, 3], KK_4_3=KK[5, 4], KK_4_4=KK[5, 5], KK_4_5=KK[5, 6], KK_4_6=KK[5, 7], KK_4_7=KK[5, 8], KK_4_8=KK[5, 9], KK_4_9=KK[5, 10], KK_4_10=KK[5, 11], KK_4_11=KK[5, 12], KK_4_12=KK[5, 13], KK_4_13=KK[5, 14], KK_4_14=KK[5, 15], KK_4_15=KK[5, 16], KK_4_16=KK[5, 17], KK_5_0=KK[6, 1], KK_5_1=KK[6, 2], KK_5_2=KK[6, 3], KK_5_3=KK[6, 4], KK_5_4=KK[6, 5], KK_5_5=KK[6, 6], KK_5_6=KK[6, 7], KK_5_7=KK[6, 8], KK_5_8=KK[6, 9], KK_5_9=KK[6, 10], KK_5_10=KK[6, 11], KK_5_11=KK[6, 12], KK_5_12=KK[6, 13], KK_5_13=KK[6, 14], KK_5_14=KK[6, 15], KK_5_15=KK[6, 16], KK_5_16=KK[6, 17], KK_6_0=KK[7, 1], KK_6_1=KK[7, 2], KK_6_2=KK[7, 3], KK_6_3=KK[7, 4], KK_6_4=KK[7, 5], KK_6_5=KK[7, 6], KK_6_6=KK[7, 7], KK_6_7=KK[7, 8], KK_6_8=KK[7, 9], KK_6_9=KK[7, 10], KK_6_10=KK[7, 11], KK_6_11=KK[7, 12], KK_6_12=KK[7, 13], KK_6_13=KK[7, 14], KK_6_14=KK[7, 15], KK_6_15=KK[7, 16], KK_6_16=KK[7, 17], KK_7_0=KK[8, 1], KK_7_1=KK[8, 2], KK_7_2=KK[8, 3], KK_7_3=KK[8, 4], KK_7_4=KK[8, 5], KK_7_5=KK[8, 6], KK_7_6=KK[8, 7], KK_7_7=KK[8, 8], KK_7_8=KK[8, 9], KK_7_9=KK[8, 10], KK_7_10=KK[8, 11], KK_7_11=KK[8, 12], KK_7_12=KK[8, 13], KK_7_13=KK[8, 14], KK_7_14=KK[8, 15], KK_7_15=KK[8, 16], KK_7_16=KK[8, 17], KK_8_0=KK[9, 1], KK_8_1=KK[9, 2], KK_8_2=KK[9, 3], KK_8_3=KK[9, 4], KK_8_4=KK[9, 5], KK_8_5=KK[9, 6], KK_8_6=KK[9, 7], KK_8_7=KK[9, 8], KK_8_8=KK[9, 9], KK_8_9=KK[9, 10], KK_8_10=KK[9, 11], KK_8_11=KK[9, 12], KK_8_12=KK[9, 13], KK_8_13=KK[9, 14], KK_8_14=KK[9, 15], KK_8_15=KK[9, 16], KK_8_16=KK[9, 17], KK_9_0=KK[10, 1], KK_9_1=KK[10, 2], KK_9_2=KK[10, 3], KK_9_3=KK[10, 4], KK_9_4=KK[10, 5], KK_9_5=KK[10, 6], KK_9_6=KK[10, 7], KK_9_7=KK[10, 8], KK_9_8=KK[10, 9], KK_9_9=KK[10, 10], KK_9_10=KK[10, 11], KK_9_11=KK[10, 12], KK_9_12=KK[10, 13], KK_9_13=KK[10, 14], KK_9_14=KK[10, 15], KK_9_15=KK[10, 16], KK_9_16=KK[10, 17], KK_10_0=KK[11, 1], KK_10_1=KK[11, 2], KK_10_2=KK[11, 3], KK_10_3=KK[11, 4], KK_10_4=KK[11, 5], KK_10_5=KK[11, 6], KK_10_6=KK[11, 7], KK_10_7=KK[11, 8], KK_10_8=KK[11, 9], KK_10_9=KK[11, 10], KK_10_10=KK[11, 11], KK_10_11=KK[11, 12], KK_10_12=KK[11, 13], KK_10_13=KK[11, 14], KK_10_14=KK[11, 15], KK_10_15=KK[11, 16], KK_10_16=KK[11, 17], KK_11_0=KK[12, 1], KK_11_1=KK[12, 2], KK_11_2=KK[12, 3], KK_11_3=KK[12, 4], KK_11_4=KK[12, 5], KK_11_5=KK[12, 6], KK_11_6=KK[12, 7], KK_11_7=KK[12, 8], KK_11_8=KK[12, 9], KK_11_9=KK[12, 10], KK_11_10=KK[12, 11], KK_11_11=KK[12, 12], KK_11_12=KK[12, 13], KK_11_13=KK[12, 14], KK_11_14=KK[12, 15], KK_11_15=KK[12, 16], KK_11_16=KK[12, 17], KK_12_0=KK[13, 1], KK_12_1=KK[13, 2], KK_12_2=KK[13, 3], KK_12_3=KK[13, 4], KK_12_4=KK[13, 5], KK_12_5=KK[13, 6], KK_12_6=KK[13, 7], KK_12_7=KK[13, 8], KK_12_8=KK[13, 9], KK_12_9=KK[13, 10], KK_12_10=KK[13, 11], KK_12_11=KK[13, 12], KK_12_12=KK[13, 13], KK_12_13=KK[13, 14], KK_12_14=KK[13, 15], KK_12_15=KK[13, 16], KK_12_16=KK[13, 17], KK_13_0=KK[14, 1], KK_13_1=KK[14, 2], KK_13_2=KK[14, 3], KK_13_3=KK[14, 4], KK_13_4=KK[14, 5], KK_13_5=KK[14, 6], KK_13_6=KK[14, 7], KK_13_7=KK[14, 8], KK_13_8=KK[14, 9], KK_13_9=KK[14, 10], KK_13_10=KK[14, 11], KK_13_11=KK[14, 12], KK_13_12=KK[14, 13], KK_13_13=KK[14, 14], KK_13_14=KK[14, 15], KK_13_15=KK[14, 16], KK_13_16=KK[14, 17], KK_14_0=KK[15, 1], KK_14_1=KK[15, 2], KK_14_2=KK[15, 3], KK_14_3=KK[15, 4], KK_14_4=KK[15, 5], KK_14_5=KK[15, 6], KK_14_6=KK[15, 7], KK_14_7=KK[15, 8], KK_14_8=KK[15, 9], KK_14_9=KK[15, 10], KK_14_10=KK[15, 11], KK_14_11=KK[15, 12], KK_14_12=KK[15, 13], KK_14_13=KK[15, 14], KK_14_14=KK[15, 15], KK_14_15=KK[15, 16], KK_14_16=KK[15, 17], KK_15_0=KK[16, 1], KK_15_1=KK[16, 2], KK_15_2=KK[16, 3], KK_15_3=KK[16, 4], KK_15_4=KK[16, 5], KK_15_5=KK[16, 6], KK_15_6=KK[16, 7], KK_15_7=KK[16, 8], KK_15_8=KK[16, 9], KK_15_9=KK[16, 10], KK_15_10=KK[16, 11], KK_15_11=KK[16, 12], KK_15_12=KK[16, 13], KK_15_13=KK[16, 14], KK_15_14=KK[16, 15], KK_15_15=KK[16, 16], KK_15_16=KK[16, 17], KK_16_0=KK[17, 1], KK_16_1=KK[17, 2], KK_16_2=KK[17, 3], KK_16_3=KK[17, 4], KK_16_4=KK[17, 5], KK_16_5=KK[17, 6], KK_16_6=KK[17, 7], KK_16_7=KK[17, 8], KK_16_8=KK[17, 9], KK_16_9=KK[17, 10], KK_16_10=KK[17, 11], KK_16_11=KK[17, 12], KK_16_12=KK[17, 13], KK_16_13=KK[17, 14], KK_16_14=KK[17, 15], KK_16_15=KK[17, 16], KK_16_16=KK[17, 17], zz_0 = zz_0, zz_1 = zz_1, zz_2 = zz_2, zz_3 = zz_3, zz_4 = zz_4, zz_5 = zz_5, zz_6 = zz_6, zz_7 = zz_7, zz_8 = zz_8, zz_9 = zz_9, zz_10 = zz_10, zz_11 = zz_11, zz_12 = zz_12, zz_13 = zz_13, zz_14 = zz_14, zz_15 = zz_15, zz_16 = zz_16)
res <- ode( c( ss_0=S0, ee_0=E0, ii_0=I0, aa_0=A0, ss_1=S1, ee_1=E1, ii_1=I1, aa_1=A1, ss_2=S2, ee_2=E2, ii_2=I2, aa_2=A2, ss_3=S3, ee_3=E3, ii_3=I3, aa_3=A3, ss_4=S4, ee_4=E4, ii_4=I4, aa_4=A4, ss_5=S5, ee_5=E5, ii_5=I5, aa_5=A5, ss_6=S6, ee_6=E6, ii_6=I6, aa_6=A6, ss_7=S7, ee_7=E7, ii_7=I7, aa_7=A7, ss_8=S8, ee_8=E8, ii_8=I8, aa_8=A8, ss_9=S9, ee_9=E9, ii_9=I9, aa_9=A9, ss_10=S10, ee_10=E10, ii_10=I10, aa_10=A10, ss_11=S11, ee_11=E11, ii_11=I11, aa_11=A11, ss_12=S12, ee_12=E12, ii_12=I12, aa_12=A12, ss_13=S13, ee_13=E13, ii_13=I13, aa_13=A13, ss_14=S14, ee_14=E14, ii_14=I14, aa_14=A14, ss_15=S15, ee_15=E15, ii_15=I15, aa_15=A15, ss_16=S16, ee_16=E16, ii_16=I16, aa_16=A16), times, model, params, method=rk4 )

tt <- res[,1]
cases <- res[,4]+res[,5]+res[,8]+res[,9]+res[,12]+res[,13]+res[,16]+res[,17]+res[,20]+res[,21]+res[,24]+res[,25]+res[,28]+res[,29]+res[,32]+res[,33]+res[,36]+res[,37]+res[,40]+res[,41]+res[,44]+res[,45]+res[,48]+res[,49]+res[,52]+res[,53]+res[,56]+res[,57]+res[,60]+res[,61]+res[,64]+res[,65]+res[,68]+res[,69]
cases <- cases[tt==floor(tt)]*gamma
tt <- tt[tt==floor(tt)]+1
cumul <- cumsum(cases)
dat <- read.table("DeryaSE-summary-Daily.dat", header=FALSE)
dat$Cases <- dat$V2 + dat$V3 + dat$V4 + dat$V5
dat$Cumul <- cumsum(dat$Cases)

plot(  tt, cases, type='l', col='red', lwd=2, ylim=c(0,NN), xlab='Time', ylab='Incidence in population' )
lines( tt, cumul, col='red', lt=2 )
lines( dat$V1, dat$Cases, col='blue', lt=1 )
lines( dat$V1, dat$Cumul, col='blue', lt=2 )
legend( 'topright', legend=c('I'), bg=rgb(1,1,1), lwd=2, col=c(3,2) )
