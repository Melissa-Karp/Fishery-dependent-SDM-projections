### Fishery-dependent SDM (L^3) Estimation Model Code
## Function to run BRTs on data from OM, with all Covariates (temp, chla, mld)
## Function returns only the fitted and predicted values
## modified my M.Karp 1/15/21

#############################
#   Full Models- with Chla and space-time component #
#############################
#
###Random sampling
#
if("ran" %in% sampling) {
  print("Fitting BRT-Ran_te")
  brt_R_P_te <- gbm.step(data=dat_hist_random,
                        gbm.x = c(1:3, 22, 24, 25),
                        gbm.y = 'pres',
                        family = "bernoulli",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  brt_R_N_te <- gbm.step(data=dat_hist_random[dat_hist_random$abundance>0,], 
                        gbm.x = c(1:3, 22, 24, 25),
                        gbm.y = 'log_abundance',
                        family = "gaussian",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  presRx_te <- predict(brt_R_P_te, dat_hist, n.trees=brt_R_P_te$gbm.call$best.trees, type="response")
  abundRx_te <- exp(predict(brt_R_N_te, dat_hist, n.trees=brt_R_N_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_Ran_te <- presRx_te * abundRx_te
  
  presRx_te <- predict(brt_R_P_te, dat_fcast, n.trees=brt_R_P_te$gbm.call$best.trees, type="response")
  abundRx_te <- exp(predict(brt_R_N_te, dat_fcast, n.trees=brt_R_N_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Ran_te <- presRx_te * abundRx_te
}

#
### Preferential (target) sampling
#
#Target sampling-0.5
if("tar_0.5" %in% sampling){
  print("Fitting BRT-Pref0.5_te")
  brt_T_P_1_te <- gbm.step(data=dat_hist_Tar_1,
                          gbm.x = c(1:3, 22, 24, 25),
                          gbm.y = 'pres',
                          family = "bernoulli",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  brt_T_N_1_te <- gbm.step(data=dat_hist_Tar_1[dat_hist_Tar_1$abundance>0,],
                          gbm.x = c(1:3, 22, 24, 25),
                          gbm.y = 'log_abundance',
                          family = "gaussian",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  presTx_1_te <- predict(brt_T_P_1_te, dat_hist, n.trees=brt_T_P_1_te$gbm.call$best.trees, type="response")
  abundTx_1_te <- exp(predict(brt_T_N_1_te, dat_hist, n.trees=brt_T_N_1_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_Tar_0.5_te <- presTx_1_te * abundTx_1_te
  
  presTx_1_te <- predict(brt_T_P_1_te, dat_fcast, n.trees=brt_T_P_1_te$gbm.call$best.trees, type="response")
  abundTx_1_te <- exp(predict(brt_T_N_1_te, dat_fcast, n.trees=brt_T_N_1_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Tar_0.5_te <- presTx_1_te * abundTx_1_te
}

#Target sampling-0.6
if("tar_0.6" %in% sampling){
  print("Fitting BRT-Pref0.6_te")
  brt_T_P_2_te <- gbm.step(data=dat_hist_Tar_2,
                          gbm.x = c(1:3, 22, 24, 25),
                          gbm.y = 'pres',
                          family = "bernoulli",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  brt_T_N_2_te <- gbm.step(data=dat_hist_Tar_2[dat_hist_Tar_2$abundance>0,],
                          gbm.x = c(1:3, 22, 24, 25),
                          gbm.y = 'log_abundance',
                          family = "gaussian",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  presTx_2_te <- predict(brt_T_P_2_te, dat_hist, n.trees=brt_T_P_2_te$gbm.call$best.trees, type="response")
  abundTx_2_te <- exp(predict(brt_T_N_2_te, dat_hist, n.trees=brt_T_N_2_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_Tar_0.6_te <- presTx_2_te* abundTx_2_te
  
  presTx_2_te <- predict(brt_T_P_2_te, dat_fcast, n.trees=brt_T_P_2_te$gbm.call$best.trees, type="response")
  abundTx_2_te <- exp(predict(brt_T_N_2_te, dat_fcast, n.trees=brt_T_N_2_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Tar_0.6_te <- presTx_2_te * abundTx_2_te
}

#Target sampling-0.7
if("tar_0.7" %in% sampling){
  print("Fitting BRT-Pref0.7_te")
  brt_T_P_3_te <- gbm.step(data=dat_hist_Tar_3,
                          gbm.x = c(1:3, 22, 24, 25),
                          gbm.y = 'pres',
                          family = "bernoulli",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  brt_T_N_3_te <- gbm.step(data=dat_hist_Tar_3[dat_hist_Tar_3$abundance>0,],
                          gbm.x = c(1:3, 22, 24, 25),
                          gbm.y = 'log_abundance',
                          family = "gaussian",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  presTx_3_te <- predict(brt_T_P_3_te, dat_hist, n.trees=brt_T_P_3_te$gbm.call$best.trees, type="response")
  abundTx_3_te <- exp(predict(brt_T_N_3_te, dat_hist, n.trees=brt_T_N_3_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_Tar_0.7_te <- presTx_3_te * abundTx_3_te
  
  presTx_3_te <- predict(brt_T_P_3_te, dat_fcast, n.trees=brt_T_P_3_te$gbm.call$best.trees, type="response")
  abundTx_3_te <- exp(predict(brt_T_N_3_te, dat_fcast, n.trees=brt_T_N_3_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Tar_0.7_te <- presTx_3_te * abundTx_3_te
}

#Target sampling-0.8
if("tar_0.8" %in% sampling){
  print("Fitting BRT-Pref0.8_te")
  brt_T_P_4_te <- gbm.step(data=dat_hist_Tar_4,
                          gbm.x = c(1:3, 22, 24, 25),
                          gbm.y = 'pres',
                          family = "bernoulli",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  brt_T_N_4_te <- gbm.step(data=dat_hist_Tar_4[dat_hist_Tar_4$abundance>0,],
                          gbm.x = c(1:3, 22, 24, 25),
                          gbm.y = 'log_abundance',
                          family = "gaussian",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  presTx_4_te <- predict(brt_T_P_4_te, dat_hist, n.trees=brt_T_P_4_te$gbm.call$best.trees, type="response")
  abundTx_4_te <- exp(predict(brt_T_N_4_te, dat_hist, n.trees=brt_T_N_4_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_Tar_0.8_te <- presTx_4_te * abundTx_4_te
  
  presTx_4_te <- predict(brt_T_P_4_te, dat_fcast, n.trees=brt_T_P_4_te$gbm.call$best.trees, type="response")
  abundTx_4_te <- exp(predict(brt_T_N_4_te, dat_fcast, n.trees=brt_T_N_4_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Tar_0.8_te <- presTx_4_te * abundTx_4_te
}

#Target sampling-0.9
if("tar_0.9" %in% sampling){
  print("Fitting BRT-Pref0.9_te")
  brt_T_P_5_te <- gbm.step(data=dat_hist_Tar_5,
                          gbm.x = c(1:3, 22, 24, 25),
                          gbm.y = 'pres',
                          family = "bernoulli",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  brt_T_N_5_te <- gbm.step(data=dat_hist_Tar_5[dat_hist_Tar_5$abundance>0,],
                          gbm.x = c(1:3, 22, 24, 25),
                          gbm.y = 'log_abundance',
                          family = "gaussian",
                          tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                          plot.main=FALSE, verbose = FALSE)
  
  presTx_5_te <- predict(brt_T_P_5_te, dat_hist, n.trees=brt_T_P_5_te$gbm.call$best.trees, type="response")
  abundTx_5_te <- exp(predict(brt_T_N_5_te, dat_hist, n.trees=brt_T_N_5_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_Tar_0.9_te <- presTx_5_te * abundTx_5_te
  
  presTx_5_te <- predict(brt_T_P_5_te, dat_fcast, n.trees=brt_T_P_5_te$gbm.call$best.trees, type="response")
  abundTx_5_te <- exp(predict(brt_T_N_5_te, dat_fcast, n.trees=brt_T_N_5_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Tar_0.9_te <- presTx_5_te * abundTx_5_te
}

#
### Distance from Port Sampling
#
#Northern Offshore
if("npo" %in% sampling){
  print("Fitting BRT-NPO_te")
  brt_dist_P_npo_te <- gbm.step(data=dat_hist_Dist_npo,
                               gbm.x = c(1:3, 22, 24, 25),
                               gbm.y = 'pres',
                               family = "bernoulli",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_npo_te <- gbm.step(data=dat_hist_Dist_npo[dat_hist_Dist_npo$abundance>0,],
                               gbm.x = c(1:3, 22, 24, 25),
                               gbm.y = 'log_abundance',
                               family = "gaussian",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  presDx_npo_te <- predict(brt_dist_P_npo_te, dat_hist, n.trees=brt_dist_P_npo_te$gbm.call$best.trees, type="response")
  abundDx_npo_te <- exp(predict(brt_dist_N_npo_te, dat_hist, n.trees=brt_dist_N_npo_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_npo_te <- presDx_npo_te * abundDx_npo_te
  
  presDx_npo_te <- predict(brt_dist_P_npo_te, dat_fcast, n.trees=brt_dist_P_npo_te$gbm.call$best.trees, type="response")
  abundDx_npo_te <- exp(predict(brt_dist_N_npo_te, dat_fcast, n.trees=brt_dist_N_npo_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_npo_te <- presDx_npo_te * abundDx_npo_te
}

#Northern Nearshore
if("npn" %in% sampling){
  print("Fitting BRT-NPN_te")
  brt_dist_P_npn_te <- gbm.step(data=dat_hist_Dist_npn,
                               gbm.x = c(1:3, 22, 24, 25),
                               gbm.y = 'pres',
                               family = "bernoulli",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_npn_te <- gbm.step(data=dat_hist_Dist_npn[dat_hist_Dist_npn$abundance>0,],
                               gbm.x = c(1:3, 22, 24, 25),
                               gbm.y = 'log_abundance',
                               family = "gaussian",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  presDx_npn_te <- predict(brt_dist_P_npn_te, dat_hist, n.trees=brt_dist_P_npn_te$gbm.call$best.trees, type="response")
  abundDx_npn_te <- exp(predict(brt_dist_N_npn_te, dat_hist, n.trees=brt_dist_N_npn_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_npn_te <- presDx_npn_te * abundDx_npn_te
  
  presDx_npn_te <- predict(brt_dist_P_npn_te, dat_fcast, n.trees=brt_dist_P_npn_te$gbm.call$best.trees, type="response")
  abundDx_npn_te <- exp(predict(brt_dist_N_npn_te, dat_fcast, n.trees=brt_dist_N_npn_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_npn_te <- presDx_npn_te * abundDx_npn_te
}  

#Middle Offshore
if("mpo" %in% sampling){
  print("Fitting BRT-MPO_te")
  brt_dist_P_mpo_te <- gbm.step(data=dat_hist_Dist_mpo,
                               gbm.x = c(1:3, 22, 24, 25),
                               gbm.y = 'pres',
                               family = "bernoulli",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_mpo_te <- gbm.step(data=dat_hist_Dist_mpo[dat_hist_Dist_mpo$abundance>0,],
                               gbm.x = c(1:3, 22, 24, 25),
                               gbm.y = 'log_abundance',
                               family = "gaussian",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  presDx_mpo_te <- predict(brt_dist_P_mpo_te, dat_hist, n.trees=brt_dist_P_mpo_te$gbm.call$best.trees, type="response")
  abundDx_mpo_te <- exp(predict(brt_dist_N_mpo_te, dat_hist, n.trees=brt_dist_N_mpo_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_mpo_te <- presDx_mpo_te * abundDx_mpo_te
  
  presDx_mpo_te <- predict(brt_dist_P_mpo_te, dat_fcast, n.trees=brt_dist_P_mpo_te$gbm.call$best.trees, type="response")
  abundDx_mpo_te <- exp(predict(brt_dist_N_mpo_te, dat_fcast, n.trees=brt_dist_N_mpo_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_mpo_te <- presDx_mpo_te * abundDx_mpo_te
}

#Middle nearshore
if("mpn" %in% sampling){
  print("Fitting BRT-MPN_te")
  brt_dist_P_mpn_te <- gbm.step(data=dat_hist_Dist_mpn,
                               gbm.x = c(1:3, 22, 24, 25),
                               gbm.y = 'pres',
                               family = "bernoulli",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_mpn_te <- gbm.step(data=dat_hist_Dist_mpn[dat_hist_Dist_mpn$abundance>0,],
                               gbm.x = c(1:3, 22, 24, 25),
                               gbm.y = 'log_abundance',
                               family = "gaussian",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  presDx_mpn_te <- predict(brt_dist_P_mpn_te, dat_hist, n.trees=brt_dist_P_mpn_te$gbm.call$best.trees, type="response")
  abundDx_mpn_te <- exp(predict(brt_dist_N_mpn_te, dat_hist, n.trees=brt_dist_N_mpn_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_mpn_te <- presDx_mpn_te * abundDx_mpn_te
  
  presDx_mpn_te <- predict(brt_dist_P_mpn_te, dat_fcast, n.trees=brt_dist_P_mpn_te$gbm.call$best.trees, type="response")
  abundDx_mpn_te <- exp(predict(brt_dist_N_mpn_te, dat_fcast, n.trees=brt_dist_N_mpn_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_mpn_te <- presDx_mpn_te * abundDx_mpn_te
}

#Southern Offshore
if("spo" %in% sampling){
  print("Fitting BRT-SPO_te")
  brt_dist_P_spo_te <- gbm.step(data=dat_hist_Dist_spo,
                               gbm.x = c(1:3, 22, 24, 25),
                               gbm.y = 'pres',
                               family = "bernoulli",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_spo_te <- gbm.step(data=dat_hist_Dist_spo[dat_hist_Dist_spo$abundance>0,],
                               gbm.x = c(1:3, 22, 24, 25),
                               gbm.y = 'log_abundance',
                               family = "gaussian",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  presDx_spo_te <- predict(brt_dist_P_spo_te, dat_hist, n.trees=brt_dist_P_spo_te$gbm.call$best.trees, type="response")
  abundDx_spo_te <- exp(predict(brt_dist_N_spo_te, dat_hist, n.trees=brt_dist_N_spo_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_spo_te <- presDx_spo_te * abundDx_spo_te
  
  presDx_spo_te <- predict(brt_dist_P_spo_te, dat_fcast, n.trees=brt_dist_P_spo_te$gbm.call$best.trees, type="response")
  abundDx_spo_te <- exp(predict(brt_dist_N_spo_te, dat_fcast, n.trees=brt_dist_N_spo_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_spo_te <- presDx_spo_te * abundDx_spo_te
}

#Southern Nearshore
if("spn" %in% sampling){
  print("Fitting BRT-SPN_te")
  brt_dist_P_spn_te <- gbm.step(data=dat_hist_Dist_spn,
                               gbm.x = c(1:3, 22, 24, 25),
                               gbm.y = 'pres',
                               family = "bernoulli",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_spn_te <- gbm.step(data=dat_hist_Dist_spn[dat_hist_Dist_spn$abundance>0,],
                               gbm.x = c(1:3, 22, 24, 25),
                               gbm.y = 'log_abundance',
                               family = "gaussian",
                               tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                               plot.main=FALSE, verbose = FALSE)
  
  presDx_spn_te <- predict(brt_dist_P_spn_te, dat_hist, n.trees=brt_dist_P_spn_te$gbm.call$best.trees, type="response")
  abundDx_spn_te <- exp(predict(brt_dist_N_spn_te, dat_hist, n.trees=brt_dist_N_spn_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_spn_te <- presDx_spn_te * abundDx_spn_te
  
  presDx_spn_te <- predict(brt_dist_P_spn_te, dat_fcast, n.trees=brt_dist_P_spn_te$gbm.call$best.trees, type="response")
  abundDx_spn_te <- exp(predict(brt_dist_N_spn_te, dat_fcast, n.trees=brt_dist_N_spn_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_spn_te <- presDx_spn_te * abundDx_spn_te
}

#All Offshore
if("allo" %in% sampling){
  print("Fitting BRT-ALLO_te")
  brt_dist_P_allo_te <- gbm.step(data=dat_hist_Dist_allo,
                                gbm.x = c(1:3, 22, 24, 25),
                                gbm.y = 'pres',
                                family = "bernoulli",
                                tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_allo_te <- gbm.step(data=dat_hist_Dist_allo[dat_hist_Dist_allo$abundance>0,],
                                gbm.x = c(1:3, 22, 24, 25),
                                gbm.y = 'log_abundance',
                                family = "gaussian",
                                tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                plot.main=FALSE, verbose = FALSE)
  
  presDx_allo_te <- predict(brt_dist_P_allo_te, dat_hist, n.trees=brt_dist_P_allo_te$gbm.call$best.trees, type="response")
  abundDx_allo_te <- exp(predict(brt_dist_N_allo_te, dat_hist, n.trees=brt_dist_N_allo_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_allo_te <- presDx_allo_te * abundDx_allo_te
  
  presDx_allo_te <- predict(brt_dist_P_allo_te, dat_fcast, n.trees=brt_dist_P_allo_te$gbm.call$best.trees, type="response")
  abundDx_allo_te <- exp(predict(brt_dist_N_allo_te, dat_fcast, n.trees=brt_dist_N_allo_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_allo_te <- presDx_allo_te * abundDx_allo_te
}

#All Nearshore
if("alln" %in% sampling){
  print("Fitting BRT-ALLN_te")
  brt_dist_P_alln_te <- gbm.step(data=dat_hist_Dist_alln,
                                gbm.x = c(1:3, 22, 24, 25),
                                gbm.y = 'pres',
                                family = "bernoulli",
                                tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                plot.main=FALSE, verbose = FALSE)
  
  brt_dist_N_alln_te <- gbm.step(data=dat_hist_Dist_alln[dat_hist_Dist_alln$abundance>0,],
                                gbm.x = c(1:3, 22, 24, 25),
                                gbm.y = 'log_abundance',
                                family = "gaussian",
                                tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                                plot.main=FALSE, verbose = FALSE)
  
  presDx_alln_te <- predict(brt_dist_P_alln_te, dat_hist, n.trees=brt_dist_P_alln_te$gbm.call$best.trees, type="response")
  abundDx_alln_te <- exp(predict(brt_dist_N_alln_te, dat_hist, n.trees=brt_dist_N_alln_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_Dist_alln_te <- presDx_alln_te * abundDx_alln_te
  
  presDx_alln_te <- predict(brt_dist_P_alln_te, dat_fcast, n.trees=brt_dist_P_alln_te$gbm.call$best.trees, type="response")
  abundDx_alln_te <- exp(predict(brt_dist_N_alln_te, dat_fcast, n.trees=brt_dist_N_alln_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_Dist_alln_te <- presDx_alln_te * abundDx_alln_te
}
#
## Bycatch + Opt Target Sampling
#
if("BY_2" %in% sampling){
  print("Fitting BRT-BY_te")
  brt_B_P_te <- gbm.step(data=dat_hist_BY_2,
                        gbm.x = c(1:3, 22, 24, 25),
                        gbm.y = 'pres',
                        family = "bernoulli",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  brt_B_N_te <- gbm.step(data=dat_hist_BY_2[dat_hist_BY_2$abundance>0,], 
                        gbm.x = c(1:3, 22, 24, 25),
                        gbm.y = 'log_abundance',
                        family = "gaussian",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
  
  presBx_te <- predict(brt_B_P_te, dat_hist, n.trees=brt_B_P_te$gbm.call$best.trees, type="response")
  abundBx_te <- exp(predict(brt_B_N_te, dat_hist, n.trees=brt_B_N_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_BY_te <- presBx_te * abundBx_te
  
  presBx_te <- predict(brt_B_P_te, dat_fcast, n.trees=brt_B_P_te$gbm.call$best.trees, type="response")
  abundBx_te <- exp(predict(brt_B_N_te, dat_fcast, n.trees=brt_B_N_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_BY_te<- presBx_te * abundBx_te
}

#
###Closed Areas + Opt Target Species
#
# Small Closed Area
if("CA_sm" %in% sampling){
  print("Fitting BRT-CASM_te")
  brt_CA_P_sm_te <- gbm.step(data=dat_hist_CA_sm,
                            gbm.x = c(1:3, 22, 24, 25),
                            gbm.y = 'pres',
                            family = "bernoulli",
                            tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                            plot.main=FALSE, verbose = FALSE)
  
  brt_CA_N_sm_te <- gbm.step(data=dat_hist_CA_sm[dat_hist_CA_sm$abundance>0,], 
                            gbm.x = c(1:3, 22, 24, 25),
                            gbm.y = 'log_abundance',
                            family = "gaussian",
                            tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                            plot.main=FALSE, verbose = FALSE)
  
  presCAx_sm_te <- predict(brt_CA_P_sm_te, dat_hist, n.trees=brt_CA_P_sm_te$gbm.call$best.trees, type="response")
  abundCAx_sm_te <- exp(predict(brt_CA_N_sm_te, dat_hist, n.trees=brt_CA_N_sm_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_CA_sm_te <- presCAx_sm_te * abundCAx_sm_te
  
  presCAx_sm_te <- predict(brt_CA_P_sm_te, dat_fcast, n.trees=brt_CA_P_sm_te$gbm.call$best.trees, type="response")
  abundCAx_sm_te <- exp(predict(brt_CA_N_sm_te, dat_fcast, n.trees=brt_CA_N_sm_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_CA_sm_te<- presCAx_sm_te * abundCAx_sm_te
}

# Medium Closed Area
if("CA_med" %in% sampling){
  print("Fitting BRT-CAMED_te")
  brt_CA_P_med_te <- gbm.step(data=dat_hist_CA_med,
                             gbm.x = c(1:3, 22, 24, 25),
                             gbm.y = 'pres',
                             family = "bernoulli",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  brt_CA_N_med_te <- gbm.step(data=dat_hist_CA_med[dat_hist_CA_med$abundance>0,], 
                             gbm.x = c(1:3, 22, 24, 25),
                             gbm.y = 'log_abundance',
                             family = "gaussian",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  presCAx_med_te <- predict(brt_CA_P_med_te, dat_hist, n.trees=brt_CA_P_med_te$gbm.call$best.trees, type="response")
  abundCAx_med_te <- exp(predict(brt_CA_N_med_te, dat_hist, n.trees=brt_CA_N_med_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_CA_med_te <- presCAx_med_te * abundCAx_med_te
  
  presCAx_med_te <- predict(brt_CA_P_med_te, dat_fcast, n.trees=brt_CA_P_med_te$gbm.call$best.trees, type="response")
  abundCAx_med_te <- exp(predict(brt_CA_N_med_te, dat_fcast, n.trees=brt_CA_N_med_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_CA_med_te<- presCAx_med_te * abundCAx_med_te
}

# Large Closed Area
if("CA_lar" %in% sampling){
  print("Fitting BRT-CALAR_te")
  brt_CA_P_lar_te <- gbm.step(data=dat_hist_CA_lar,
                             gbm.x = c(1:3, 22, 24, 25),
                             gbm.y = 'pres',
                             family = "bernoulli",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  brt_CA_N_lar_te <- gbm.step(data=dat_hist_CA_lar[dat_hist_CA_lar$abundance>0,], 
                             gbm.x = c(1:3, 22, 24, 25),
                             gbm.y = 'log_abundance',
                             family = "gaussian",
                             tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                             plot.main=FALSE, verbose = FALSE)
  
  presCAx_lar_te <- predict(brt_CA_P_lar_te, dat_hist, n.trees=brt_CA_P_lar_te$gbm.call$best.trees, type="response")
  abundCAx_lar_te <- exp(predict(brt_CA_N_lar_te, dat_hist, n.trees=brt_CA_N_lar_te$gbm.call$best.trees, type="response"))
  dat_hist$brt_CA_lar_te <- presCAx_lar_te * abundCAx_lar_te
  
  presCAx_lar_te <- predict(brt_CA_P_lar_te, dat_fcast, n.trees=brt_CA_P_lar_te$gbm.call$best.trees, type="response")
  abundCAx_lar_te <- exp(predict(brt_CA_N_lar_te, dat_fcast, n.trees=brt_CA_N_lar_te$gbm.call$best.trees, type="response"))
  dat_fcast$brt_CA_lar_te<- presCAx_lar_te * abundCAx_lar_te
}
