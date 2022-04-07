### Fishery-dependent SDM (L^3) Estimation Model Code
## Function to run GAMs on data from OM
## Function returns only the fitted and predicted values
#modified my M.Karp 1/04/20
library(mgcv)
library(pROC)
library(ROCR)

#############################
#   Full Models- with Chla  #
#############################

##### Random sampling ####
#
  if("ran" %in% sampling){
    print("Fitting GAM-Ran")
    gam_Ran_P <- gam(pres~s(temp) + s(mld) + s(chl_surface), data=dat_hist_random, family=binomial) #note - chl-surface is log of chl-surface
    #plot(gam_Ran_P, pages=1)
    gam_Ran_N <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface), data=dat_hist_random[dat_hist_random$abundance>0,], family=gaussian)
    #plot(gam_Ran_N, pages=1)
    
    dat_hist$presxR <- predict(gam_Ran_P, dat_hist, type="response")
    abundxR <- predict(gam_Ran_N, dat_hist, type="response") 
    dat_hist$gam_Ran <- dat_hist$presxR * exp(abundxR)  #predicted catch
    
    dat_fcast$presxR <- predict(gam_Ran_P, dat_fcast, type="response")
    abundxR <- predict(gam_Ran_N, dat_fcast, type="response")
    dat_fcast$gam_Ran <- dat_fcast$presxR * exp(abundxR)
  }

#
##### Preferential Sampling #######
#
  #Pref sampling - 0.5
  if("tar_0.5" %in% sampling) {
    print("Fitting GAM-Pref0.5")
    gam_Tar_P_1 <- gam(pres~s(temp)+ s(mld) + s(chl_surface), data=dat_hist_Tar_1, family=binomial)
    #plot(gam_Tar_P_1, pages=1)
    gam_Tar_N_1 <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface), data=dat_hist_Tar_1[dat_hist_Tar_1$abundance>0,], family=gaussian)
    #plot(gam_Tar_N_1, pages=1)
    
    dat_hist$presxT_1 <- predict(gam_Tar_P_1, dat_hist, type="response")
    abundxT_1 <- predict(gam_Tar_N_1, dat_hist, type="response")
    dat_hist$gam_Tar_0.5 <- dat_hist$presxT_1 * exp(abundxT_1)
    
    dat_fcast$presxT_1 <- predict(gam_Tar_P_1, dat_fcast, type="response")
    abundxT_1 <- predict(gam_Tar_N_1, dat_fcast, type="response")
    dat_fcast$gam_Tar_0.5 <- dat_fcast$presxT_1 * exp(abundxT_1)
  }
  #Pref sampling -0.6
  if("tar_0.6" %in% sampling) {
    print("Fitting GAM-Pref0.6")
    gam_Tar_P_2 <- gam(pres~s(temp)+ s(mld) + s(chl_surface), data=dat_hist_Tar_2, family=binomial)
    #plot(gam_Tar_P_2, pages=1)
    gam_Tar_N_2 <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface), data=dat_hist_Tar_2[dat_hist_Tar_2$abundance>0,], family=gaussian)
    #plot(gam_Tar_N_2, pages=1)
    
    dat_hist$presxT_2 <- predict(gam_Tar_P_2, dat_hist, type="response")
    abundxT_2 <- predict(gam_Tar_N_2, dat_hist, type="response")
    dat_hist$gam_Tar_0.6 <- dat_hist$presxT_2 * exp(abundxT_2)
    
    dat_fcast$presxT_2 <- predict(gam_Tar_P_2, dat_fcast, type="response")
    abundxT_2 <- predict(gam_Tar_N_2, dat_fcast, type="response")
    dat_fcast$gam_Tar_0.6 <- dat_fcast$presxT_2 * exp(abundxT_2)
  }
  #Pref sampling -0.7
  if("tar_0.7" %in% sampling) {
    print("Fitting GAM-Pref0.7")
    gam_Tar_P_3 <- gam(pres~s(temp)+ s(mld) + s(chl_surface), data=dat_hist_Tar_3, family=binomial)
    #plot(gam_Tar_P_3, pages=1)
    gam_Tar_N_3 <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface), data=dat_hist_Tar_3[dat_hist_Tar_3$abundance>0,], family=gaussian)
    #plot(gam_Tar_N_3, pages=1)
    
    dat_hist$presxT_3 <- predict(gam_Tar_P_3, dat_hist, type="response")
    abundxT_3 <- predict(gam_Tar_N_3, dat_hist, type="response")
    dat_hist$gam_Tar_0.7 <- dat_hist$presxT_3 * exp(abundxT_3)
    
    dat_fcast$presxT_3 <- predict(gam_Tar_P_3, dat_fcast, type="response")
    abundxT_3 <- predict(gam_Tar_N_3, dat_fcast, type="response")
    dat_fcast$gam_Tar_0.7 <- dat_fcast$presxT_3 * exp(abundxT_3)
  }
  #Pref sampling -0.8
  if("tar_0.8" %in% sampling) {
    print("Fitting GAM-Pref0.8")
    gam_Tar_P_4 <- gam(pres~s(temp)+ s(mld) + s(chl_surface), data=dat_hist_Tar_4, family=binomial)
    #plot(gam_Tar_P_4, pages=1)
    gam_Tar_N_4 <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface), data=dat_hist_Tar_4[dat_hist_Tar_4$abundance>0,], family=gaussian)
    #plot(gam_Tar_N_4, pages=1)
    
    dat_hist$presxT_4 <- predict(gam_Tar_P_4, dat_hist, type="response")
    abundxT_4 <- predict(gam_Tar_N_4, dat_hist, type="response")
    dat_hist$gam_Tar_0.8 <- dat_hist$presxT_4 * exp(abundxT_4)
    
    dat_fcast$presxT_4 <- predict(gam_Tar_P_4, dat_fcast, type="response")
    abundxT_4 <- predict(gam_Tar_N_4, dat_fcast, type="response")
    dat_fcast$gam_Tar_0.8 <- dat_fcast$presxT_4 * exp(abundxT_4)
  }
  #Pref sampling -0.9
  if("tar_0.9" %in% sampling) {
    print("Fitting GAM-Pref0.9")
    gam_Tar_P_5 <- gam(pres~s(temp)+ s(mld) + s(chl_surface), data=dat_hist_Tar_5, family=binomial)
    #plot(gam_Tar_P_5, pages=1)
    gam_Tar_N_5 <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface), data=dat_hist_Tar_5[dat_hist_Tar_5$abundance>0,], family=gaussian)
    #plot(gam_Tar_N_5, pages=1)
    
    dat_hist$presxT_5 <- predict(gam_Tar_P_5, dat_hist, type="response")
    abundxT_5 <- predict(gam_Tar_N_5, dat_hist, type="response")
    dat_hist$gam_Tar_0.9 <- dat_hist$presxT_5 * exp(abundxT_5)
    
    dat_fcast$presxT_5 <- predict(gam_Tar_P_5, dat_fcast, type="response")
    abundxT_5 <- predict(gam_Tar_N_5, dat_fcast, type="response")
    dat_fcast$gam_Tar_0.9 <- dat_fcast$presxT_5 * exp(abundxT_5)
  }

#
##### Distance from Port Sampling #######
#
  #Dist sampling - NPO
  if("npo" %in% sampling){
    print("Fitting GAM-NPO")
    gam_Dist_P_npo <- gam(pres~s(temp) + s(mld) + s(chl_surface), data=dat_hist_Dist_npo, family=binomial)
    #plot(gam_Dist_P_npo, pages=1)
    gam_Dist_N_npo <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface), data=dat_hist_Dist_npo[dat_hist_Dist_npo$abundance>0,], family=gaussian)
    #plot(gam_Dist_N_npo, pages=1)
    
    dat_hist$presxD_npo <- predict(gam_Dist_P_npo, dat_hist, type="response")
    abundxD_npo <- predict(gam_Dist_N_npo, dat_hist, type="response")
    dat_hist$gam_Dist_npo <- dat_hist$presxD_npo * exp(abundxD_npo)
    
    dat_fcast$presxD_npo <- predict(gam_Dist_P_npo, dat_fcast, type="response")
    abundxD_npo <- predict(gam_Dist_N_npo, dat_fcast, type="response")
    dat_fcast$gam_Dist_npo <- dat_fcast$presxD_npo * exp(abundxD_npo)
  }
  #Dist sampling - NPN
  if("npn" %in% sampling){
    print("Fitting GAM-NPN")
    gam_Dist_P_npn <- gam(pres~s(temp) + s(mld) + s(chl_surface), data=dat_hist_Dist_npn, family=binomial)
    #plot(gam_Dist_P_npn, pages=1)
    gam_Dist_N_npn <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface), data=dat_hist_Dist_npn[dat_hist_Dist_npn$abundance>0,], family=gaussian)
    #plot(gam_Dist_N_npn, pages=1)
    
    dat_hist$presxD_npn <- predict(gam_Dist_P_npn, dat_hist, type="response")
    abundxD_npn <- predict(gam_Dist_N_npn, dat_hist, type="response")
    dat_hist$gam_Dist_npn <- dat_hist$presxD_npn * exp(abundxD_npn)
    
    dat_fcast$presxD_npn <- predict(gam_Dist_P_npn, dat_fcast, type="response")
    abundxD_npn <- predict(gam_Dist_N_npn, dat_fcast, type="response")
    dat_fcast$gam_Dist_npn <- dat_fcast$presxD_npn * exp(abundxD_npn)
  }
  #Dist sampling - MPO
  if("mpo" %in% sampling){
    print("Fitting GAM-MPO")
    gam_Dist_P_mpo <- gam(pres~s(temp) + s(mld) + s(chl_surface), data=dat_hist_Dist_mpo, family=binomial)
    #plot(gam_Dist_P_mpo, pages=1)
    gam_Dist_N_mpo <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface), data=dat_hist_Dist_mpo[dat_hist_Dist_mpo$abundance>0,], family=gaussian)
    #plot(gam_Dist_N_mpo, pages=1)
    
    dat_hist$presxD_mpo <- predict(gam_Dist_P_mpo, dat_hist, type="response")
    abundxD_mpo <- predict(gam_Dist_N_mpo, dat_hist, type="response")
    dat_hist$gam_Dist_mpo <- dat_hist$presxD_mpo * exp(abundxD_mpo)
    
    dat_fcast$presxD_mpo <- predict(gam_Dist_P_mpo, dat_fcast, type="response")
    abundxD_mpo <- predict(gam_Dist_N_mpo, dat_fcast, type="response")
    dat_fcast$gam_Dist_mpo <- dat_fcast$presxD_mpo * exp(abundxD_mpo)
  }
  #Dist sampling - MPN
  if("mpn" %in% sampling){
    print("Fitting GAM-MPN")
    gam_Dist_P_mpn <- gam(pres~s(temp) + s(mld) + s(chl_surface), data=dat_hist_Dist_mpn, family=binomial)
    #plot(gam_Dist_P_mpn, pages=1)
    gam_Dist_N_mpn <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface), data=dat_hist_Dist_mpn[dat_hist_Dist_mpn$abundance>0,], family=gaussian)
    #plot(gam_Dist_N_mpn, pages=1)
    
    dat_hist$presxD_mpn <- predict(gam_Dist_P_mpn, dat_hist, type="response")
    abundxD_mpn <- predict(gam_Dist_N_mpn, dat_hist, type="response")
    dat_hist$gam_Dist_mpn <- dat_hist$presxD_mpn * exp(abundxD_mpn)
    
    dat_fcast$presxD_mpn <- predict(gam_Dist_P_mpn, dat_fcast, type="response")
    abundxD_mpn <- predict(gam_Dist_N_mpn, dat_fcast, type="response")
    dat_fcast$gam_Dist_mpn <- dat_fcast$presxD_mpn * exp(abundxD_mpn)
  }
  #Dist sampling - SPO
  if("spo" %in% sampling){
    print("Fitting GAM-SPO")
    gam_Dist_P_spo <- gam(pres~s(temp) + s(mld) + s(chl_surface), data=dat_hist_Dist_spo, family=binomial)
    #plot(gam_Dist_P_spo, pages=1)
    gam_Dist_N_spo <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface), data=dat_hist_Dist_spo[dat_hist_Dist_spo$abundance>0,], family=gaussian)
    #plot(gam_Dist_N_spo, pages=1)
    
    dat_hist$presxD_spo <- predict(gam_Dist_P_spo, dat_hist, type="response")
    abundxD_spo <- predict(gam_Dist_N_spo, dat_hist, type="response")
    dat_hist$gam_Dist_spo <- dat_hist$presxD_spo * exp(abundxD_spo)
    
    dat_fcast$presxD_spo <- predict(gam_Dist_P_spo, dat_fcast, type="response")
    abundxD_spo <- predict(gam_Dist_N_spo, dat_fcast, type="response")
    dat_fcast$gam_Dist_spo <- dat_fcast$presxD_spo * exp(abundxD_spo)
  }

  #Dist sampling - SPN
  if("spn" %in% sampling){
    print("Fitting GAM-SPN")
    gam_Dist_P_spn <- gam(pres~s(temp) + s(mld) + s(chl_surface), data=dat_hist_Dist_spn, family=binomial)
    #plot(gam_Dist_P_spn, pages=1)
    gam_Dist_N_spn <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface), data=dat_hist_Dist_spn[dat_hist_Dist_spn$abundance>0,], family=gaussian)
    #plot(gam_Dist_N_spn, pages=1)
    
    dat_hist$presxD_spn <- predict(gam_Dist_P_spn, dat_hist, type="response")
    abundxD_spn <- predict(gam_Dist_N_spn, dat_hist, type="response")
    dat_hist$gam_Dist_spn <- dat_hist$presxD_spn * exp(abundxD_spn)
    
    dat_fcast$presxD_spn <- predict(gam_Dist_P_spn, dat_fcast, type="response")
    abundxD_spn <- predict(gam_Dist_N_spn, dat_fcast, type="response")
    dat_fcast$gam_Dist_spn <- dat_fcast$presxD_spn * exp(abundxD_spn)
  }
  #Dist sampling - ALLO
  if("allo" %in% sampling){
    print("Fitting GAM-ALLO")
    gam_Dist_P_allo <- gam(pres~s(temp) + s(mld) + s(chl_surface), data=dat_hist_Dist_allo, family=binomial)
    #plot(gam_Dist_P_allo, pages=1)
    gam_Dist_N_allo <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface), data=dat_hist_Dist_allo[dat_hist_Dist_allo$abundance>0,], family=gaussian)
    #plot(gam_Dist_N_allo, pages=1)
    
    dat_hist$presxD_allo <- predict(gam_Dist_P_allo, dat_hist, type="response")
    abundxD_allo <- predict(gam_Dist_N_allo, dat_hist, type="response")
    dat_hist$gam_Dist_allo <- dat_hist$presxD_allo * exp(abundxD_allo)
    
    dat_fcast$presxD_allo <- predict(gam_Dist_P_allo, dat_fcast, type="response")
    abundxD_allo <- predict(gam_Dist_N_allo, dat_fcast, type="response")
    dat_fcast$gam_Dist_allo <- dat_fcast$presxD_allo * exp(abundxD_allo)
  }
  #Dist sampling - Alln
  if("alln" %in% sampling){
    print("Fitting GAM-ALLN")
    gam_Dist_P_alln <- gam(pres~s(temp) + s(mld) + s(chl_surface), data=dat_hist_Dist_alln, family=binomial)
    #plot(gam_Dist_P_alln, pages=1)
    gam_Dist_N_alln <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface), data=dat_hist_Dist_alln[dat_hist_Dist_alln$abundance>0,], family=gaussian)
    #plot(gam_Dist_N_alln, pages=1)
    
    dat_hist$presxD_alln <- predict(gam_Dist_P_alln, dat_hist, type="response")
    abundxD_alln <- predict(gam_Dist_N_alln, dat_hist, type="response")
    dat_hist$gam_Dist_alln <- dat_hist$presxD_alln * exp(abundxD_alln)
    
    dat_fcast$presxD_alln <- predict(gam_Dist_P_alln, dat_fcast, type="response")
    abundxD_alln <- predict(gam_Dist_N_alln, dat_fcast, type="response")
    dat_fcast$gam_Dist_alln <- dat_fcast$presxD_alln * exp(abundxD_alln)
  }
#
######BY sampling
#
  #Bycatch 0.5...
  if("BY_1" %in% sampling){
    print("Fitting GAM-BY_1")
    gam_BY_P_1 <- gam(pres~s(temp) + s(mld) + s(chl_surface), data=dat_hist_BY_1, family=binomial)
    #plot(gam_BY_P, pages=1)
    gam_BY_N_1 <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface), data=dat_hist_BY_1[dat_hist_BY_1$abundance>0,], family=gaussian)
    #plot(gam_BY_N, pages=1)
    
    dat_hist$presxB_1 <- predict(gam_BY_P_1, dat_hist, type="response")
    abundxB_1 <- predict(gam_BY_N_1, dat_hist, type="response")
    dat_hist$gam_BY_1 <- dat_hist$presxB_1 * exp(abundxB_1)
    
    dat_fcast$presxB_1 <- predict(gam_BY_P_1, dat_fcast, type="response")
    abundxB_1 <- predict(gam_BY_N_1, dat_fcast, type="response")
    dat_fcast$gam_BY_1 <- dat_fcast$presxB_1 * exp(abundxB_1)
  }

  #Bycatch 0.7...
  if("BY_2" %in% sampling){
    print("Fitting GAM-BY_2")
    gam_BY_P_2 <- gam(pres~s(temp) + s(mld) + s(chl_surface), data=dat_hist_BY_2, family=binomial)
    #plot(gam_BY_P, pages=1)
    gam_BY_N_2 <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface), data=dat_hist_BY_2[dat_hist_BY_2$abundance>0,], family=gaussian)
    #plot(gam_BY_N, pages=1)
    
    dat_hist$presxB_2 <- predict(gam_BY_P_2, dat_hist, type="response")
    abundxB_2 <- predict(gam_BY_N_2, dat_hist, type="response")
    dat_hist$gam_BY_2 <- dat_hist$presxB_2 * exp(abundxB_2)
    
    dat_fcast$presxB_2 <- predict(gam_BY_P_2, dat_fcast, type="response")
    abundxB_2 <- predict(gam_BY_N_2, dat_fcast, type="response")
    dat_fcast$gam_BY_2 <- dat_fcast$presxB_2 * exp(abundxB_2)
  }

  #Bycatch 0.9...
  if("BY_3" %in% sampling){
    print("Fitting GAM-BY_3")
    gam_BY_P_3 <- gam(pres~s(temp) + s(mld) + s(chl_surface), data=dat_hist_BY_3, family=binomial)
    #plot(gam_BY_P, pages=1)
    gam_BY_N_3 <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface), data=dat_hist_BY_3[dat_hist_BY_3$abundance>0,], family=gaussian)
    #plot(gam_BY_N, pages=1)
    
    dat_hist$presxB_3 <- predict(gam_BY_P_3, dat_hist, type="response")
    abundxB_3 <- predict(gam_BY_N_3, dat_hist, type="response")
    dat_hist$gam_BY_3 <- dat_hist$presxB_3 * exp(abundxB_3)
    
    dat_fcast$presxB_3 <- predict(gam_BY_P_3, dat_fcast, type="response")
    abundxB_3 <- predict(gam_BY_N_3, dat_fcast, type="response")
    dat_fcast$gam_BY_3 <- dat_fcast$presxB_3 * exp(abundxB_3)
  }
#
#####Closed Area Sampling
#
  #Closed Area - Small
  if("CA_sm" %in% sampling) {
    print("Fitting GAM-CASM")
    gam_CA_P_sm <- gam(pres~s(temp) + s(mld) + s(chl_surface), data=dat_hist_CA_sm, family=binomial)
    #plot(gam_CA_P_sm, pages=1)
    gam_CA_N_sm <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface), data=dat_hist_CA_sm[dat_hist_CA_sm$abundance>0,], family=gaussian)
    #plot(gam_CA_N_sm, pages=1)
    
    dat_hist$presxCA_sm <- predict(gam_CA_P_sm, dat_hist, type="response")
    abundxCA_sm <- predict(gam_CA_N_sm, dat_hist, type="response")
    dat_hist$gam_CA_sm <- dat_hist$presxCA_sm * exp(abundxCA_sm)
    
    dat_fcast$presxCA_sm <- predict(gam_CA_P_sm, dat_fcast, type="response")
    abundxCA_sm <- predict(gam_CA_N_sm, dat_fcast, type="response")
    dat_fcast$gam_CA_sm <- dat_fcast$presxCA_sm * exp(abundxCA_sm)
  }
  #Closed Area - Medium
  if("CA_med" %in% sampling) {
    print("Fitting GAM-CAMED")
    gam_CA_P_med <- gam(pres~s(temp) + s(mld) + s(chl_surface), data=dat_hist_CA_med, family=binomial)
    #plot(gam_CA_P_med, pages=1)
    gam_CA_N_med <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface), data=dat_hist_CA_med[dat_hist_CA_med$abundance>0,], family=gaussian)
    #plot(gam_CA_N_med, pages=1)
    
    dat_hist$presxCA_med <- predict(gam_CA_P_med, dat_hist, type="response")
    abundxCA_med <- predict(gam_CA_N_med, dat_hist, type="response")
    dat_hist$gam_CA_med <- dat_hist$presxCA_med * exp(abundxCA_med)
    
    dat_fcast$presxCA_med <- predict(gam_CA_P_med, dat_fcast, type="response")
    abundxCA_med <- predict(gam_CA_N_med, dat_fcast, type="response")
    dat_fcast$gam_CA_med <- dat_fcast$presxCA_med * exp(abundxCA_med)
  }
  #Closed Area - Large
  if("CA_lar" %in% sampling) {
    print("Fitting GAM-CALAR")
    gam_CA_P_lar <- gam(pres~s(temp) + s(mld) + s(chl_surface), data=dat_hist_CA_lar, family=binomial)
    #plot(gam_CA_P_lar, pages=1)
    gam_CA_N_lar <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface), data=dat_hist_CA_lar[dat_hist_CA_lar$abundance>0,], family=gaussian)
    #plot(gam_CA_N_lar, pages=1)
    
    dat_hist$presxCA_lar <- predict(gam_CA_P_lar, dat_hist, type="response")
    abundxCA_lar <- predict(gam_CA_N_lar, dat_hist, type="response")
    dat_hist$gam_CA_lar <- dat_hist$presxCA_lar * exp(abundxCA_lar)
    
    dat_fcast$presxCA_lar <- predict(gam_CA_P_lar, dat_fcast, type="response")
    abundxCA_lar <- predict(gam_CA_N_lar, dat_fcast, type="response")
    dat_fcast$gam_CA_lar <- dat_fcast$presxCA_lar * exp(abundxCA_lar)
  }


#### PLOTS ####
#NOTE: if you are running individual models (not all of them) then will need to
# "#" the below out and turn on the relevant lines for the specific models
# being run above for it to work 

#Env Cov 1
  #Pres-Abs
  par(mfrow=c(3,5))
  plot(gam_Ran_P, select=1, main="Random, Pres-Abs", cex.main=0.8, scale=0)
  plot(gam_Tar_P_1, select=1, main="Preferential, Pres-Abs", cex.main=0.8, scale=0)
 # plot(gam_Tar_P_2, select=1, main="gam_Tar_0.6, Pres-Abs", scale=0)
  #plot(gam_Tar_P_3, select=1, main="gam_Tar_0.7, Pres-Abs", scale=0)
 # plot(gam_Tar_P_4, select=1, main="gam_Tar_0.8, Pres-Abs", scale=0)
  #plot(gam_Tar_P_5, select=1, main="gam_Tar_0.9, Pres-Abs", scale=0)
  plot(gam_Dist_P_npo, select=1, main="Port Northern Offshore, Pres-Abs", cex.main=0.8, scale=0)
  plot(gam_Dist_P_npn, select=1, main="Port Northern Nearshore, Pres-Abs", cex.main=0.8,scale=0)
  plot(gam_Dist_P_spo, select=1, main="Port Southern Offshore, Pres-Abs", cex.main=0.8,scale=0)
  plot(gam_Dist_P_spn, select=1, main="Port Southern Nearshore, Pres-Abs",cex.main=0.8, scale=0)
  plot(gam_Dist_P_mpo, select=1, main="Port Middle Offshore, Pres-Abs", cex.main=0.8,scale=0)
  plot(gam_Dist_P_mpn, select=1, main="Port Middle Nearshore, Pres-Abs",cex.main=0.8, scale=0)
  plot(gam_Dist_P_allo, select=1, main="Port All Offshore, Pres-Abs", cex.main=0.8,scale=0)
  plot(gam_Dist_P_alln, select=1, main="Port All Nearshore, Pres-Abs",cex.main=0.8, scale=0)
  plot(gam_BY_P_2, select=1, main="Bycatch, Pres-Abs", cex.main=0.8, scale=0)
  plot(gam_CA_P_sm, select=1, main="Closed Area Small, Pres-Abs", cex.main=0.8,scale=0)
  plot(gam_CA_P_med, select=1, main="Closed Area Medium, Pres-Abs", cex.main=0.8,scale=0)
  plot(gam_CA_P_lar, select=1, main="Closed Area Large, Pres-Abs", cex.main=0.8,scale=0)
  
  #Abund
  par(mfrow=c(4,4))
  plot(gam_Ran_N, select=1, main="Random, Abund", cex.main=0.8, scale=0)
  plot(gam_Tar_N_1, select=1, main="Preferential, Abund", cex.main=0.8, scale=0)
  # plot(gam_Tar_N_2, select=1, main="gam_Tar_0.6, Abund", scale=0)
  #plot(gam_Tar_N_3, select=1, main="gam_Tar_0.7, Abund", scale=0)
  # plot(gam_Tar_N_4, select=1, main="gam_Tar_0.8, Abund", scale=0)
  #plot(gam_Tar_N_5, select=1, main="gam_Tar_0.9, Abund", scale=0)
  plot(gam_Dist_N_npo, select=1, main="Port Northern Offshore, Abund", cex.main=0.8, scale=0)
  plot(gam_Dist_N_npn, select=1, main="Port Northern Nearshore, Abund", cex.main=0.8,scale=0)
  plot(gam_Dist_N_spo, select=1, main="Port Southern Offshore, Abund", cex.main=0.8,scale=0)
  plot(gam_Dist_N_spn, select=1, main="Port Southern Nearshore, Abund",cex.main=0.8, scale=0)
  plot(gam_Dist_N_mpo, select=1, main="Port Middle Offshore, Abund", cex.main=0.8,scale=0)
  plot(gam_Dist_N_mpn, select=1, main="Port Middle Nearshore, Abund",cex.main=0.8, scale=0)
  plot(gam_Dist_N_allo, select=1, main="Port All Offshore, Abund", cex.main=0.8,scale=0)
  plot(gam_Dist_N_alln, select=1, main="Port All Nearshore, Abund",cex.main=0.8, scale=0)
  plot(gam_BY_N_2, select=1, main="Bycatch, Abund", cex.main=0.8, scale=0)
  plot(gam_CA_N_sm, select=1, main="Closed Area Small, Abund", cex.main=0.8,scale=0)
  plot(gam_CA_N_med, select=1, main="Closed Area Medium, Abund", cex.main=0.8,scale=0)
  plot(gam_CA_N_lar, select=1, main="Closed Area Large, Abund", cex.main=0.8,scale=0)
  
#Env Cov 2
  #Pres-Abs
  par(mfrow=c(3,5))
  plot(gam_Ran_P, select=2, main="gam_Ran, Pres-Abs", scale=0)
  plot(gam_Tar_P_1, select=2, main="gam_Tar_0.5, Pres-Abs", scale=0)
  # plot(gam_Tar_P_2, select=2, main="gam_Tar_0.6, Pres-Abs", scale=0)
  # plot(gam_Tar_P_3, select=2, main="gam_Tar_0.7, Pres-Abs", scale=0)
  # plot(gam_Tar_P_4, select=2, main="gam_Tar_0.8, Pres-Abs", scale=0)
  # plot(gam_Tar_P_5, select=2, main="gam_Tar_0.9, Pres-Abs", scale=0)
  plot(gam_Dist_P_npo, select=2, main="gam_Dist_npo, Pres-Abs", scale=0)
  plot(gam_Dist_P_npn, select=2, main="gam_Dist_npn, Pres-Abs", scale=0)
  plot(gam_Dist_P_spo, select=2, main="gam_Dist_spo, Pres-Abs", scale=0)
  plot(gam_Dist_P_mpo, select=2, main="gam_Dist_mpo, Pres-Abs", scale=0)
  plot(gam_Dist_P_mpn, select=2, main="gam_Dist_mpn, Pres-Abs", scale=0)
  plot(gam_Dist_P_allo, select=2, main="gam_Dist_allo, Pres-Abs", scale=0)
  plot(gam_Dist_P_alln, select=2, main="gam_Dist_alln, Pres-Abs", scale=0)
  plot(gam_BY_P_2, select=2, main="gam_BY, Pres-Abs", scale=0)
  plot(gam_CA_P_sm, select=2, main="gam_CA_sm, Pres-Abs", scale=0)
  plot(gam_CA_P_med, select=2, main="gam_CA_med, Pres-Abs", scale=0)
  plot(gam_CA_P_lar, select=2, main="gam_CA_lar, Pres-Abs", scale=0)
  
  #Abund
  par(mfrow=c(3,5))
  plot(gam_Ran_N, select=2, main="gam_Ran, Abund", scale=0)
  plot(gam_Tar_N_1, select=2, main="gam_Tar_0.5, Abund", scale=0)
  #plot(gam_Tar_N_2, select=2, main="gam_Tar_0.6, Abund", scale=0)
  #plot(gam_Tar_N_3, select=2, main="gam_Tar_0.7, Abund", scale=0)
  #plot(gam_Tar_N_4, select=2, main="gam_Tar_0.8, Abund", scale=0)
  #plot(gam_Tar_N_5, select=2, main="gam_Tar_0.9, Abund", scale=0)
  plot(gam_Dist_N_npo, select=2, main="gam_Dist_npo, Abund", scale=0)
  plot(gam_Dist_N_npn, select=2, main="gam_Dist_npn, Abund", scale=0)
  plot(gam_Dist_N_spo, select=2, main="gam_Dist_spo, Abund", scale=0)
  plot(gam_Dist_N_spn, select=2, main="gam_Dist_spo, Abund", scale=0)
  plot(gam_Dist_N_mpo, select=2, main="gam_Dist_mpo, Abund", scale=0)
  plot(gam_Dist_N_mpn, select=2, main="gam_Dist_mpn, Abund", scale=0)
  plot(gam_Dist_N_allo, select=2, main="gam_Dist_allo, Abund", scale=0)
  plot(gam_Dist_N_alln, select=2, main="gam_Dist_alln, Abund", scale=0)
  plot(gam_BY_N_2, select=2, main="gam_BY, Abund", scale=0)
  plot(gam_CA_N_sm, select=2, main="gam_CA_sm, Abund", scale=0)
  plot(gam_CA_N_med, select=2, main="gam_CA_med, Abund", scale=0)
  plot(gam_CA_N_lar, select=2, main="gam_CA_lar, Abund", scale=0)
  
#Env Cov 3
  #Pres-Abs
  par(mfrow=c(3,5))
  plot(gam_Ran_P, select=3, main="gam_Ran, Pres-Abs", scale=0)
  plot(gam_Tar_P_1, select=3, main="gam_Tar_0.5, Pres-Abs", scale=0)
  # plot(gam_Tar_P_2, select=3, main="gam_Tar_0.6, Pres-Abs", scale=0)
  # plot(gam_Tar_P_3, select=3, main="gam_Tar_0.7, Pres-Abs", scale=0)
  # plot(gam_Tar_P_4, select=3, main="gam_Tar_0.8, Pres-Abs", scale=0)
  # plot(gam_Tar_P_5, select=3, main="gam_Tar_0.9, Pres-Abs", scale=0)
  plot(gam_Dist_P_npo, select=3, main="gam_Dist_npo, Pres-Abs", scale=0)
  plot(gam_Dist_P_npn, select=3, main="gam_Dist_npn, Pres-Abs", scale=0)
  plot(gam_Dist_P_spo, select=3, main="gam_Dist_spo, Pres-Abs", scale=0)
  plot(gam_Dist_P_mpo, select=3, main="gam_Dist_mpo, Pres-Abs", scale=0)
  plot(gam_Dist_P_mpn, select=3, main="gam_Dist_mpn, Pres-Abs", scale=0)
  plot(gam_Dist_P_allo, select=3, main="gam_Dist_allo, Pres-Abs", scale=0)
  plot(gam_Dist_P_alln, select=3, main="gam_Dist_alln, Pres-Abs", scale=0)
  plot(gam_BY_P_2, select=3, main="gam_BY, Pres-Abs", scale=0)
  plot(gam_CA_P_sm, select=3, main="gam_CA_sm, Pres-Abs", scale=0)
  plot(gam_CA_P_med, select=3, main="gam_CA_med, Pres-Abs", scale=0)
  plot(gam_CA_P_lar, select=3, main="gam_CA_lar, Pres-Abs", scale=0)
  
  #Abund
  par(mfrow=c(3,5))
  plot(gam_Ran_N, select=3, main="gam_Ran, Abund", scale=0)
  plot(gam_Tar_N_1, select=3, main="gam_Tar_0.5, Abund", scale=0)
  # plot(gam_Tar_N_2, select=3, main="gam_Tar_0.6, Abund", scale=0)
  # plot(gam_Tar_N_3, select=3, main="gam_Tar_0.7, Abund", scale=0)
  # plot(gam_Tar_N_4, select=3, main="gam_Tar_0.8, Abund", scale=0)
  # plot(gam_Tar_N_5, select=3, main="gam_Tar_0.9, Abund", scale=0)
  plot(gam_Dist_N_npo, select=3, main="gam_Dist_npo, Abund", scale=0)
  plot(gam_Dist_N_npn, select=3, main="gam_Dist_npn, Abund", scale=0)
  plot(gam_Dist_N_spo, select=3, main="gam_Dist_spo, Abund", scale=0)
  plot(gam_Dist_N_mpo, select=3, main="gam_Dist_mpo, Abund", scale=0)
  plot(gam_Dist_N_mpn, select=3, main="gam_Dist_mpn, Abund", scale=0)
  plot(gam_Dist_N_allo, select=3, main="gam_Dist_allo, Abund", scale=0)
  plot(gam_Dist_N_alln, select=3, main="gam_Dist_alln, Abund", scale=0)
  plot(gam_BY_N_2, select=3, main="gam_BY, Abund", scale=0)
  plot(gam_CA_N_sm, select=3, main="gam_CA_sm, Abund", scale=0)
  plot(gam_CA_N_med, select=3, main="gam_CA_med, Abund", scale=0)
  plot(gam_CA_N_lar, select=3, main="gam_CA_lar, Abund", scale=0)
  
  