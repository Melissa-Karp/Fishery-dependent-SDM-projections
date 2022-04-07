### Fishery-dependent SDM (L^3) Estimation Model Code
## Function to run GAMs on data from OM
## Function returns only the fitted and predicted values
#modified my M.Karp 1/07/20


###################################################
#   Full Models- with Chla and space-time Tensor  #
###################################################

mte <- list(c(1,.5), 1)  

#
##### Random sampling ####
#
if("ran" %in% sampling){
  print("Fitting GAM-Ran-te")
  gam_Ran_Pte <- gam(pres~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_random, family=binomial) #note - chl-surface is log of chl-surface
  #plot(gam_Ran_Pte, pages=1)
  gam_Ran_Nte <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_random[dat_hist_random$abundance>0,], family=gaussian)
  #plot(gam_Ran_Nte, pages=1)

  dat_hist$presxRte <- predict(gam_Ran_Pte, dat_hist, type="response")
  abundxRte <- predict(gam_Ran_Nte, dat_hist, type="response") 
  dat_hist$gam_Ran_te <- dat_hist$presxRte * exp(abundxRte)  #predicted catch
  
  dat_fcast$presxRte <- predict(gam_Ran_Pte, dat_fcast, type="response")
  abundxRte <- predict(gam_Ran_Nte, dat_fcast, type="response")
  dat_fcast$gam_Ran_te <- dat_fcast$presxRte * exp(abundxRte)
}

#
##### Preferential Sampling #######
#
#Pref sampling - 0.5
if("tar_0.5" %in% sampling) {
  print("Fitting GAM-Pref0.5_te")
  gam_Tar_P_1te<- gam(pres~s(temp)+ s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Tar_1, family=binomial)
  #plot(gam_Tar_P_1te, pages=1)
  gam_Tar_N_1te <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface)+ te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Tar_1[dat_hist_Tar_1$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_1te, pages=1)
  
  dat_hist$presxT_1_te <- predict(gam_Tar_P_1te, dat_hist, type="response")
  abundxT_1te <- predict(gam_Tar_N_1te, dat_hist, type="response")
  dat_hist$gam_Tar_0.5_te <- dat_hist$presxT_1_te * exp(abundxT_1te)
  
  dat_fcast$presxT_1_te <- predict(gam_Tar_P_1te, dat_fcast, type="response")
  abundxT_1te <- predict(gam_Tar_N_1te, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.5_te <- dat_fcast$presxT_1_te * exp(abundxT_1te)
}
#Pref sampling -0.6
if("tar_0.6" %in% sampling) {
  print("Fitting GAM-Pref0.6_te")
  gam_Tar_P_2te <- gam(pres~s(temp)+ s(mld) + s(chl_surface)+ te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Tar_2, family=binomial)
  #plot(gam_Tar_P_2te, pages=1)
  gam_Tar_N_2te <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface)+ te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Tar_2[dat_hist_Tar_2$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_2te, pages=1)
  
  dat_hist$presxT_2te <- predict(gam_Tar_P_2te, dat_hist, type="response")
  abundxT_2te <- predict(gam_Tar_N_2te, dat_hist, type="response")
  dat_hist$gam_Tar_0.6_te <- dat_hist$presxT_2te * exp(abundxT_2te)
  
  dat_fcast$presxT_2te <- predict(gam_Tar_P_2te, dat_fcast, type="response")
  abundxT_2te <- predict(gam_Tar_N_2te, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.6_te <- dat_fcast$presxT_2te * exp(abundxT_2te)
}
#Pref sampling -0.7
if("tar_0.7" %in% sampling) {
  print("Fitting GAM-Pref0.7_te")
  gam_Tar_P_3te <- gam(pres~s(temp)+ s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Tar_3, family=binomial)
  #plot(gam_Tar_P_3te, pages=1)
  gam_Tar_N_3te <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface)+ te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Tar_3[dat_hist_Tar_3$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_3te, pages=1)
  
  dat_hist$presxT_3te <- predict(gam_Tar_P_3te, dat_hist, type="response")
  abundxT_3te <- predict(gam_Tar_N_3te, dat_hist, type="response")
  dat_hist$gam_Tar_0.7_te <- dat_hist$presxT_3te * exp(abundxT_3te)
  
  dat_fcast$presxT_3te <- predict(gam_Tar_P_3te, dat_fcast, type="response")
  abundxT_3te <- predict(gam_Tar_N_3te, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.7_te <- dat_fcast$presxT_3te * exp(abundxT_3te)
}
#Pref sampling -0.8
if("tar_0.8" %in% sampling) {
  print("Fitting GAM-Pref0.8_te")
  gam_Tar_P_4te <- gam(pres~s(temp)+ s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Tar_4, family=binomial)
  #plot(gam_Tar_P_4te, pages=1)
  gam_Tar_N_4te <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Tar_4[dat_hist_Tar_4$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_4te, pages=1)
  
  dat_hist$presxT_4te <- predict(gam_Tar_P_4te, dat_hist, type="response")
  abundxT_4te <- predict(gam_Tar_N_4te, dat_hist, type="response")
  dat_hist$gam_Tar_0.8_te <- dat_hist$presxT_4te * exp(abundxT_4te)
  
  dat_fcast$presxT_4te <- predict(gam_Tar_P_4te, dat_fcast, type="response")
  abundxT_4te <- predict(gam_Tar_N_4te, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.8_te <- dat_fcast$presxT_4te * exp(abundxT_4te)
}
#Pref sampling -0.9
if("tar_0.9" %in% sampling) {
  print("Fitting GAM-Pref0.9_te")
  gam_Tar_P_5te <- gam(pres~s(temp)+ s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Tar_5, family=binomial)
  #plot(gam_Tar_P_5te, pages=1)
  gam_Tar_N_5te <- gam(log_abundance~s(temp)+ s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Tar_5[dat_hist_Tar_5$abundance>0,], family=gaussian)
  #plot(gam_Tar_N_5te, pages=1)
  
  dat_hist$presxT_5te <- predict(gam_Tar_P_5te, dat_hist, type="response")
  abundxT_5te <- predict(gam_Tar_N_5te, dat_hist, type="response")
  dat_hist$gam_Tar_0.9_te <- dat_hist$presxT_5te * exp(abundxT_5te)
  
  dat_fcast$presxT_5te <- predict(gam_Tar_P_5te, dat_fcast, type="response")
  abundxT_5te <- predict(gam_Tar_N_5te, dat_fcast, type="response")
  dat_fcast$gam_Tar_0.9_te <- dat_fcast$presxT_5te * exp(abundxT_5te)
}

#
##### Distance from Port Sampling #######
#
#Dist sampling - NPO
if("npo" %in% sampling){
  print("Fitting GAM-NPO_te")
  gam_Dist_P_npote <- gam(pres~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Dist_npo, family=binomial)
  #plot(gam_Dist_P_npote, pages=1)
  gam_Dist_N_npote <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Dist_npo[dat_hist_Dist_npo$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_npote, pages=1)
  
  dat_hist$presxD_npote <- predict(gam_Dist_P_npote, dat_hist, type="response")
  abundxD_npote <- predict(gam_Dist_N_npote, dat_hist, type="response")
  dat_hist$gam_Dist_npo_te <- dat_hist$presxD_npote * exp(abundxD_npote)
  
  dat_fcast$presxD_npote <- predict(gam_Dist_P_npote, dat_fcast, type="response")
  abundxD_npote <- predict(gam_Dist_N_npote, dat_fcast, type="response")
  dat_fcast$gam_Dist_npo_te <- dat_fcast$presxD_npote * exp(abundxD_npote)
}
#Dist sampling - NPN
if("npn" %in% sampling){
  print("Fitting GAM-NPN_te")
  gam_Dist_P_npnte <- gam(pres~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Dist_npn, family=binomial)
  #plot(gam_Dist_P_npnte, pages=1)
  gam_Dist_N_npnte <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Dist_npn[dat_hist_Dist_npn$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_npnte, pages=1)
  
  dat_hist$presxD_npnte <- predict(gam_Dist_P_npnte, dat_hist, type="response")
  abundxD_npnte <- predict(gam_Dist_N_npnte, dat_hist, type="response")
  dat_hist$gam_Dist_npn_te <- dat_hist$presxD_npnte * exp(abundxD_npnte)
  
  dat_fcast$presxD_npnte <- predict(gam_Dist_P_npnte, dat_fcast, type="response")
  abundxD_npnte <- predict(gam_Dist_N_npnte, dat_fcast, type="response")
  dat_fcast$gam_Dist_npn_te <- dat_fcast$presxD_npnte * exp(abundxD_npnte)
}
#Dist sampling - MPO
if("mpo" %in% sampling){
  print("Fitting GAM-MPO_te")
  gam_Dist_P_mpote <- gam(pres~s(temp) + s(mld) + s(chl_surface)+ te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Dist_mpo, family=binomial)
  #plot(gam_Dist_P_mpote, pages=1)
  gam_Dist_N_mpote <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface)+ te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Dist_mpo[dat_hist_Dist_mpo$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_mpote, pages=1)
  
  dat_hist$presxD_mpote <- predict(gam_Dist_P_mpote, dat_hist, type="response")
  abundxD_mpote <- predict(gam_Dist_N_mpote, dat_hist, type="response")
  dat_hist$gam_Dist_mpo_te <- dat_hist$presxD_mpote * exp(abundxD_mpote)
  
  dat_fcast$presxD_mpote <- predict(gam_Dist_P_mpote, dat_fcast, type="response")
  abundxD_mpote <- predict(gam_Dist_N_mpote, dat_fcast, type="response")
  dat_fcast$gam_Dist_mpo_te <- dat_fcast$presxD_mpote * exp(abundxD_mpote)
}
#Dist sampling - MPN
if("mpn" %in% sampling){
  print("Fitting GAM-MPN_te")
  gam_Dist_P_mpnte <- gam(pres~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Dist_mpn, family=binomial)
  #plot(gam_Dist_P_mpnte, pages=1)
  gam_Dist_N_mpnte <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Dist_mpn[dat_hist_Dist_mpn$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_mpnte, pages=1)
  
  dat_hist$presxD_mpnte <- predict(gam_Dist_P_mpnte, dat_hist, type="response")
  abundxD_mpnte <- predict(gam_Dist_N_mpnte, dat_hist, type="response")
  dat_hist$gam_Dist_mpn_te <- dat_hist$presxD_mpnte * exp(abundxD_mpnte)
  
  dat_fcast$presxD_mpnte <- predict(gam_Dist_P_mpnte, dat_fcast, type="response")
  abundxD_mpnte <- predict(gam_Dist_N_mpnte, dat_fcast, type="response")
  dat_fcast$gam_Dist_mpn_te <- dat_fcast$presxD_mpnte * exp(abundxD_mpnte)
}
#Dist sampling - SPO
if("spo" %in% sampling){
  print("Fitting GAM-SPO_te")
  gam_Dist_P_spote <- gam(pres~s(temp) + s(mld) + s(chl_surface)+ te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Dist_spo, family=binomial)
  #plot(gam_Dist_P_spote, pages=1)
  gam_Dist_N_spote <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface)+ te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Dist_spo[dat_hist_Dist_spo$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_spote, pages=1)
  
  dat_hist$presxD_spote <- predict(gam_Dist_P_spote, dat_hist, type="response")
  abundxD_spote <- predict(gam_Dist_N_spote, dat_hist, type="response")
  dat_hist$gam_Dist_spo_te <- dat_hist$presxD_spote * exp(abundxD_spote)
  
  dat_fcast$presxD_spote <- predict(gam_Dist_P_spote, dat_fcast, type="response")
  abundxD_spote <- predict(gam_Dist_N_spote, dat_fcast, type="response")
  dat_fcast$gam_Dist_spo_te <- dat_fcast$presxD_spote * exp(abundxD_spote)
}
#Dist sampling - SPN
if("spn" %in% sampling){
  print("Fitting GAM-SPN_te")
  gam_Dist_P_spnte <- gam(pres~s(temp) + s(mld) + s(chl_surface)+ te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Dist_spn, family=binomial)
  #plot(gam_Dist_P_spnte, pages=1)
  gam_Dist_N_spnte <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface)+ te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Dist_spn[dat_hist_Dist_spn$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_spnte, pages=1)

  dat_hist$presxD_spnte <- predict(gam_Dist_P_spnte, dat_hist, type="response")
  abundxD_spnte <- predict(gam_Dist_N_spnte, dat_hist, type="response")
  dat_hist$gam_Dist_spn_te <- dat_hist$presxD_spnte * exp(abundxD_spnte)
  
  dat_fcast$presxD_spnte <- predict(gam_Dist_P_spnte, dat_fcast, type="response")
  abundxD_spnte <- predict(gam_Dist_N_spnte, dat_fcast, type="response")
  dat_fcast$gam_Dist_spn_te <- dat_fcast$presxD_spnte * exp(abundxD_spnte)
}
#Dist sampling - ALLO
if("allo" %in% sampling){
  print("Fitting GAM-ALLO_te")
  gam_Dist_P_allote <- gam(pres~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Dist_allo, family=binomial)
  #plot(gam_Dist_P_allote, pages=1)
  gam_Dist_N_allote <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Dist_allo[dat_hist_Dist_allo$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_allote, pages=1)
  
  dat_hist$presxD_allote <- predict(gam_Dist_P_allote, dat_hist, type="response")
  abundxD_allote <- predict(gam_Dist_N_allote, dat_hist, type="response")
  dat_hist$gam_Dist_allo_te <- dat_hist$presxD_allote * exp(abundxD_allote)
  
  dat_fcast$presxD_allote <- predict(gam_Dist_P_allote, dat_fcast, type="response")
  abundxD_allote <- predict(gam_Dist_N_allote, dat_fcast, type="response")
  dat_fcast$gam_Dist_allo_te <- dat_fcast$presxD_allote * exp(abundxD_allote)
}
#Dist sampling - Alln
if("alln" %in% sampling){
  print("Fitting GAM-ALLN_te")
  gam_Dist_P_allnte <- gam(pres~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Dist_alln, family=binomial)
  #plot(gam_Dist_P_allnte, pages=1)
  gam_Dist_N_allnte <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_Dist_alln[dat_hist_Dist_alln$abundance>0,], family=gaussian)
  #plot(gam_Dist_N_allnte, pages=1)
  
  dat_hist$presxD_allnte <- predict(gam_Dist_P_allnte, dat_hist, type="response")
  abundxD_allnte <- predict(gam_Dist_N_allnte, dat_hist, type="response")
  dat_hist$gam_Dist_alln_te <- dat_hist$presxD_allnte * exp(abundxD_allnte)
  
  dat_fcast$presxD_allnte <- predict(gam_Dist_P_allnte, dat_fcast, type="response")
  abundxD_allnte <- predict(gam_Dist_N_allnte, dat_fcast, type="response")
  dat_fcast$gam_Dist_alln_te <- dat_fcast$presxD_allnte * exp(abundxD_allnte)
}
#
######BY sampling
#
if("BY_2" %in% sampling){
  print("Fitting GAM-BY_te")
  gam_BY_Pte <- gam(pres~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_BY_2, family=binomial)
  #plot(gam_BY_Pte, pages=1)
  gam_BY_Nte <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_BY_2[dat_hist_BY_2$abundance>0,], family=gaussian)
  #plot(gam_BY_Nte, pages=1)
  
  dat_hist$presxBte <- predict(gam_BY_Pte, dat_hist, type="response")
  abundxBte <- predict(gam_BY_Nte, dat_hist, type="response")
  dat_hist$gam_BY_te <- dat_hist$presxBte * exp(abundxBte)
  
  dat_fcast$presxBte <- predict(gam_BY_Pte, dat_fcast, type="response")
  abundxBte <- predict(gam_BY_Nte, dat_fcast, type="response")
  dat_fcast$gam_BY_te <- dat_fcast$presxBte * exp(abundxBte)
}

#
#####Closed Area Sampling
#
#Closed Area - Small
if("CA_sm" %in% sampling) {
  print("Fitting GAM-CASM_te")
  gam_CA_P_smte <- gam(pres~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_CA_sm, family=binomial)
  #plot(gam_CA_P_smte, pages=1)
  gam_CA_N_smte <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_CA_sm[dat_hist_CA_sm$abundance>0,], family=gaussian)
  #plot(gam_CA_N_smte, pages=1)
  
  dat_hist$presxCA_smte <- predict(gam_CA_P_smte, dat_hist, type="response")
  abundxCA_smte <- predict(gam_CA_N_smte, dat_hist, type="response")
  dat_hist$gam_CA_sm_te <- dat_hist$presxCA_smte * exp(abundxCA_smte)
  
  dat_fcast$presxCA_smte <- predict(gam_CA_P_smte, dat_fcast, type="response")
  abundxCA_smte <- predict(gam_CA_N_smte, dat_fcast, type="response")
  dat_fcast$gam_CA_sm_te <- dat_fcast$presxCA_smte * exp(abundxCA_smte)
}
#Closed Area - Medium
if("CA_med" %in% sampling) {
  print("Fitting GAM-CAMED_te")
  gam_CA_P_medte <- gam(pres~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_CA_med, family=binomial)
  #plot(gam_CA_P_medte, pages=1)
  gam_CA_N_medte <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_CA_med[dat_hist_CA_med$abundance>0,], family=gaussian)
  #plot(gam_CA_N_medte, pages=1)
  
  dat_hist$presxCA_medte <- predict(gam_CA_P_medte, dat_hist, type="response")
  abundxCA_medte <- predict(gam_CA_N_medte, dat_hist, type="response")
  dat_hist$gam_CA_med_te <- dat_hist$presxCA_medte * exp(abundxCA_medte)
  
  dat_fcast$presxCA_medte <- predict(gam_CA_P_medte, dat_fcast, type="response")
  abundxCA_medte <- predict(gam_CA_N_medte, dat_fcast, type="response")
  dat_fcast$gam_CA_med_te <- dat_fcast$presxCA_medte * exp(abundxCA_medte)
}
#Closed Area - Large
if("CA_lar" %in% sampling) {
  print("Fitting GAM-CALAR_te")
  gam_CA_P_larte <- gam(pres~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_CA_lar, family=binomial)
  #plot(gam_CA_P_larte, pages=1)
  gam_CA_N_larte <- gam(log_abundance~s(temp) + s(mld) + s(chl_surface) + te(lat, lon, year, d=c(2,1), bs=c('ds','tp'), m=mte, k=c(30,10)), data=dat_hist_CA_lar[dat_hist_CA_lar$abundance>0,], family=gaussian)
  #plot(gam_CA_N_larte, pages=1)
  
  dat_hist$presxCA_larte <- predict(gam_CA_P_larte, dat_hist, type="response")
  abundxCA_larte <- predict(gam_CA_N_larte, dat_hist, type="response")
  dat_hist$gam_CA_lar_te <- dat_hist$presxCA_larte * exp(abundxCA_larte)
  
  dat_fcast$presxCA_larte <- predict(gam_CA_P_larte, dat_fcast, type="response")
  abundxCA_larte <- predict(gam_CA_N_larte, dat_fcast, type="response")
  dat_fcast$gam_CA_lar_te <- dat_fcast$presxCA_larte * exp(abundxCA_larte)
}