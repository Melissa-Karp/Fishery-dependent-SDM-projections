
setwd("~/DisMAP project/Location, Location, Location/Location Workshop")

#################################
##       Load Libraries        ##
#################################
library(dplyr)
library(mgcv)
library(ggplot2)
library(viridis)
library(raster)
library(gbm)
library(dismo)
library(tidyverse)
library(reshape2)


#################################
##       Data Simulation       ##
#################################

source("~/DisMAP project/Location, Location, Location/Location Workshop/SimulatedWorld_ROMS_FishDep_Final_updatedPrefs.R") #load ROMS simulation function
#dir <- "~/DisMAP project/Location, Location, Location/Location Workshop/ROMS/hadley_11_2_21" #directory where hadley ROMS data is stored (on dropbox, email steph for access)

dat <- SimulateWorld_ROMS_FishDepFun_Final(dir=dir, nsamples = 100) #takes a few mins
names(dat)[names(dat) == 'sst'] <- 'temp' #matching ROMS names. Quick temporary fix.
############################
#    Prepare Data Sets     #
############################

#Create dataframe with historical/forecast data
dat_hist_1 <- dat[dat$year<=2010,]
dat_hist<-dat_hist_1[dat_hist_1$year>1984,] ## so dat_hist used to fit the model is 1985-2010
dat_fcast <- dat[dat$year>2010,] #forecast using 2011-2100

dat_hist$log_abundance <- log(dat_hist$abundance)

dat_hist_random<-dat_hist[dat_hist$random_sampled>0,]
dat_hist_Tar_1<-dat_hist[dat_hist$pref_sampled>0,]
dat_hist_Dist_npo<-dat_hist[dat_hist$dist_sampled_npo>0,]
dat_hist_Dist_npn<-dat_hist[dat_hist$dist_sampled_npn>0,]
dat_hist_Dist_mpo<-dat_hist[dat_hist$dist_sampled_mpo>0,]
dat_hist_Dist_mpn<-dat_hist[dat_hist$dist_sampled_mpn>0,]
dat_hist_Dist_spo<-dat_hist[dat_hist$dist_sampled_spo>0,]
dat_hist_Dist_spn<-dat_hist[dat_hist$dist_sampled_spn>0,]
dat_hist_Dist_allo<-dat_hist[dat_hist$dist_sampled_allo>0,]
dat_hist_Dist_alln<-dat_hist[dat_hist$dist_sampled_alln>0,]
dat_hist_BY_2<-dat_hist[dat_hist$BY_sampled_2>0,]
dat_hist_CA_sm<-dat_hist[dat_hist$Closed_sampled_1>0,]
dat_hist_CA_med<-dat_hist[dat_hist$Closed_sampled_2>0,]
dat_hist_CA_lar<-dat_hist[dat_hist$Closed_sampled_3>0,]

###########################
#   MODEL FITTING         #
###########################
#sampling scenario options: ran, tar_0.5, npo, npn, mpo,
# mpn, spo, spn, allo, alln, BY_2, CA_sm, CA_med, CA_lar
#total of 14 different sampling scenarios/rules X 2 model algorithms 

#NOTE: if you are running individual models (not all of them) then will need to
# go into each fitting code and "#" the "PlOTS" section out and turn on the relevant 
#lines for the specific models to plot just those models being run for it to work 

## Environment Only 
 ### GAMS #####
  #### GAMS - Full Model #####
    sampling <- c("ran", "tar_0.5", "npo", "npn", 
                  "mpo", "mpn", "spo", "spn",
                  "allo", "alln", "BY_2", 
                  "CA_sm", "CA_med", "CA_lar")
    source("~/DisMAP project/Location, Location, Location/Location Workshop/Fitting_GAMs.R") 
    
    save(gam_Ran_N, gam_Ran_P, gam_Tar_N_1, gam_Tar_P_1, gam_Dist_N_npn, gam_Dist_P_npn, gam_Dist_N_npo,
     gam_Dist_P_npo, gam_Dist_N_mpn, gam_Dist_P_mpn, gam_Dist_N_mpo, gam_Dist_P_mpo, gam_Dist_N_spn,
     gam_Dist_P_spn, gam_Dist_N_spo, gam_Dist_P_spo, gam_Dist_N_alln, gam_Dist_P_alln, gam_Dist_N_allo,
     gam_Dist_P_allo, gam_BY_N_2, gam_BY_P_2, gam_CA_N_sm, gam_CA_P_sm, gam_CA_N_med,
     gam_CA_P_med, gam_CA_N_lar, gam_CA_P_lar, file="GAM_E_Full_data_11_8_21_had.RData")

  ###BRTS ####  
  #### Boosted Regression Trees (BRTs) - Full Model #####
    sampling <- c("ran", "tar_0.5", "npo", "npn", "mpo", 
                  "mpn", "spo", "spn","allo", "alln", "BY_2", 
                  "CA_sm", "CA_med", "CA_lar")
    source("~/DisMAP project/Location, Location, Location/Location Workshop/Fitting_BRTs.R") 
    
    save(brt_R_N, brt_R_P, brt_T_N_1, brt_T_P_1, brt_dist_N_npn, brt_dist_P_npn, brt_dist_N_npo,
         brt_dist_P_npo, brt_dist_N_mpn, brt_dist_P_mpn, brt_dist_N_mpo, brt_dist_P_mpo, brt_dist_N_spn,
         brt_dist_P_spn, brt_dist_N_spo, brt_dist_P_spo, brt_dist_N_alln, brt_dist_P_alln, brt_dist_N_allo,
         brt_dist_P_allo, brt_B_N_2, brt_B_P_2, brt_CA_N_sm, brt_CA_P_sm, brt_CA_N_med, brt_CA_P_med,
         brt_CA_N_lar, brt_CA_P_lar, file="BRTs_E_Full_data_11_8_21_had.RData")

  ## Add spacetime variables
    #GAM-spacetime
    sampling <- c("ran", "tar_0.5", "npo", "npn", "mpo", 
                  "mpn", "spo", "spn","allo", "alln", "BY_2", 
                  "CA_sm", "CA_med", "CA_lar")
    source("~/DisMAP project/Location, Location, Location/Location Workshop/GAM_SpaceTime_Config3.R") 
    
    save(gam_Ran_Nte, gam_Ran_Pte, gam_Tar_N_1te, gam_Tar_P_1te, gam_Dist_N_npnte, gam_Dist_P_npnte, gam_Dist_N_npote,
         gam_Dist_P_npote, gam_Dist_N_mpnte, gam_Dist_P_mpnte, gam_Dist_N_mpote, gam_Dist_P_mpote, gam_Dist_N_spnte,
         gam_Dist_P_spnte, gam_Dist_N_spote, gam_Dist_P_spote, gam_Dist_N_allnte, gam_Dist_P_allnte, gam_Dist_N_allote,
         gam_Dist_P_allote, gam_BY_Nte, gam_BY_Pte, gam_CA_N_smte, gam_CA_P_smte, gam_CA_N_medte,
         gam_CA_P_medte, gam_CA_N_larte, gam_CA_P_larte, file="GAM_EST_Full_data_11_8_21_had.RData")
    
    #BRT-spacetime
    sampling <- c("ran", "tar_0.5", "npo", "npn", "mpo", 
                  "mpn", "spo", "spn","allo", "alln", "BY_2", 
                  "CA_sm", "CA_med", "CA_lar")
    source("~/DisMAP project/Location, Location, Location/Location Workshop/BRT_spacetime.R") 
   
    save(brt_R_N_te, brt_R_P_te, brt_T_N_1_te, brt_T_P_1_te, brt_dist_N_npn_te, brt_dist_P_npn_te, brt_dist_N_npo_te,
         brt_dist_P_npo_te, brt_dist_N_mpn_te, brt_dist_P_mpn_te, brt_dist_N_mpo_te, brt_dist_P_mpo_te, brt_dist_N_spn_te,
         brt_dist_P_spn_te, brt_dist_N_spo_te, brt_dist_P_spo_te, brt_dist_N_alln_te, brt_dist_P_alln_te, brt_dist_N_allo_te,
         brt_dist_P_allo_te, brt_B_N_te, brt_B_P_te, brt_CA_N_sm_te, brt_CA_P_sm_te, brt_CA_N_med_te, brt_CA_P_med_te,
         brt_CA_N_lar_te, brt_CA_P_lar_te, file="BRTs_EST_Full_data_11_8_21_had.RData")
    
  #### save the Rdata to load later ####
  saveRDS(dat_hist, "dat_hist_results_full_11_8_21_had_Stephs_EST.rds")  # * 'full' [full models] or 'temp' [temp-only models]
  saveRDS(dat_fcast, "dat_fcast_results_full_11_8_21_had_Stepsh_EST.rds")
  # all_mods <- c("gam_Ran", "gam_Tar_0.5", "gam_Dist_npo", "gam_Dist_npn",
  #               "gam_Dist_mpo", "gam_Dist_mpn", "gam_Dist_spo", 
  #               "gam_Dist_spn", "gam_Dist_allo", "gam_Dist_alln", 
  #               "gam_BY_2", "gam_CA_sm", "gam_CA_med", "gam_CA_lar",
  #               
  #               "brt_Ran","brt_Tar_0.5", "brt_Dist_npo","brt_Dist_npn",
  #               "brt_Dist_mpo","brt_Dist_mpn","brt_Dist_spo",
  #               "brt_Dist_spn", "brt_Dist_allo","brt_Dist_alln", 
  #               "brt_BY_2", "brt_CA_sm", "brt_CA_med", "brt_CA_lar")
  # 
  # save(list = ls(pattern = paste0(all_mods, collapse="|")),
  #      file = "saved_models_11_8_21_had_Stephs.RData")  # * 'full' or 'temp'; this will save all models of these names, even when outdated, if they're in the environment


  ##"gam_Ran_te","gam_Tar_0.5_te", "gam_Dist_npo_te", "gam_Dist_npn_te","gam_Dist_mpo_te", 
 ## "gam_Dist_mpn_te","gam_Dist_spo_te", "gam_Dist_spn_te", "gam_Dist_allo_te", 
  ##"gam_Dist_alln_te", "gam_BY_te", "gam_CA_sm_te", "gam_CA_med_te", "gam_CA_lar_te",
  
  ##"brt_Ran_te","brt_Tar_0.5_te", "brt_Dist_npo_te","brt_Dist_npn_te","brt_Dist_mpo_te","brt_Dist_mpn_te","brt_Dist_spo_te",
  ##"brt_Dist_spn_te","brt_Dist_allo_te","brt_Dist_alln_te", "brt_BY_te", "brt_CA_sm_te",
  ##"brt_CA_med_te", "brt_CA_lar_te"
