library(tidyverse)
library(dplyr)
library(philentropy) #K-L distance (binning method)
library(statip) #Hellinger function 
#library(ggrepel) 
library(RColorBrewer)
library(gridExtra)
library(readxl)
install.packages("ggridges")
library(ggridges)
library(viridis)

#Load the data
dat_hist <- readRDS("~/DisMAP project/Location, Location, Location/Location Workshop/dat_hist_results_full_11_8_21_had_Stephs.rds")
dat_fcast <- readRDS("~/DisMAP project/Location, Location, Location/Location Workshop/dat_fcast_results_full_11_8_21_had_Stepsh.rds")

dat_hist <- dat_hist %>% 
  mutate(all_sampled = 1) #add a dummy column for ALL sampled if we want to use that 

dat_fcast <- dat_fcast %>% 
  mutate(all_sampled = 1) #add a dummy column for ALL sampled if we want to use that 

#This function is taken directly from here: https://stackoverflow.com/questions/15436702/estimate-cohens-d-for-effect-size
# Modified on 5April2021 to not take absolute value of mean difference
#It looks correct to me and results match cohen.d() from effsize package, but using this means you don't have to install an extra package
cohens_d <- function(x, y) {
  lx <- length(x)- 1
  ly <- length(y)- 1
  md  <- mean(x) - mean(y)        ## mean difference (numerator)
  csd <- lx * var(x) + ly * var(y)
  csd <- csd/(lx + ly)
  csd <- sqrt(csd)                     ## common sd computation
  
  cd  <- md/csd                        ## cohen's d
  
  return(cd)
}

#Slightly altered from Steven's code. My change is just setting n = 1024 instead of 1000 because documentation says "it almost always makes sense to specify n as a power of two"
norm_vec <- function(x) sqrt(sum(x^2))

hell_dist <- function (p, q, from, to, n = 1024) {
  P <- density(p, kernel = "gaussian", from = from, to = to, n = n)
  p <- P$y
  p <- p / sum(p)
  Q <- density(q, kernel = "gaussian", from = from, to = to, n = n)
  q <- Q$y
  q <- q / sum(q)
  hd <- norm_vec(sqrt(p) - sqrt(q)) / sqrt(2)
  hd
}

#Function that returns Cohen's d, K-L distance, and Hellinger distance for any given pair of comparisons
#Predictor is the name of the column with the variable we want to compare, sampregime is the sampling regime column name
#Returns a list with the data frame of results as the first element and a ggplot as the second element


compare_dat <- function(predictor, #Can be "temp" "zoo_200", "mld" or "chl_surface"
                        sampregime1, #Can be any of the sampling regimes (e.g."random_sampled", "pref_sampled_1") or "all" for all data
                        sampregime2,
                        data1, #Data frame with sampregime1 (e.g dat_hist, dat_fcast)
                        data2) #Data frame with sampregime2 
{
  
  #dat <- dat_hist
  
  #vectors of the data we are comparing (sampling regime 1 and 2)
  var1 <- data1[ ,predictor][data1[ ,sampregime1] == 1]
  var2 <- data2[ ,predictor][data2[ ,sampregime2] == 1]
  
  #Calculate Cohen's D
  cd <- cohens_d(var1,
                 var2)
  
  #Do we want a descriptor of the effect size for reference?
  #Note: Cohen suggested that d=0.2 be considered a 'small' effect size, 0.5 represents a 'medium' effect size and 0.8 a 'large' effect size
  cd_effect <- ifelse(cd <= 0.2, "small",
                      ifelse(cd > 0.2 & cd <= 0.5, "medium", 
                             "large"))
  
  #Calculate K-L distance with a binned approach using density
  
  n_bins <- 1024 #arbitrary
  
  bins <- seq(floor(min(c(var1, var2))), ceiling(max(c(var1, var2))), length = n_bins)
  
  var1_dens <- density(var1, n=n_bins, from=bins[1], to=bins[n_bins])
  var2_dens <- density(var2, n=n_bins, from=bins[1], to=bins[n_bins])
  
  # density output does not sum to 1, take sum of the density vector and scale so sums to 1
  dens1 <- tibble(x=var1_dens$x,
                  dens1=var1_dens$y) %>%
    mutate(total_dens1=sum(dens1),
           rel_dens1=dens1/total_dens1)
  
  dens2 <- tibble(x=var2_dens$x,
                  dens2=var2_dens$y) %>%
    mutate(total_dens2=sum(dens2),
           rel_dens2=dens2/total_dens2)
  
  
  dens_join <- full_join(dens1, dens2, by = "x") %>% 
    replace(., is.na(.), 0)
  
  kld_bin <- suppressMessages(philentropy::KL(rbind(dens_join$rel_dens1, dens_join$rel_dens2)) %>% as.numeric())
  
  #Hellinger distance (continuous version). Integration fails for comparing chl for all vs. dist_sampled_mpn, so I have it return an NA
  hell_dist_cont <- tryCatch(statip::hellinger(var1, var2, lower = -Inf, upper = Inf), error=function(err) NA) 
  
  #Discrete Hellinger distance, using same bins as K-L distance
  hell_dist_discr <- hell_dist(var1, var2, from=bins[1], to=bins[n_bins])
  
  #Output
  out <- data.frame("sampling_regime_1" = sampregime1,
                    "sampling_regime_2" = sampregime2,
                    "data_1" = deparse(substitute(data1)),
                    "data_2" = deparse(substitute(data2)),
                    "predictor" = predictor,
                    "cohens_d" = cd, 
                    "cohens_d_effect" = cd_effect, 
                    "kullback-leibler_dist" = kld_bin,
                    "hellinger_dist" = hell_dist_discr,
                    "hellinger_dist_discr" = hell_dist_discr)
  
  #plot
  df1 <- data.frame(var = var1, 
                    regime = sampregime1)
  
  df2 <- data.frame(var = var2, 
                    regime = sampregime2)
  
  df <- rbind(df1, df2)
  
  plotlabel <- paste0("Cohen's D = ", round(cd, 3), ", ", cd_effect, "\n",
                      "Kullback-Leibler distance = ", round(kld_bin, 3),"\n",
                      "Hellinger distance = ", round(hell_dist_discr, 3))
  
  histplot <- ggplot(df, aes(x = var, group = regime, fill = regime)) +
    #geom_histogram(position = "dodge", bins = 30)+
    stat_bin(aes(y = ..density..), position = 'dodge', bins = 30)+
    annotate(geom = "label", x = Inf, y = Inf, label = plotlabel,  vjust = "inward", hjust = "inward")+
    theme_bw()+
    ggtitle(paste0(predictor, ", ", deparse(substitute(data1)), " ", sampregime1, " ", " ", "vs.\n ", deparse(substitute(data2)), " ", sampregime2, " ")) +
    labs(x = predictor)
  
  return(list(out,histplot))
  
}

test<-compare_dat("temp",
                  "all_sampled",
                  "pref_sampled",
                  dat_hist,
                  dat_hist)


#Apply to all combinations of random + fishery dependent sampling, predictors

fishdep_regimes <- names(dat_hist)[grepl("sampled", names(dat_hist)) &
                                     #!grepl("random_", names(dat_hist)) &
                                     !grepl("all_", names(dat_hist))]

#Doing historic and future loops separately out of laziness. 

predictors <- c("temp", "zoo_200", "mld", "chl_surface")


hist_comps <- expand.grid(c("random_sampled", "all_sampled", "pres"), fishdep_regimes, predictors, stringsAsFactors = FALSE) %>%
  rename("regime1" = Var1, "regime2" = Var2, "predictor" = Var3)%>%
  filter(paste0(regime1, regime2) != "random_sampledrandom_sampled") %>% 
  #filter(paste0(regime1, regime2) != "presrandom_sampled") %>% 
  #bind_rows(tibble(regime1 = rep("all_sampled", length(predictors)), regime2 = rep("pres", length(predictors)), predictor = predictors)) %>% 
  arrange(regime1, regime2)

hist_df_res <- list()
hist_plot_res <- list()

for(i in 1:nrow(hist_comps))
{
  res <- compare_dat(predictor = hist_comps$predictor[i],
                     sampregime1 = hist_comps$regime1[i],
                     sampregime2 = hist_comps$regime2[i],
                     data1 = dat_hist,
                     data2 = dat_hist)
  
  hist_df_res[[i]] <- res[[1]]
  hist_plot_res[[i]] <- res[[2]]
  
}

#put results into a data frame
hist_res_df <- bind_rows(hist_df_res) %>%  #not sure why this changes kullback_leibler_dist to kullback.leibler_dist??
  rename(kullback_leibler_dist = kullback.leibler_dist) %>% 
  mutate(samp_grp = case_when(grepl("pref", sampling_regime_2) ~ "pref",
                              grepl("dist", sampling_regime_2) ~ "dist",
                              grepl("Closed", sampling_regime_2) ~ "closed",
                              grepl("BY", sampling_regime_2) ~ "bycatch",
                              grepl("random", sampling_regime_2) ~ "random")) %>% 
  arrange(samp_grp, sampling_regime_2) %>% 
  mutate(samp_reg = sampling_regime_2) %>% #Display names for plotting
  mutate(samp_reg = gsub("_sampled_", " ", samp_reg),
         samp_reg = gsub("BY", "BY", samp_reg),
         samp_reg = ifelse(samp_reg == "random_sampled", "Random", samp_reg),
         samp_reg = gsub("dist", "Dist", samp_reg),
         samp_reg = ifelse(samp_reg == "pref_sampled", "Pref", samp_reg)) %>% 
  rename(`Sampling regime` = samp_reg) %>% 
  mutate(disp_predictor = case_when(predictor == "temp" ~ "Temperature",
                                    predictor == "zoo_200" ~"Zooplankton",
                                    predictor == "mld" ~ "Mixed layer depth",
                                    predictor == "chl_surface" ~ "Chlorophyll")) %>% 
  # left_join(future_rmse %>% rename(sampling_regime_2 = `Sampling Regime`,
  #                                  `Mean RMSE` = RMSE_hist), by = "sampling_regime_2") %>% 
  mutate(Comparison = case_when(sampling_regime_1 == "all_sampled" ~ "All",
                                sampling_regime_1 == "random_sampled" ~ "Random",
                                sampling_regime_1 == "pres" ~ "Presence"))

# #generate colors for different groups
# pref_pal <- colorRampPalette(brewer.pal(9, "PuBu"))
# dist_pal <- colorRampPalette(brewer.pal(9, "YlOrRd"))
# closed_pal <- colorRampPalette(brewer.pal(9, "Greens"))
# 
# #Make a darker version
# pref_pal2 <- colorRampPalette(brewer.pal(9, "PuBu")[4:8])
# dist_pal2 <- colorRampPalette(brewer.pal(9, "YlOrRd")[4:9])
# closed_pal2 <- colorRampPalette(brewer.pal(9, "Greens")[4:8])
# 
# #do we want the gray  "vs all" in there now? I think so...
# #Plot all comparisons on same axes -- looks a little too busy
# ggplot(hist_res_df, 
#        aes(x = cohens_d, y = hellinger_dist_discr)) +
#   geom_vline(xintercept = 0, color = "gray70") + #for referenc
#   geom_point(aes(color = `Sampling regime`, shape = Comparison, size = `Mean RMSE`)) +
#   scale_color_manual(values = c("black", closed_pal2(3), dist_pal2(8), pref_pal2(5),"gray70"))+
#   #scale_shape_manual(values = c(18, rep(17, 3), rep(15, 8), rep(19, 5), 4))+
#   facet_wrap(~ disp_predictor) +
#   theme_bw()+
#   xlab("Cohen's D")+
#   ylab("Hellinger distance")+
#   ggtitle("Fishery dependent data compared to\nall, randomly-sampled, and species-presence data")+
#   NULL
# 
# #Try facet_grid -- looks better
# hist_bias_grid_plot <- ggplot(hist_res_df, 
#                               aes(x = cohens_d, y = hellinger_dist_discr)) +
#   geom_vline(xintercept = 0, color = "gray70") + #for referenc
#   geom_point(aes(color = `Sampling regime`, size = `Mean RMSE`)) +
#   scale_color_manual(values = c("black", closed_pal2(3), dist_pal2(8), pref_pal2(5),"gray70"))+
#   #scale_shape_manual(values = c(18, rep(17, 3), rep(15, 8), rep(19, 5), 4))+
#   facet_grid(Comparison ~ disp_predictor) +
#   theme_bw()+
#   xlab("Cohen's D")+
#   ylab("Hellinger distance")+
#   ggtitle("Fishery dependent data compared to\nall, randomly-sampled, and species-presence data")+
#   NULL
# 
# ggsave(hist_bias_grid_plot, filename = "hist_bias_grid_plot.png")


################Compare fishery dependent sampling to future all conditions

#First, put the data into time bins:
dat_fcast_bins <- dat_fcast %>% 
  mutate(per_2011_2039 = ifelse(year >= 2011 & year <= 2039, 1, 0)) %>% 
  mutate(per_2040_2069 = ifelse(year >= 2040 & year <= 2069, 1, 0)) %>%
  mutate(per_2070_2100 = ifelse(year >= 2070 & year <= 2100, 1, 0)) 

future_comps <- expand.grid(c("per_2011_2039", "per_2040_2069", "per_2070_2100"), fishdep_regimes, predictors, stringsAsFactors = FALSE) %>% 
  rename("regime1" = Var1, "regime2" = Var2, "predictor" = Var3) %>% 
  arrange(regime1, regime2)

future_df_res <- list()
future_plot_res <- list()

for(i in 1:nrow(future_comps))
{
  res <- compare_dat(predictor = future_comps$predictor[i],
                     sampregime1 = future_comps$regime1[i],
                     sampregime2 = future_comps$regime2[i],
                     data1 = dat_fcast_bins,
                     data2 = dat_hist)
  
  future_df_res[[i]] <- res[[1]]
  future_plot_res[[i]] <- res[[2]]
  
}

#put results into a data frame
future_res_df <- bind_rows(future_df_res) %>%  #not sure why this changes kullback_leibler_dist to kullback.leibler_dist??
  rename(kullback_leibler_dist = kullback.leibler_dist) %>% 
  mutate(samp_grp = case_when(grepl("pref", sampling_regime_2) ~ "pref",
                              grepl("dist", sampling_regime_2) ~ "dist",
                              grepl("Closed", sampling_regime_2) ~ "closed",
                              grepl("BY", sampling_regime_2) ~ "bycatch",
                              grepl("random", sampling_regime_2) ~ "random")) %>% 
  arrange(samp_grp, sampling_regime_2) %>% 
  mutate(samp_reg = sampling_regime_2) %>% #Display names for plotting
  mutate(samp_reg = gsub("_sampled_", " ", samp_reg),
         samp_reg = gsub("BY", "BY", samp_reg),
         samp_reg = ifelse(samp_reg == "random_sampled", "Random", samp_reg),
         samp_reg = gsub("dist", "Dist", samp_reg),
         samp_reg = ifelse(samp_reg == "pref_sampled", "Pref", samp_reg)) %>% 
  rename(`Sampling regime` = samp_reg) %>% 
  mutate(disp_predictor = case_when(predictor == "temp" ~ "Temperature",
                                    predictor == "zoo_200" ~"Zooplankton",
                                    predictor == "mld" ~ "Mixed layer depth",
                                    predictor == "chl_surface" ~ "Chlorophyll")) %>% 
  mutate(Comparison = case_when(sampling_regime_1 == "per_2011_2039" ~ "2011-2039",
                                sampling_regime_1 == "per_2040_2069" ~ "2040-2069",
                                sampling_regime_1 == "per_2070_2100" ~ "2070-2100"))


###STOP HERE AND GO TO RMSE_vs_HD_fishdep.R to PLOT FACET GRAPHS!!! ####

#############################
#Bias and RMSE scatterplots
#############################

#Plot all on 1 panel
future_bias_all_plot <- ggplot(future_res_df, 
                               aes(x = cohens_d, y = hellinger_dist_discr)) +
  geom_vline(xintercept = 0, color = "gray70") + #for referenc
  geom_point(aes(color = `Sampling regime`, shape = Comparison)) +
  scale_color_manual(values = c("black", closed_pal2(3), dist_pal2(8), pref_pal2(5),"gray70"))+
  #scale_shape_manual(values = c(18, rep(17, 3), rep(15, 8), rep(19, 5), 4))+
  facet_wrap(~ disp_predictor) +
  theme_bw()+
  xlab("Cohen's D")+
  ylab("Hellinger distance")+
  ggtitle("Fishery dependent data compared to all future data")+
  NULL

ggsave(future_bias_all_plot, filename = "future_bias_all_plot.png")

#Plot different time periods separately

#2011-2039
future_bias_2011_2039_plot <- ggplot(future_res_df %>% filter(Comparison == "2011-2039"), 
                                     aes(x = cohens_d, y = hellinger_dist_discr)) +
  geom_vline(xintercept = 0, color = "gray70") + #for referenc
  geom_point(aes(color = `Sampling regime`, size = `Mean RMSE`)) +
  scale_size_continuous(limits = c(min(future_res_df$`Mean RMSE`), max(future_res_df$`Mean RMSE`)))+
  scale_color_manual(values = c("black", closed_pal2(3), dist_pal2(8), pref_pal2(5),"gray70"))+
  facet_wrap(~ disp_predictor) +
  theme_bw()+
  xlab("Cohen's D")+
  ylab("Hellinger distance")+
  ggtitle("Fishery dependent data compared to 2011-2039 data")+
  xlim(min(future_res_df$cohens_d), max(future_res_df$cohens_d))+
  ylim(min(future_res_df$hellinger_dist_discr), max(future_res_df$hellinger_dist_discr))+
  NULL

ggsave(future_bias_2011_2039_plot, filename = "future_bias_2011_2039_plot.png")

#2040-2069
future_bias_2040_2069_plot <- ggplot(future_res_df %>% filter(Comparison == "2040-2069"), 
                                     aes(x = cohens_d, y = hellinger_dist_discr)) +
  geom_vline(xintercept = 0, color = "gray70") + #for referenc
  geom_point(aes(color = `Sampling regime`, size = `Mean RMSE`)) +
  scale_size_continuous(limits = c(min(future_res_df$`Mean RMSE`), max(future_res_df$`Mean RMSE`)))+
  scale_color_manual(values = c("black", closed_pal2(3), dist_pal2(8), pref_pal2(5),"gray70"))+
  facet_wrap(~ disp_predictor) +
  theme_bw()+
  xlab("Cohen's D")+
  ylab("Hellinger distance")+
  ggtitle("Fishery dependent data compared to 2040-2069 data")+
  xlim(min(future_res_df$cohens_d), max(future_res_df$cohens_d))+
  ylim(min(future_res_df$hellinger_dist_discr), max(future_res_df$hellinger_dist_discr))+
  NULL

ggsave(future_bias_2040_2069_plot, filename = "future_bias_2040_2069_plot.png")

#2070-2100
future_bias_2070_2100_plot <- ggplot(future_res_df %>% filter(Comparison == "2070-2100"), 
                                     aes(x = cohens_d, y = hellinger_dist_discr)) +
  geom_vline(xintercept = 0, color = "gray70") + #for referenc
  geom_point(aes(color = `Sampling regime`, size = `Mean RMSE`)) +
  scale_size_continuous(limits = c(min(future_res_df$`Mean RMSE`), max(future_res_df$`Mean RMSE`)))+
  scale_color_manual(values = c("black", closed_pal2(3), dist_pal2(8), pref_pal2(5),"gray70"))+
  facet_wrap(~ disp_predictor) +
  theme_bw()+
  xlab("Cohen's D")+
  ylab("Hellinger distance")+
  ggtitle("Fishery dependent data compared to 2070-2100 data")+
  xlim(min(future_res_df$cohens_d), max(future_res_df$cohens_d))+
  ylim(min(future_res_df$hellinger_dist_discr), max(future_res_df$hellinger_dist_discr))+
  NULL

ggsave(future_bias_2070_2100_plot, filename = "future_bias_2070_2100_plot.png")

#Can also try one big grid?
future_bias_grid_plot <- ggplot(future_res_df, 
                                aes(x = cohens_d, y = hellinger_dist_discr)) +
  geom_vline(xintercept = 0, color = "gray70") + #for referenc
  geom_point(aes(color = `Sampling regime`, size = `Mean RMSE`)) +
  scale_color_manual(values = c("black", closed_pal2(3), dist_pal2(8), pref_pal2(5),"gray70"))+
  facet_grid(Comparison ~ disp_predictor) +
  theme_bw()+
  xlab("Cohen's D")+
  ylab("Hellinger distance")+
  ggtitle("Fishery dependent data compared to future data")+
  NULL

ggsave(future_bias_grid_plot, filename = "future_bias_grid_plot.png")


##########################################
#Ridge plots to compare data distributions
#########################################

#Put historic data in long format
dat_hist_long <- dplyr::select(dat_hist, all_sampled, pres, random_sampled, pref_sampled, dist_sampled_npo, dist_sampled_npn, dist_sampled_mpo, dist_sampled_mpn, dist_sampled_spo, dist_sampled_spn, dist_sampled_allo, dist_sampled_alln, BY_sampled_2, Closed_sampled_1, Closed_sampled_2, Closed_sampled_3, temp, zoo_200,  mld, chl_surface) %>% 
  gather(reg, sampled, all_sampled: Closed_sampled_3) %>% 
  filter(sampled == 1) %>%
  mutate(samp_reg = reg)%>%  #Display names for plotting
  mutate(samp_reg = gsub("_sampled_", " ", samp_reg),
         samp_reg = ifelse(samp_reg == "BY_sampled_2", "BY 2", samp_reg),
         samp_reg = ifelse(samp_reg == "all_sampled", "All", samp_reg),
         samp_reg = gsub("dist", "Dist", samp_reg),
         samp_reg = ifelse(samp_reg=="random_sampled", "Random", samp_reg),
         samp_reg = ifelse(samp_reg=="pref_sampled", "Pref", samp_reg)) %>%
  # arrange(desc(samp_reg)) %>%
  rename(`Sampling regime` = samp_reg)

#Just pull out bias results comparing All to other sampling regimes
hist_res_df2 <- filter(hist_res_df, sampling_regime_1 == "all_sampled") %>% 
  rename(reg = sampling_regime_2)

future_res_df2 <- future_res_df %>% 
  rename(reg = sampling_regime_2)

#get bias results to add in 
# dat_hist_long_temp <- dat_hist_long %>% 
#   left_join(dplyr::select(filter(hist_res_df2, predictor == "temp"), `Sampling regime`, cohens_d, hellinger_dist_discr), by = "Sampling regime") %>% 
#   filter(reg !="pres") %>% 
#   filter(`Sampling regime` != "random_sampled") %>% 
#   rename(`Cohen's D` = cohens_d,
#          `Hellinger distance` = hellinger_dist_discr)

dat_hist_bias <- dat_hist_long %>% 
  gather(predictor, value, temp:chl_surface) %>% 
  left_join(dplyr::select(hist_res_df2, `Sampling regime`, cohens_d, hellinger_dist_discr, predictor), by = c("Sampling regime","predictor")) %>% 
  filter(reg !="pres") %>% 
  #filter(reg != "random_sampled") %>% 
  rename(`Cohen's D` = cohens_d,
         `Hellinger distance` = hellinger_dist_discr) %>% 
  mutate(Predictor = case_when (predictor == "chl_surface" ~ "Chlorophyll",
                                predictor == "mld" ~ "Mixed layer depth",
                                predictor == "temp" ~ "Temperature",
                                predictor == "zoo_200" ~ "Zooplankton"))

#Get the future data that we'll also be comparing to
dat_fcast_long <- dat_fcast_bins %>% 
  dplyr::select(per_2011_2039, per_2040_2069, per_2070_2100, temp, zoo_200,  mld, chl_surface) %>% 
  gather(reg, sampled, per_2011_2039: per_2070_2100) %>% 
  filter(sampled == 1) %>% 
  mutate(`Sampling regime` = case_when(reg == "per_2011_2039" ~ "All 2011-2039",
                                       reg == "per_2040_2069" ~ "All 2040-2069",
                                       reg == "per_2070_2100" ~ "All 2070-2100"))

#We also want to show future all vs (historic) fishery dependent. Combine future all data with historic fishdep, then add in bias
dat_future_bias <- dat_fcast_long %>% 
  gather(predictor, value, temp:chl_surface) %>% 
  bind_rows(dat_hist_long %>% 
              gather(predictor, value, temp:chl_surface)) %>% 
  filter(!`Sampling regime` %in% c("pres", "All")) %>%  #I don't THINK we want to plot these....
  left_join(dplyr::select(future_res_df2, `Sampling regime`, cohens_d, hellinger_dist_discr, predictor, comparison_period = sampling_regime_1), by = c("Sampling regime","predictor")) %>% 
  rename(`Cohen's D` = cohens_d,
         `Hellinger distance` = hellinger_dist_discr) %>% 
  mutate(Predictor = case_when (predictor == "chl_surface" ~ "Chlorophyll",
                                predictor == "mld" ~ "Mixed layer depth",
                                predictor == "temp" ~ "Temperature",
                                predictor == "zoo_200" ~ "Zooplankton"))

#save the max and min values of the bias metrics to keep a common scale across plots?
#Not sure whether this should be common across all predictors...
mincd <- min(c(hist_res_df$cohens_d, future_res_df$cohens_d))
maxcd <- max(c(hist_res_df$cohens_d, future_res_df$cohens_d))
minhd <- min(c(hist_res_df$hellinger_dist_discr, future_res_df$hellinger_dist_discr))
maxhd <- max(c(hist_res_df$hellinger_dist_discr, future_res_df$hellinger_dist_discr))

#Create historic temperature ridgeplots####################

hist_temp_p1<- ggplot(dat_hist_bias %>% filter(predictor == "temp"), aes(x = value, y = `Sampling regime`, fill = `Cohen's D`))+
  geom_density_ridges_gradient(scale = 4) +
  theme_bw()+
  theme(legend.position="top") +
  scale_fill_gradient2(limits = c(mincd, maxcd))+ 
  xlab("Temperature")


hist_temp_p2<-ggplot(dat_hist_bias %>% filter(predictor == "temp"), aes(x = value, y = `Sampling regime`, fill = `Hellinger distance`))+
  geom_density_ridges_gradient(scale = 4) + 
  theme_bw()+
  theme(legend.position="top") +
  scale_fill_viridis(limits = c(minhd, maxhd))+
  xlab("Temperature")

hist_temp_plot <- grid.arrange(hist_temp_p1, hist_temp_p2, nrow = 1)

ggsave(hist_temp_plot, filename = "historic_temp_ridgeplot.png")

#Create future temperature ridgeplots#####################
future_2011_2039_temp_p1<- ggplot(dat_future_bias %>% filter(predictor == "temp" & 
                                                               (comparison_period == "per_2011_2039" | `Sampling regime` == "All 2011-2039")), aes(x = value, y = `Sampling regime`, fill = `Cohen's D`))+
  geom_density_ridges_gradient(scale = 4) +
  theme_bw()+
  theme(legend.position="top") +
  scale_fill_gradient2(limits = c(mincd, maxcd))+ #
  xlab("Temperature")

future_2011_2039_temp_p2<- ggplot(dat_future_bias %>% filter(predictor == "temp" & 
                                                               (comparison_period == "per_2011_2039" | `Sampling regime` == "All 2011-2039")), aes(x = value, y = `Sampling regime`, fill = `Hellinger distance`))+
  geom_density_ridges_gradient(scale = 4) +
  theme_bw()+
  theme(legend.position="top") +
  scale_fill_viridis(limits = c(minhd, maxhd))+ # 
  xlab("Temperature")

future1_temp_plot <- grid.arrange(future_2011_2039_temp_p1, future_2011_2039_temp_p2, nrow = 1)

ggsave(future1_temp_plot, filename = "future_2011_2039_temp_ridgeplot.png")

####Let's try all predictors in one plot
###Historic
hist_all_p1<- ggplot(dat_hist_bias, aes(x = value, y = `Sampling regime`, fill = `Cohen's D`))+
  geom_density_ridges_gradient(scale = 4) +
  theme_bw()+
  theme(legend.position="top", legend.title=element_text(size=10)) +
  theme(axis.text = element_text(size =6, face="bold"))+
  scale_fill_gradient2(limits = c(mincd, maxcd))+ 
  facet_wrap(~Predictor, scales = "free_x", ncol = 1)+
  xlab("Value")

hist_all_p2<- ggplot(dat_hist_bias, aes(x = value, y = `Sampling regime`, fill = `Hellinger distance`))+
  geom_density_ridges_gradient(scale = 4) +
  theme_bw()+
  theme(legend.position="top", legend.title=element_text(size=10)) +
  theme(axis.text = element_text(size =6, face="bold"))+
  scale_fill_viridis(limits = c(minhd, maxhd))+ 
  facet_wrap(~Predictor, scales = "free_x", ncol = 1)+
  xlab("Value")

hist_all_plot <- grid.arrange(hist_all_p1, hist_all_p2, nrow = 1)

ggsave(hist_all_plot, filename = "historic_all_predictors_ridgeplot.png", height = 12, width =8)

###Future "per_2011_2039" 
future_2011_2039_all_p1 <- ggplot(dat_future_bias %>% filter(comparison_period == "per_2011_2039" | `Sampling regime` == "All 2011-2039"), aes(x = value, y = `Sampling regime`, fill = `Cohen's D`))+
  geom_density_ridges_gradient(scale = 4) +
  theme_bw()+
  theme(legend.position="top", legend.title=element_text(size=10)) +
  theme(axis.text = element_text(size =6, face="bold"))+
  scale_fill_gradient2(limits = c(mincd, maxcd))+ 
  facet_wrap(~Predictor, scales = "free_x", ncol = 1)+
  xlab("Value")

future_2011_2039_all_p2 <- ggplot(dat_future_bias %>% filter(comparison_period == "per_2011_2039" | `Sampling regime` == "All 2011-2039"), aes(x = value, y = `Sampling regime`, fill = `Hellinger distance`))+
  geom_density_ridges_gradient(scale = 4) +
  theme_bw()+
  theme(legend.position="top", legend.title=element_text(size=10)) +
  theme(axis.text = element_text(size =6, face="bold"))+
  scale_fill_viridis(limits = c(minhd, maxhd))+ 
  facet_wrap(~Predictor, scales = "free_x", ncol = 1)+
  xlab("Value")

future_2011_2039_all_plot <- grid.arrange(future_2011_2039_all_p1, future_2011_2039_all_p2, nrow = 1)

ggsave(future_2011_2039_all_plot, filename = "future_2011_2039_all_predictors_ridgeplot.png", height = 12, width = 8)

###Future "per_2040_2069" 
future_2040_2069_all_p1 <- ggplot(dat_future_bias %>% filter(comparison_period == "per_2040_2069" | `Sampling regime` == "All 2040-2069"), aes(x = value, y = `Sampling regime`, fill = `Cohen's D`))+
  geom_density_ridges_gradient(scale = 4) +
  theme_bw()+
  theme(legend.position="top", legend.title=element_text(size=10)) +
  theme(axis.text = element_text(size =6, face="bold"))+
  scale_fill_gradient2(limits = c(mincd, maxcd))+ 
  facet_wrap(~Predictor, scales = "free_x", ncol = 1)+
  xlab("Value")

future_2040_2069_all_p2 <- ggplot(dat_future_bias %>% filter(comparison_period == "per_2040_2069" | `Sampling regime` == "All 2040-2069"), aes(x = value, y = `Sampling regime`, fill = `Hellinger distance`))+
  geom_density_ridges_gradient(scale = 4) +
  theme_bw()+
  theme(legend.position="top", legend.title=element_text(size=10)) +
  theme(axis.text = element_text(size =6, face="bold"))+
  scale_fill_viridis(limits = c(minhd, maxhd))+ 
  facet_wrap(~Predictor, scales = "free_x", ncol = 1)+
  xlab("Value")

future_2040_2069_all_plot <- grid.arrange(future_2040_2069_all_p1, future_2040_2069_all_p2, nrow = 1)

ggsave(future_2040_2069_all_plot, filename = "future_2040_2069_all_predictors_ridgeplot.png", height = 12, width = 8)

###Future "per_2070_2100"
future_2070_2100_all_p1 <- ggplot(dat_future_bias %>% filter(comparison_period == "per_2070_2100" | `Sampling regime` == "All 2070-2100"), aes(x = value, y = `Sampling regime`, fill = `Cohen's D`))+
  geom_density_ridges_gradient(scale = 4) +
  theme_bw()+
  theme(legend.position="top", legend.title=element_text(size=10)) +
  theme(axis.text = element_text(size =6, face="bold"))+
  scale_fill_gradient2(limits = c(mincd, maxcd))+ 
  facet_wrap(~Predictor, scales = "free_x", ncol = 1)+
  xlab("Value")

future_2070_2100_all_p2 <- ggplot(dat_future_bias %>% filter(comparison_period == "per_2070_2100" | `Sampling regime` == "All 2070-2100"), aes(x = value, y = `Sampling regime`, fill = `Hellinger distance`))+
  geom_density_ridges_gradient(scale = 4) +
  theme_bw()+
  theme(legend.position="top", legend.title=element_text(size=10)) +
  theme(axis.text = element_text(size =6, face="bold"))+
  scale_fill_viridis(limits = c(minhd, maxhd))+ 
  facet_wrap(~Predictor, scales = "free_x", ncol = 1)+
  xlab("Value")

future_2070_2100_all_plot <- grid.arrange(future_2070_2100_all_p1, future_2070_2100_all_p2, nrow = 1)

ggsave(future_2070_2100_all_plot, filename = "future_2070_2100_all_predictors_ridgeplot.png", height = 12, width = 8)
