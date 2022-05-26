# Fishery-dependent-SDM-projections
Simulation looking at the impact of using fishery-dependent data in SDM and forecasting distributions under future climate change

**Code authors**: Melissa Karp, Stephanie Brodie, James Smith, Owen Liu, Kate Richerson, Becca Selden

**Relevant paper:**
Karp et al. Projecting species distributions using fishery-dependent data.   

**File Descriptions:**
1. SimulatedWorld_ROMS_FishDep_Final_updatedPres.R: function to simulate species presence and abundance in time and space with respect to environmental habitat preferences, and the different fishery-dependent sampling scenarios.

2. ModelComparison_FishSuitability_v_10_28_21.R: this code uses the function above to generate data, then builds GAM and BRT models, and makes predictions into the future 2011-2100

3. Fitting_BRTs.R: code to fit the boosted regression tree models to the simulated data with just environmental covariates. This is called in 2 above. 

4. Fitting_GAMs.R: code to fit the generalized additive models to the simulated data with just environmental covariates. This is called in 2 above. 

5. BRT_spacetime.R: code to fit the boosted regression tree models to the simulated data with environmental and space and time covariates. This is called in 2 above.

6. GAM_SpaceTime_Config3.R: code to fit the generalized additive models to the simulated data with environmental and space and time covariates. This is called in 2 above.

7. calculate_kl_ks_hd_and_cohens_d_All.R: code to calculate the climate bias and climate novelty (e.g., Hellinger Distance and Cohen's d) 

8. DistancetoPort.R: This code calculates the distance from every cell in the ROMS extent to 5 different fishing ports along the US West Coast in CA, OR, and WA. 

