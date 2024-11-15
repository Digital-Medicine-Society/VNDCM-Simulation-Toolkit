#### Load libraries ####
library(SimDesign)
library(mirt)
library(lavaan)
library(tidySEM)
library(missMethods)
library(parallel)
library(rsimsum)
#### ####

#### FINAL VERSION USED FOR MANUSCRIPT. #### 
#Attempt 5 - Including a new MLR(DHT_bar,wCRO + wPRO).
Design_5 <- createDesign(N = c(35,100,200,1000), meas_error_mag = c(0.5,1,1.5,2),
                         n_assess = c(1,3,5,7),missing_method = c('None', 'MNAR','MCAR'), missing_rate = c(0,0.1,0.25,0.4),
                         subset = !((missing_method == 'None' & missing_rate != 0) | (missing_method != 'None' & missing_rate == 0)))

Generate_with_latent_vals <- function(condition,fixed_objects=NULL) {
  
  #Use Attach if it's quicker
  Attach(condition)
  #### DHT_Data ####
  
  fluct_sd = 0.3
  
  meth_filt_sd = 0.5
  base_rate = 10000
  latent_effect = 1250
  
  thetas<- Generate_Daily_Thetas(N,fluct_sd)
  #with(condition,Generate_Daily_Thetas(N,fluct_sd))
  
  DHT_raw <- Generate_DHT_Data_Full_with_filtered_thetas(thetas,N,fluct_sd,meth_filt_sd,base_rate,latent_effect,meas_error_mag)
  #with(condition,Generate_DHT_Data_Full(thetas,N,fluct_sd,meth_filt_sd,base_rate,latent_effect,meas_error_mag))
  
  #apply missingness condition
  if (missing_method == "MCAR") {
    #(condition$missing_method == "MCAR") {
    
    DHT_data <- Apply_MCAR_Missing_Data_with_latents(DHT_raw,missing_rate)
    #DHT_data<- sapply(DHT_raw[8:15],add_missing,rate=missing_rate)
    #sapply(DHT_raw,add_missing,rate=condition$missing_rate)
    
  } else if (missing_method == "MNAR") {
    #(condition$missing_method == "MNAR") {
    
    DHT_data <- data.frame(DHT_raw[1:7],
                           delete_MNAR_1_to_x(DHT_raw[8:14],
                                              missing_rate,
                                              #condition$missing_rate,
                                              c("DHT_day_1", "DHT_day_2", "DHT_day_3",
                                                "DHT_day_4", "DHT_day_5", "DHT_day_6", "DHT_day_7"),
                                              x=4)
    )
    
    
  } else {
    
    DHT_data <- DHT_raw
    
  }
  
  #Vector of days to include in assessment (DHT-side only; COAs always use 7 days)
  days_to_include <- c((8-n_assess):7)
  #with(condition,c((8-n_assess):7))
  
  # Calc the DHT_mean data based on which assessment days are being included in the data collection ####
  if (n_assess != 1) {
    
    DHT_Means <- rowMeans(DHT_data[,7+days_to_include],na.rm = TRUE)
  } else {
    DHT_Means <- DHT_data[,14]
  }
  
  # Calc the method-filtered latent traits means based on which assessment days
  #are being included in the data collection
  if (n_assess != 1) {
    
    filtered_DHT_latent_Means <- rowMeans(DHT_data[,days_to_include],na.rm = TRUE)
  } else {
    filtered_DHT_latent_Means <- DHT_data[,7]
  }
  
  
  #### ####
  
  #Generate the mean latent trait over the assessment period for the Weekly PROS.
  #We're using 7 previous days in every case, no matter the number of DHT assessment days included,
  #as this is the typical period for COA's.
  days_to_include <- c(1:7)
  thetas_bar <- Generate_avg_thetas(thetas,days_to_include)
  #hist(thetas_bar)
  
  #### Main PRO Data ####
  
  #Set perception filter
  per_filt_sd = 1       #0.75 #0.5       #Varying this affects the strength of the day-to-day correlations
  # THE PERCEPTION FILTER IS NOT ASSUMED TO RANDOMLY CHANGE FROM DAY TO DAY.
  # THIS IS A 'FIXED' FILTER.
  
  #Set the b2 parameters (the ‘middle’ location parameters)
  b2  <- c(-1.0, -0.8, -0.8, -0.4, -0.4, 0, 0, 0.4, 0.4, 0.8, 0.8, 1.0)
  #Set reliability
  reliability<-0.8
  
  Weekly_PRO <-Simulate_4_12_IRT_data_return_latents(thetas_bar,per_filt_sd,N,b2,reliability)
  
  ####          ####
  
  #### Weekly ClinRO Data ####
  
  per_Filt_sd = 1 #0.5#0.25
  #b0 = c(-1.0, -0.6, -0.2, 0, 0.2, 0.6, 1.0)
  b0 = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75) #Brought the tails in
  # and made the transitions uniform to give n=35 a shot of converging with CFA
  reliability = 0.7     #Need a more nuanced way of using this to determine the alpha parameters
  #Currently coming out at 0.9 - fixed
  
  ClinRO_Data <- Simulate_5_7_ClinRo_IRT_Data_return_latents(thetas_bar,per_Filt_sd,N,b0,reliability)
  #ClinRO_Data <- Simulate_5_7_ClinRo_IRT_Data_Full(thetas_bar,per_Filt_sd,N,b0,reliability)
  #with(condition, Simulate_5_7_ClinRo_IRT_Data_Full(thetas_bar,per_Filt_sd,N,b0,reliability))
  
  ####               ####
  
  #### Daily Single-item PRO (reliability = 0.6??? No, alpha makes no sense for single item)####
  per_filt_sd = 0.5 #0.75
  
  # Create iMIC distribution and other thresholds
  imic <- rnorm(N, 0.5, 0.075) #with(condition,rnorm(N, 0.5, 0.075))
  thrd1 <- rnorm(N, -1.5, 0.075)#with(condition,rnorm(N, -1.5, 0.075))
  thrd2 <- rnorm(N, -0.5, 0.075)#with(condition,rnorm(N, -0.5, 0.075))
  thrd3 <- imic
  thrd4 <- rnorm(N, 1.5, 0.075)#with(condition,rnorm(N, 1.5, 0.075))
  
  thrds<- data.frame(thrd1,thrd2,thrd3,thrd4)
  #as.data.frame(cbind(thrd1,thrd2,thrd3,thrd4))
  
  #Sim the data
  Daily_PRO <- Sim_Daily_Single_Item_return_latents(thetas,N,thrds,per_filt_sd)
  #Sim_Daily_Single_Item(thetas,N,thrds,per_filt_sd)
  #with(condition,Sim_Daily_Single_Item(thetas,N,thrds,per_filt_sd))
  
  
  ##########################
  #REVISIT THIS TOMORROW - vary both of included days in mean and CFA model, or neither
  ##########################
  #Vector of days to include in assessment (DHT-side only; COAs always use 7 days)
  days_to_include <- c(1:7)
  #c((8-n_assess):7)
  #with(condition,c((8-n_assess):7))
  
  
  #Calc the Daily PRO mean data based on which assessment days are being included in the data collection
  if (n_assess != 1) {
    #(condition$n_assess != 1) {
    Daily_PRO_means <- rowMeans(Daily_PRO[,7+days_to_include])
  } else {
    Daily_PRO_means <- Daily_PRO[,14]
  }
  #For comparison
  #Full_Week_dPRO_means<- rowMeans(Daily_PRO[,c(1:7)])
  ####              ####
  
  # Calc the perception-filtered latent traits means based on which assessment days
  #are being included in the data collection
  if (n_assess != 1) {
    #(condition$n_assess != 1) {
    filtered_dPRO_latent_Means <- rowMeans(Daily_PRO[,days_to_include],na.rm = TRUE)
  } else {
    filtered_dPRO_latent_Means <- Daily_PRO[,7]
  }
  
  
  #Combine all generated data into data frame, rename columns and return
  dat<-data.frame(thetas, DHT_data,filtered_DHT_latent_Means,DHT_Means,
                  thetas_bar,Weekly_PRO,ClinRO_Data,Daily_PRO,filtered_dPRO_latent_Means,
                  Daily_PRO_means,subject=c(1:N))
  # dat<-data.frame(thetas, DHT_data,DHT_Means,
  #                 Weekly_PRO,ClinRO_Data,Daily_PRO,Daily_PRO_means)
  
  names(dat)[1:7]<-c("theta_day_1","theta_day_2","theta_day_3","theta_day_4",
                     "theta_day_5","theta_day_6","theta_day_7")
  names(dat)[8:14]<-c("filtered_theta_day_1","filtered_theta_day_2","filtered_theta_day_3",
                      "filtered_theta_day_4", "filtered_theta_day_5","filtered_theta_day_6",
                      "filtered_theta_day_7")
  
  return(dat)
}

Analyse_with_latents_additional_MLR <- function(condition, dat, fixed_objects = NULL) {
  
  #Pearson correlation with measured variables
  Avg <- cor(dat$DHT_Means,dat$Weekly_PRO,use="na.or.complete")
  
  #Pearson correlation with filtered latent variables
  Avg_filtered_latents <- with(dat,cor(filtered_DHT_latent_Means,filtered_latents,use="na.or.complete"))
  
  #### single regressions ####
  
  #DHT average data against the weekly_PRO#
  wPRO_regr <- lm(DHT_Means~Weekly_PRO,data=dat)
  wPRO_regr_R2 <- summary(wPRO_regr)$r.squared
  #as.numeric(summary(PRO_regr)$r.squared)
  
  wPRO_latent_regr<- lm(filtered_DHT_latent_Means~filtered_latents,data=dat)
  wPRO_latent_regr_R2<- summary(wPRO_latent_regr)$r.squared
  
  #### Multiple regressions ####
  
  #Use mean of Daily PRO as a single variable
  regr_all <- lm(DHT_Means~Weekly_PRO+Weekly_ClinRO+Daily_PRO_means,data=dat)
  
  #filtered latent traits, use mean of dPRO filtered latent as a single variable
  regr_all_filtered_latents <- lm(filtered_DHT_latent_Means~filtered_latents+filtered_latents.1
                                  +filtered_dPRO_latent_Means,data=dat)
  
  #Use each day of daily PRO as a separate variable
  #Revisit tomorrow - Should all 7 dPRO days be used when n_assess <7?
  regr_all_daybyday <- lm(DHT_Means~Weekly_PRO+Weekly_ClinRO
                          +PRO_day_1+PRO_day_2+PRO_day_3+PRO_day_4+PRO_day_5+PRO_day_6+PRO_day_7,data=dat)
  
  #filtered latent traits, use each day of daily PRO filtered latent as a separate variable
  regr_all_filtered_latents_daybyday <- lm(filtered_DHT_latent_Means~filtered_latents+filtered_latents.1
                                           +day_1.2+day_2.2+day_3.2+day_4.2+day_5.2+day_6.2+day_7.2,data=dat)
  
  #Include just the weekly COAs
  MLR_weekly_COAs <- lm(DHT_Means~Weekly_PRO+Weekly_ClinRO,data=dat)
  
  #filtered latent traits, include just the weekly COAs
  MLR_weekly_COAs_filtered_latents <- lm(filtered_DHT_latent_Means~filtered_latents+filtered_latents.1,data=dat)
  
  
  
  ####Calculate R-squareds ####
  #Daily PRO as a single variable
  #Measured values
  all_adjR2 <- summary(regr_all)$adj.r.squared
 
  #filtered latents
  all_filtered_latents_adjR2 <- summary(regr_all_filtered_latents)$adj.r.squared
  
  #Each day of daily PRO as a separate variable
  all_daybyday_adjR2 <- summary(regr_all_daybyday)$adj.r.squared
  
  #filtered latents
  all_filtered_latents_daybyday_adjR2 <- summary(regr_all_filtered_latents_daybyday)$adj.r.squared
  
  MLR_weekly_COAs_adjR2 <- summary(MLR_weekly_COAs)$adj.r.squared
  
  MLR_weekly_COAs_filtered_latents_adjR2 <- summary(MLR_weekly_COAs_filtered_latents)$adj.r.squared
  #### ####
  
  
  #### CFA ####
  
  #Scale COA data and convert to ordinal data
  dat_scaled_ordinal <- dat
  dat_scaled_ordinal[,15:21] <- dat_scaled_ordinal[,15:21]/1000 #linear scale
  
  dat_scaled_ordinal[,26:37] <- lapply(dat_scaled_ordinal[,26:37],factor, order=TRUE, levels=c(0,1,2,3))
  dat_scaled_ordinal[,40:46] <- lapply(dat_scaled_ordinal[,40:46],factor, order=TRUE, levels=c(0,1,2,3,4))
  dat_scaled_ordinal[,55:61] <- lapply(dat_scaled_ordinal[,55:61],factor, order=TRUE, levels=c(0,1,2,3,4))
  
  #### Model definition ####
  
  if (condition$n_assess == 7) {
    ####Aim 1####
    
    sim.model <- '
               DM   =~ DHT_day_1 + DHT_day_2 + DHT_day_3 + DHT_day_4 + DHT_day_5 + DHT_day_6 + DHT_day_7
               wPRO =~ Item_1 + Item_2 + Item_3 + Item_4 + Item_5 + Item_6 + Item_7 + Item_8 + Item_9 + Item_10 + Item_11 + Item_12'
    
  } else if (condition$n_assess == 5) {
    
    ####Aim 1####
    
    sim.model <- '
               DM   =~ DHT_day_3 + DHT_day_4 + DHT_day_5 + DHT_day_6 + DHT_day_7
               wPRO =~ Item_1 + Item_2 + Item_3 + Item_4 + Item_5 + Item_6 + Item_7 + Item_8 + Item_9 + Item_10 + Item_11 + Item_12'
    
  } else if (condition$n_assess == 3) {
    
    ####Aim 1####
    
    sim.model <- '
               DM   =~ DHT_day_5 + DHT_day_6 + DHT_day_7
               wPRO =~ Item_1 + Item_2 + Item_3 + Item_4 + Item_5 + Item_6 + Item_7 + Item_8 + Item_9 + Item_10 + Item_11 + Item_12'
    
  } else if (condition$n_assess == 1) {
    
    ####Aim 1####
    
    sim.model <- '
               DM   =~ DHT_day_7
               wPRO =~ Item_1 + Item_2 + Item_3 + Item_4 + Item_5 + Item_6 + Item_7 + Item_8 + Item_9 + Item_10 + Item_11 + Item_12'
    
  }
  
  #### ####
  
  #### Fit the models ####
  # Aim 1 #
  
  try(
    {
      fitt <- cfa(sim.model, data = dat_scaled_ordinal, estimator = 'uls',missing='pairwise')
      
      lav_fitt<-lavInspect(fitt,"std.all")
      fitt_meas <- fitmeasures(fitt)
    }
    # ,silent=TRUE
  )
  
  #### ####
  
  ret<-c(
    
    "Avg" = Avg,
    "Avg_filtered_latents" = Avg_filtered_latents,
    
    "wPRO_regr_R2" = wPRO_regr_R2,
    "wPRO_latent_regr_R2" =  wPRO_latent_regr_R2,
    
    "all_adjR2"=all_adjR2,
    "all_filtered_latents_adjR2" = all_filtered_latents_adjR2,
    "all_daybyday_adjR2"=all_daybyday_adjR2,
    "all_filtered_latents_daybyday_adjR2" =  all_filtered_latents_daybyday_adjR2,
    "MLR_weekly_COAs_adjR2" = MLR_weekly_COAs_adjR2,
    "MLR_weekly_COAs_filtered_latents_adjR2" = MLR_weekly_COAs_filtered_latents_adjR2,
    
    "Aim_1_cor"= ifelse(exists("fitt"),lav_fitt$psi["wPRO","DM"],NA),
    "Aim_1_cfi"=ifelse(exists("fitt"),fitt_meas["cfi"],NA),
    "Aim_1_tli"=ifelse(exists("fitt"),fitt_meas["tli"],NA),
    "Aim_1_rmsea"=ifelse(exists("fitt"),fitt_meas["rmsea"],NA),
    "Aim_1_srmr"=ifelse(exists("fitt"),fitt_meas["srmr"],NA)
    
  )#,
  
  
  return(ret)
}

Summarise_with_latents_ignore_large_CFA_vals <- function(condition, results, fixed_objects = NULL) {
  
  ret<- 
    c(
      "mean_cor"=mean(results$Avg,na.rm=TRUE),
      "Emp_SE_cor"=sqrt(var(results$Avg,na.rm=TRUE)),
      "mean_filtered_latents_cor" = mean(results$Avg_filtered_latents,na.rm=TRUE),
      "EmpSE_filtered_latents_cor" = sqrt(var(results$Avg_filtered_latents,na.rm=TRUE)),
      "Mean_Emp_bias_cor" = mean(results$Avg - results$Avg_filtered_latents,na.rm=TRUE),
      
      "wPRO_regr_R2_Mean"=mean(results$wPRO_regr_R2,na.rm=TRUE),
      "wPRO_regr_R2_SE"=sqrt(var(results$wPRO_regr_R2,na.rm=TRUE)),
      "wPRO_latent_regr_R2_Mean"=mean(results$wPRO_latent_regr_R2,na.rm=TRUE),
      "wPRO_latent_regr_R2_SE"=sqrt(var(results$wPRO_latent_regr_R2,na.rm=TRUE)),
      "wPRO_regr_R2_Mean_Emp_Bias"= mean(results$wPRO_regr_R2 - results$wPRO_latent_regr_R2,na.rm=TRUE),
      
      "all_adjR2_Mean"=mean(results$all_adjR2,na.rm=TRUE),
      "all_adjR2_SE"=sqrt(var(results$all_adjR2,na.rm=TRUE)),
      "all_filtered_latents_adjR2_Mean"=mean(results$all_filtered_latents_adjR2,na.rm=TRUE),
      "all_filtered_latents_adjR2_SE"=sqrt(var(results$all_filtered_latents_adjR2,na.rm=TRUE)),
      "all_adjR2_Mean_Emp_Bias"= mean(results$all_adjR2 - results$all_filtered_latents_adjR2,na.rm=TRUE),
      "all_adjR2_Mean_Emp_Bias_targeting_daybyday"= mean(results$all_adjR2 - results$all_filtered_latents_daybyday_adjR2,na.rm=TRUE),
      
      
      "all_daybyday_adjR2_Mean"=mean(results$all_daybyday_adjR2,na.rm=TRUE),
      "all_daybyday_adjR2_SE"=sqrt(var(results$all_daybyday_adjR2,na.rm=TRUE)),
      "all_filtered_latents_daybyday_adjR2_Mean"=mean(results$all_filtered_latents_daybyday_adjR2,na.rm=TRUE),
      "all_filtered_latents_daybyday_adjR2_SE"=sqrt(var(results$all_filtered_latents_daybyday_adjR2,na.rm=TRUE)),
      "all_daybyday_adjR2_Mean_Emp_Bias" = mean(results$all_daybyday_adjR2 - results$all_filtered_latents_daybyday_adjR2,na.rm=TRUE),
      
      
      "MLR_weekly_COAs_adjR2_Mean"=mean(results$MLR_weekly_COAs_adjR2,na.rm=TRUE),
      "MLR_weekly_COAs_adjR2_SE"=sqrt(var(results$MLR_weekly_COAs_adjR2,na.rm=TRUE)),
      "MLR_weekly_COAs_filtered_latents_adjR2_Mean"=mean(results$MLR_weekly_COAs_filtered_latents_adjR2,na.rm=TRUE),
      "MLR_weekly_COAs_filtered_latents_adjR2_SE"=sqrt(var(results$MLR_weekly_COAs_filtered_latents_adjR2,na.rm=TRUE)),
      "MLR_weekly_COAs_adjR2_Mean_Emp_Bias" = mean(results$MLR_weekly_COAs_adjR2 - results$MLR_weekly_COAs_filtered_latents_adjR2,na.rm=TRUE),
      
      
      "cfa1_cor_Mean" = mean(results$Aim_1_cor,na.rm=TRUE),
      "cfa1_cor_SE" = sqrt(var(results$Aim_1_cor[abs(results$Aim_1_cor)<=2],na.rm=TRUE)),#This should be <=1 !!!
      #"cfa1_cor_SE" = sqrt(var(results$Aim_1_cor,na.rm=TRUE)),
      "cfa1_Mean_Emp_bias" = mean(results$Aim_1_cor - results$Avg_filtered_latents,na.rm=TRUE),
      
      "cfa1_cfi_Mean" = mean(results$Aim_1_cfi,na.rm=TRUE),
      "cfa1_cfi_SE" = sqrt(var(results$Aim_1_cfi,na.rm=TRUE)),
      "cfa1_cfi_acceptable_rate" = sum(results$Aim_1_cfi>=0.9, na.rm = TRUE)/sum(!is.na(results$Aim_1_cfi)),   #exclude >1 cor conditions 
      
      "cfa1_tli_Mean" = mean(results$Aim_1_tli,na.rm=TRUE),
      "cfa1_tli_SE" = sqrt(var(results$Aim_1_tli,na.rm=TRUE)),
      "cfa1_tli_acceptable_rate" = sum(results$Aim_1_tli>=0.9, na.rm = TRUE)/sum(!is.na(results$Aim_1_tli)),
      
      "cfa1_rmsea_Mean" = mean(results$Aim_1_rmsea,na.rm=TRUE),
      "cfa1_rmsea_SE" = sqrt(var(results$Aim_1_rmsea,na.rm=TRUE)),
      "cfa1_rmsea_acceptable_rate" = sum(results$Aim_1_rmsea<0.08, na.rm=TRUE)/sum(!is.na(results$Aim_1_rmsea)),
      
      "cfa1_srmr_Mean" = mean(results$Aim_1_srmr,na.rm=TRUE),
      "cfa1_srmr_SE" = sqrt(var(results$Aim_1_srmr,na.rm=TRUE)),
      "cfa1_srmr_acceptable_rate" = sum(results$Aim_1_srmr<0.08, na.rm=TRUE)/sum(!is.na(results$Aim_1_srmr)),
      
      
      "rel_precision_CFA1_vs_Pearson" = 100*((sqrt(var(results$Avg,na.rm=TRUE))/sqrt(var(results$Aim_1_cor,na.rm=TRUE)))^2-1)
      
    )
  
  
  return(ret)
  
}

Final_Results_Run_5 <- runSimulation(design=Design_5, replications=10, #500,
                                           generate=Generate_with_latent_vals, analyse=Analyse_with_latents_additional_MLR, summarise=Summarise_with_latents_ignore_large_CFA_vals,
                                           seed=c(86532:(86532+nrow(Design_5)-1)),
                                           save_results = TRUE,parallel = TRUE,packages = c('mirt','missMethods','lavaan'),beep=TRUE,
                                           control = list(print_RAM=FALSE), progress=FALSE
)

