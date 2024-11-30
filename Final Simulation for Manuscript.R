#### Load libraries ####
library(SimDesign)
library(mirt)
library(lavaan)
library(tidySEM)
library(missMethods)
library(parallel)
library(rsimsum)
#### ####

#### This code represents the final code employed in the accompanying manuscript. ####
#### Please see that manuscript and supplementary materials in this Toolkit for full details on the methodology employed. #### 

#This simulation is conducted using the SimDesign package: https://cran.r-project.org/web/packages/SimDesign/index.html

##### Create simulation design condition matrix. Conditions are:
#N = Sample size
#meas_error_mag = Magnitude of the measurement error in the digital measure (as a multiple of latent_effect, defined later)
#n_assess = number of repeated assessments being included in the analysis
#missing_method = the digital measure data missingness method (For more details on the data missingness methods, see the manuscript supplementary materials).
#missing_rate = Proprtion of data missingness in the digital measure data ####
Design <- createDesign(N = c(35,100,200,1000), 
                         meas_error_mag = c(0.5,1,1.5,2),
                         n_assess = c(1,3,5,7),
                         missing_method = c('None', 'MNAR','MCAR'), 
                         missing_rate = c(0,0.1,0.25,0.4),
                         subset = !((missing_method == 'None' & missing_rate != 0) | (missing_method != 'None' & missing_rate == 0)))


#####Data generation mechanism ####
Generate <- function(condition,fixed_objects=NULL) {
  

  Attach(condition)
  
  ##Generate the digital measure data ##
  
  fluct_sd = 0.3 #Daily fluctuation in an individual's latent trait. Fluctuations are Gaussian noise with this value for the SD, and mean 0.
  
  meth_filt_sd = 0.5 #Method filter - Controls the ability of the digital measure to observe the latent trait. Gaussian noise with this value dor the SD and mean zero.
  base_rate = 10000 #the hypothesized mean of the digital measure for an individual from the population with mean physical ability
  latent_effect = 1250 #the proportional effect of an individual's latent physical ability on their expected digital measure count

  #Simulate the values of the latent trait for each individual on each day of the study
  thetas<- Generate_Daily_Thetas(N,fluct_sd) 

  #Generate the digital measure data for each individual, which depend on the ability of the digital measure to observe an individual's latent trait
  DHT_raw <- Generate_DHT_Data_Full(thetas,N,fluct_sd,meth_filt_sd,base_rate,latent_effect,meas_error_mag)
  
  #Generate data missingness in the digital measure data
  if (missing_method == "MCAR") {
       
    DHT_data <- Apply_MCAR_Missing_Data_with_latents(DHT_raw,missing_rate)
     
  } else if (missing_method == "MNAR") {
        
    DHT_data <- data.frame(DHT_raw[1:7],
                           delete_MNAR_1_to_x(DHT_raw[8:14],
                                              missing_rate,
                                              c("DHT_day_1", "DHT_day_2", "DHT_day_3",
                                                "DHT_day_4", "DHT_day_5", "DHT_day_6", "DHT_day_7"),
                                              x=4)
    )
    
    
  } else {
    
    DHT_data <- DHT_raw
    
  }
  
  #Vector representing the specific days of digital measure data to include in assessment (DHT-side only; COAs always use 7 days)
  days_to_include <- c((8-n_assess):7)
    
  # Calculate each individual's mean digital measure response, based on which assessment days are being included in the assessment
  if (n_assess != 1) {
    
    DHT_Means <- rowMeans(DHT_data[,7+days_to_include],na.rm = TRUE)
  } else {
    DHT_Means <- DHT_data[,14]
  }
  
  # Calculate the method-filtered latent traits means based on which assessment days
  # are being included in the assessment
  if (n_assess != 1) {
    
    filtered_DHT_latent_Means <- rowMeans(DHT_data[,days_to_include],na.rm = TRUE)
  } else {
    filtered_DHT_latent_Means <- DHT_data[,7]
  }
  
  #### ####
  
  #Generate the mean latent trait over the assessment period.
  #This value is used in the Weekly COA data generaion. We're always using 7 days for the mean, 
  #no matter the number of digital measure assessment days included, as this is the typical recall period for COAs.
  days_to_include <- c(1:7)
  thetas_bar <- Generate_avg_thetas(thetas,days_to_include)

  
  ####Generate the COA data ####
  #Weekly COA data is generated using item response theory.

  
  #### Generate the primary reference measure data - 
  #### a 12-item, 4-response PRO with a recall period of one week. ####
  
  #Set perception filter - the imperfect ability of an individual to perceive their mean latent trait value
  #when recalling their activities during the previous seven days.
  per_filt_sd = 1   
    
  #Set the difficulty threshold parameters
  b2  <- c(-1.0, -0.8, -0.8, -0.4, -0.4, 0, 0, 0.4, 0.4, 0.8, 0.8, 1.0)
  
  #Set reliability
  reliability<-0.8

  #Generate the PRO data. This is a total score for each individual, scaled to a 0-100 scale
  Weekly_PRO <-Simulate_4_12_IRT_data_return_latents(thetas_bar,per_filt_sd,N,b2,reliability)
  
  ####    ####
  
  ####Generate the secondary reference measure data - 
  #### a 7-item, 5-response ClinRO with a recall period of one week. ####

  #Set perception filter - the imperfect ability of an individual to perceive their mean latent trait value
  #when recalling their activities during the previous seven days.
  per_Filt_sd = 1 

  #Set the difficulty threshold parameters
  b0 = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75) 

  #Set reliabiity
  reliability = 0.7    

  #Generate the ClinRO data
  ClinRO_Data <- Simulate_5_7_ClinRo_IRT_Data_return_latents(thetas_bar,per_Filt_sd,N,b0,reliability)
    
  ####    ####
  
  #### Generate secondary reference measure data - a Daily Single-item five-point PRO 
  #### representing a patient's global impression of severity. ####
  
  # This data is generated using concepts in Griffiths, Pip et al. 
  #“A confirmatory factor analysis approach was found to accurately estimate the reliability of transition ratings.” 
  #Journal of clinical epidemiology vol. 141 (2022): 36-45. doi:10.1016/j.jclinepi.2021.08.029

  #Set perception filter - the imperfect ability of an individual to perceive their mean latent trait value
  #when recalling their activities during the previous seven days.
  per_filt_sd = 0.5 
  
  # Create iMIC distribution and other thresholds
  imic <- rnorm(N, 0.5, 0.075) 
  thrd1 <- rnorm(N, -1.5, 0.075)
  thrd2 <- rnorm(N, -0.5, 0.075)
  thrd3 <- imic
  thrd4 <- rnorm(N, 1.5, 0.075)
  
  thrds<- data.frame(thrd1,thrd2,thrd3,thrd4)
   
  #Generate the daily PRO data
  Daily_PRO <- Sim_Daily_Single_Item_return_latents(thetas,N,thrds,per_filt_sd)

  
  
  #Vector representing the specific days of daily PRO measure data to include in assessment (DHT-side only; weekly COAs always use 7 days)
  days_to_include <- c(1:7)
   
  # Calculate each individual's mean daily PRO reference measure response, based on which assessment days are being included in the assessment.
  if (n_assess != 1) {
    Daily_PRO_means <- rowMeans(Daily_PRO[,7+days_to_include])
  } else {
    Daily_PRO_means <- Daily_PRO[,14]
  }

  ####    ####

  # Calculate the perception-filtered latent traits means based on which assessment days
  #are being included in the assessment
  if (n_assess != 1) {
    filtered_dPRO_latent_Means <- rowMeans(Daily_PRO[,days_to_include],na.rm = TRUE)
  } else {
    filtered_dPRO_latent_Means <- Daily_PRO[,7]
  }
  
  
  #Combine all generated data into a data frame, rename columns and return it.
  dat<-data.frame(thetas, DHT_data,filtered_DHT_latent_Means,DHT_Means,
                  thetas_bar,Weekly_PRO,ClinRO_Data,Daily_PRO,filtered_dPRO_latent_Means,
                  Daily_PRO_means,subject=c(1:N))

  
  names(dat)[1:7]<-c("theta_day_1","theta_day_2","theta_day_3","theta_day_4",
                     "theta_day_5","theta_day_6","theta_day_7")
  names(dat)[8:14]<-c("filtered_theta_day_1","filtered_theta_day_2","filtered_theta_day_3",
                      "filtered_theta_day_4", "filtered_theta_day_5","filtered_theta_day_6",
                      "filtered_theta_day_7")
  
  return(dat)
}

#############################################################################################

#### Data analysis function ####
Analyse <- function(condition, dat, fixed_objects = NULL) {
  
  #Calculate Pearson correlation coefficient between the digital measure and primary reference measure data
  PCC <- cor(dat$DHT_Means,dat$Weekly_PRO,use="na.or.complete")
  
  #Calculate Pearson correlation coefficient for filtered latent data for the digital measure and primary reference measure
  PCC_filtered_latents <- with(dat,cor(filtered_DHT_latent_Means,filtered_latents,use="na.or.complete"))
  
  #### Simple Linear Regression Models ####
  
  #Mean digital measure data as outcome, Weekly PRO scaled total score as the predictor.
  wPRO_regr <- lm(DHT_Means~Weekly_PRO,data=dat)
  wPRO_regr_R2 <- summary(wPRO_regr)$r.squared  #Calculate r-squared statistic

  #Equivalent model for the method/perception filtered latent variables
  wPRO_latent_regr<- lm(filtered_DHT_latent_Means~filtered_latents,data=dat)
  wPRO_latent_regr_R2<- summary(wPRO_latent_regr)$r.squared  #Calculate r-squared statistic

 
  #### Multiple linear regression models ####
  #### all models use the mean digital measure data as the outcome. ####

  #Include all reference measures as predictors. Weekly COAs as scaled total scores,
  #and weekly mean of the daily reference measures.
  regr_all <- lm(DHT_Means~Weekly_PRO+Weekly_ClinRO+Daily_PRO_means,data=dat)
  
  #Equivalent model for the method/perception filtered latent variables
  regr_all_filtered_latents <- lm(filtered_DHT_latent_Means~filtered_latents+filtered_latents.1
                                  +filtered_dPRO_latent_Means,data=dat)


  
  #As above, but each individual day of daily PRO data is included as a separate predictor.
  regr_all_daybyday <- lm(DHT_Means~Weekly_PRO+Weekly_ClinRO
                          +PRO_day_1+PRO_day_2+PRO_day_3+PRO_day_4+PRO_day_5+PRO_day_6+PRO_day_7,data=dat)
  
  #Equivalent model for the method/perception filtered latent variables
  regr_all_filtered_latents_daybyday <- lm(filtered_DHT_latent_Means~filtered_latents+filtered_latents.1
                                           +day_1.2+day_2.2+day_3.2+day_4.2+day_5.2+day_6.2+day_7.2,data=dat)



  
  #Include just the weekly COAs as predictors, in the same way as above, and exclude the daily PRO.
  MLR_weekly_COAs <- lm(DHT_Means~Weekly_PRO+Weekly_ClinRO,data=dat)
  
 #Equivalent model for the method/perception filtered latent variables
  MLR_weekly_COAs_filtered_latents <- lm(filtered_DHT_latent_Means~filtered_latents+filtered_latents.1,data=dat)
  
  
  
  ####Calculate adjusted R-squareds for the MLR models####
  ## Daily PRO as a single variable ##
  
  #Measured values
  all_adjR2 <- summary(regr_all)$adj.r.squared
   #filtered latents
  all_filtered_latents_adjR2 <- summary(regr_all_filtered_latents)$adj.r.squared

  
  ## Each day of daily PRO as a separate variable ##
  all_daybyday_adjR2 <- summary(regr_all_daybyday)$adj.r.squared
  #filtered latents
  all_filtered_latents_daybyday_adjR2 <- summary(regr_all_filtered_latents_daybyday)$adj.r.squared

  #Just the weekly COAs, no daily PRO
  MLR_weekly_COAs_adjR2 <- summary(MLR_weekly_COAs)$adj.r.squared
  #filtered latents
  MLR_weekly_COAs_filtered_latents_adjR2 <- summary(MLR_weekly_COAs_filtered_latents)$adj.r.squared
  
  #### ####
  
  
  #### CFA models ####
  
  #Scale COA data to match the weekly PRO scale, and convert to ordinal data.
  dat_scaled_ordinal <- dat
  dat_scaled_ordinal[,15:21] <- dat_scaled_ordinal[,15:21]/1000 #linear scaling
  
  dat_scaled_ordinal[,26:37] <- lapply(dat_scaled_ordinal[,26:37],factor, order=TRUE, levels=c(0,1,2,3))
  dat_scaled_ordinal[,40:46] <- lapply(dat_scaled_ordinal[,40:46],factor, order=TRUE, levels=c(0,1,2,3,4))
  dat_scaled_ordinal[,55:61] <- lapply(dat_scaled_ordinal[,55:61],factor, order=TRUE, levels=c(0,1,2,3,4))

  
  #### Model definitions ####
  #### Two-factor mdoel with correlated factors
  #### Digital measure data as one factor, primary reference measure (i.e., the weekly PRO) data as the other factor. ####
  #### A different model is required based on the number of repeated assessments included from the digital measure data. ####
  
  if (condition$n_assess == 7) {
       
    sim.model <- '
               DM   =~ DHT_day_1 + DHT_day_2 + DHT_day_3 + DHT_day_4 + DHT_day_5 + DHT_day_6 + DHT_day_7
               wPRO =~ Item_1 + Item_2 + Item_3 + Item_4 + Item_5 + Item_6 + Item_7 + Item_8 + Item_9 + Item_10 + Item_11 + Item_12'
    
  } else if (condition$n_assess == 5) {
   
    sim.model <- '
               DM   =~ DHT_day_3 + DHT_day_4 + DHT_day_5 + DHT_day_6 + DHT_day_7
               wPRO =~ Item_1 + Item_2 + Item_3 + Item_4 + Item_5 + Item_6 + Item_7 + Item_8 + Item_9 + Item_10 + Item_11 + Item_12'
    
  } else if (condition$n_assess == 3) {
     
    sim.model <- '
               DM   =~ DHT_day_5 + DHT_day_6 + DHT_day_7
               wPRO =~ Item_1 + Item_2 + Item_3 + Item_4 + Item_5 + Item_6 + Item_7 + Item_8 + Item_9 + Item_10 + Item_11 + Item_12'
    
  } else if (condition$n_assess == 1) {
    
    
    sim.model <- '
               DM   =~ DHT_day_7
               wPRO =~ Item_1 + Item_2 + Item_3 + Item_4 + Item_5 + Item_6 + Item_7 + Item_8 + Item_9 + Item_10 + Item_11 + Item_12'
    
  }
  
  #### ####

  
  #### Fit the models. ####
  
  try(
    {
      fitt <- cfa(sim.model, data = dat_scaled_ordinal, estimator = 'uls',missing='pairwise')
      
      lav_fitt<-lavInspect(fitt,"std.all")
      fitt_meas <- fitmeasures(fitt)
    }
   
  )
  
  #### ####


  #### Collate the agreement statistics from each method, and return. #### 
  ret<-c(

    #PCCs
    "PCC" = PCC,
    "PCC_filtered_latents" = PCC_filtered_latents,

    #SLR R-squareds
    "wPRO_regr_R2" = wPRO_regr_R2,
    "wPRO_latent_regr_R2" =  wPRO_latent_regr_R2,

    #MLR adjusted R-squareds
    "all_adjR2"=all_adjR2,
    "all_filtered_latents_adjR2" = all_filtered_latents_adjR2,
    "all_daybyday_adjR2"=all_daybyday_adjR2,
    "all_filtered_latents_daybyday_adjR2" =  all_filtered_latents_daybyday_adjR2,
    "MLR_weekly_COAs_adjR2" = MLR_weekly_COAs_adjR2,
    "MLR_weekly_COAs_filtered_latents_adjR2" = MLR_weekly_COAs_filtered_latents_adjR2,

    #CFA factor correlation
    "Aim_1_cor"= ifelse(exists("fitt"),lav_fitt$psi["wPRO","DM"],NA),

    #CFA model fit statistics
    "Aim_1_cfi"=ifelse(exists("fitt"),fitt_meas["cfi"],NA),  #Comparative Fit Index
    "Aim_1_tli"=ifelse(exists("fitt"),fitt_meas["tli"],NA),  #Tucker-Lewis Index
    "Aim_1_rmsea"=ifelse(exists("fitt"),fitt_meas["rmsea"],NA),  #Root Mean Square Error of Approximation
    "Aim_1_srmr"=ifelse(exists("fitt"),fitt_meas["srmr"],NA)  #Standardized Root Mean Square Residual
    
  )
  
  return(ret)
}

#######################################################################################################

#### Calculate the simulation performance measures. #####
Summarise <- function(condition, results, fixed_objects = NULL) {
  
  ret<- 
    c(
      "mean_cor"=mean(results$PCC,na.rm=TRUE),
      "Emp_SE_cor"=sqrt(var(results$PCC,na.rm=TRUE)),
      "mean_filtered_latents_cor" = mean(results$PCC_filtered_latents,na.rm=TRUE),
      "EmpSE_filtered_latents_cor" = sqrt(var(results$PCC_filtered_latents,na.rm=TRUE)),
      "Mean_Emp_bias_cor" = mean(results$PCC - results$PCC_filtered_latents,na.rm=TRUE),
      
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
      "cfa1_cor_SE" = sqrt(var(results$Aim_1_cor[abs(results$Aim_1_cor)<=2],na.rm=TRUE)),
      "cfa1_Mean_Emp_bias" = mean(results$Aim_1_cor - results$PCC_filtered_latents,na.rm=TRUE),
      
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

      #relative % increase in precision (CFA vs PCC
      "rel_precision_CFA1_vs_Pearson" = 100*((sqrt(var(results$PCC,na.rm=TRUE))/sqrt(var(results$Aim_1_cor,na.rm=TRUE)))^2-1)
      
    )
  
  
  return(ret)
  
}

#### Run the simulation. ####
#### Refer to the simDesign package for more details on the runSimulation function.
#### 10 replications to quickly test the principle. For full simulation, use replications = 500. 
#### Note that 500 replications will take several hours at least, depending on the processing power of your machine.
Final_Results <- runSimulation(design=Design, replications=10, 
                                           generate=Generate, analyse=Analyse, summarise=Summarise,
                                           seed=c(86532:(86532+nrow(Design)-1)),
                                           save_results = TRUE,parallel = TRUE,packages = c('mirt','missMethods','lavaan'),beep=TRUE,
                                           control = list(print_RAM=FALSE), progress=FALSE
)

