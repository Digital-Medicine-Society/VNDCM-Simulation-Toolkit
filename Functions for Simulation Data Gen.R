
#### Generate latent trait values ####
#Yes#
Generate_Daily_Thetas <- function(n_subj, fluct_sd) {
  
  #Create starting state of latent trait
  day_1 <- rnorm(n_subj)
  day_1 <- scale(day_1)
  
  #Generate latent trait on subsequent days
  day_2 <- day_1 + rnorm(n_subj,0,fluct_sd)
  day_2 <- scale(day_2)
  
  day_3 <- day_2 + rnorm(n_subj,0,fluct_sd)
  day_3 <- scale(day_3)
  
  day_4 <- day_3 + rnorm(n_subj,0,fluct_sd)
  day_4 <- scale(day_4)
  
  day_5 <- day_4 + rnorm(n_subj,0,fluct_sd)
  day_5 <- scale(day_5)
  
  day_6 <- day_5 + rnorm(n_subj,0,fluct_sd)
  day_6 <- scale(day_6)
  
  day_7 <- day_6 + rnorm(n_subj,0,fluct_sd)
  day_7 <- scale(day_7)
  
  return(data.frame(day_1,day_2,day_3,day_4,day_5,day_6,day_7))
  
}
#### ####

#Yes#
#### Multi-purpose ####
Apply_filter <- function(thetas,filter_sd, n_subj) {
  #usable for DHT and RM data.
  #Usable for the theta_bars as well as the daily thetas. 

  #Set method/perception filter
  # THE METHOD/PERCEPTION FILTER IS NOT ASSUMED TO RANDOMLY CHANGE FROM DAY TO DAY.
  # THIS IS A 'FIXED' FILTER.
  # THEREFORE, THE SAME FILTER IS TO BE APPLIED EVERY DAY.
  #filter <- rnorm(n_subj,0,filter_sd)   
  
   #Generate filtered latent trait
  #filtered_thetas <- thetas + filter
  
  return(thetas + rnorm(n_subj,0,filter_sd))
  #return(filtered_thetas)
  
}


#### DHT Data Gen ####
#Yes#
Generate_DHT_Data_Full <- function(thetas,n_subj,fluct_sd,method_filter_sd,base_rate,latent_effect,meas_err_sd_ratio) {
  
  filtered_thetas <- Apply_filter(thetas,method_filter_sd,n_subj)
  
  poisson_means <- Generate_Poisson_Means(filtered_thetas,n_subj,base_rate,latent_effect,meas_err_sd_ratio)
  
  Full_data <- data.frame(filtered_thetas,Generate_DHT_Data(poisson_means))
  
  
  return(Full_data)
  
}

#Yes#
Generate_Poisson_Means <- function(latent_thetas,n_subj,base_rate,latent_effect,meas_err_sd_ratio) {
  
  meas_err_sd <- meas_err_sd_ratio * latent_effect
  

  #Generate daily Poisson means and store in data-frame
  poisson_means <- data.frame("day_1" = base_rate + latent_effect*latent_thetas$day_1 + rnorm(n_subj,0,meas_err_sd),
                              "day_2" = base_rate + latent_effect*latent_thetas$day_2 + rnorm(n_subj,0,meas_err_sd),
                              "day_3" = base_rate + latent_effect*latent_thetas$day_3 + rnorm(n_subj,0,meas_err_sd), 
                              "day_4" = base_rate + latent_effect*latent_thetas$day_4 + rnorm(n_subj,0,meas_err_sd),
                              "day_5" = base_rate + latent_effect*latent_thetas$day_5 + rnorm(n_subj,0,meas_err_sd), 
                              "day_6" = base_rate + latent_effect*latent_thetas$day_6 + rnorm(n_subj,0,meas_err_sd),
                              "day_7" = base_rate + latent_effect*latent_thetas$day_7 + rnorm(n_subj,0,meas_err_sd))
  
  #If measurement error induces a small or negative lambda, fix it at a notional small value
  poisson_means[poisson_means<100]=100
  
  return(poisson_means)
  
}


#Yes#
Generate_DHT_Data <- function(poisson_means) {

  #Create empty vectors
  pdht_day_1 = integer(length(poisson_means$day_1))
  pdht_day_2 = integer(length(poisson_means$day_2))
  pdht_day_3 = integer(length(poisson_means$day_3))
  pdht_day_4 = integer(length(poisson_means$day_4))
  pdht_day_5 = integer(length(poisson_means$day_5))
  pdht_day_6 = integer(length(poisson_means$day_6))
  pdht_day_7 = integer(length(poisson_means$day_7))
  
  #Generate simulated data using the Poisson means for each day
  for (i in 1:nrow(poisson_means)){
    pdht_day_1[i] = rpois(1,poisson_means$day_1[i])
    pdht_day_2[i] = rpois(1,poisson_means$day_2[i])
    pdht_day_3[i] = rpois(1,poisson_means$day_3[i])
    pdht_day_4[i] = rpois(1,poisson_means$day_4[i])
    pdht_day_5[i] = rpois(1,poisson_means$day_5[i])
    pdht_day_6[i] = rpois(1,poisson_means$day_6[i])
    pdht_day_7[i] = rpois(1,poisson_means$day_7[i])
  }
  
  #Copy over to dataframe
  pdht_data = data.frame("DHT_day_1" = pdht_day_1, "DHT_day_2" = pdht_day_2,
                         "DHT_day_3" = pdht_day_3, "DHT_day_4" = pdht_day_4,
                         "DHT_day_5" = pdht_day_5, "DHT_day_6" = pdht_day_6,
                         "DHT_day_7" = pdht_day_7)
  

  return(pdht_data)
  
}


#Yes#
Apply_MCAR_Missing_Data_with_latents <- function(DHT_raw,missing_rate) {
  
  if (missing_rate==0.2) {
    
    exp_rate <- 1/20     #1/19 #0.2 Missing data
    
  } else if (missing_rate == 0.5) {
    
    exp_rate <- 0.18    #.215 #0.5 missing data
    
  } else if (missing_rate==0.1) {
    
    exp_rate<-1/42
    
  } else if (missing_rate==0.25) {
    
    exp_rate<-1/15
    
  } else if (missing_rate==0.4) {
    
    exp_rate<-1/8  
    
  }
  
  #Simulate drop out from low adherence or device failure
  
  exp_vals<-rexp(nrow(DHT_raw),rate=exp_rate)
  
  for (i in 1:nrow(DHT_raw)) {
    if (round(exp_vals[i],0)<=7 & round(exp_vals[i],0)>=1){
      DHT_raw[i,(7+c(round(exp_vals[i],0))):14]=NA
    }
  }
  
  #Simulate Day 1 logistic issues
  binom_vals<- rbinom(nrow(DHT_raw),1,missing_rate)#0.16)  #0.2 missing rate
  
  
  for (i in 1:nrow(DHT_raw)) {
    if (binom_vals[i]==1) {
      DHT_raw[i,8] = NA
    }
  }
  
  # print(sum(is.na(DHT_raw))/(nrow(DHT_raw)*7))
  # for (i in 1:7) {
  # print(sum(is.na(DHT_raw)[,i])/nrow(DHT_raw))
  # }
  
  return(DHT_raw)
  
}

#### ####


Generate_avg_thetas <-function(thetas, days_to_include) {
  
  return(rowMeans(thetas[,days_to_include]))
  
}


round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}
#### ####

#### General Item Response Theory (IRT) functions ####
#Yes#
Simulate_IRT_data <- function(cf.sim,n_subj,n_responses,n_items,latent_COA) {
  
  a1 <- as.matrix(cf.sim[ , 1])
  d1 <- as.matrix(cf.sim[ , -1])
  
  dat <- simdata(a1, d1, n_subj, itemtype="graded", Theta=latent_COA)
  dat <- as.data.frame(dat)
  
  return(dat)
 
}

#### ####


#### Weekly PRO Data Gen ####

#Yes#
Simulate_4_12_IRT_data_return_latents <- function(thetas_bar,per_filt_sd, n_subj,b2,reliability) {
  n_responses <- 4
  n_items <- 12
  
  filtered_latents<-Apply_filter(thetas_bar,per_filt_sd,n_subj)
  
  cf.sim<-Generate_4_12_IRT_parameters(b2,reliability)
  
  Full_COA_data<-Simulate_IRT_data(cf.sim,n_subj,n_responses,n_items,filtered_latents)
  
  Weekly_PRO <- rowSums(Full_COA_data)*(100/((n_responses-1)*n_items)) #(100/36)
  
  return(data.frame(filtered_latents,Full_COA_data,Weekly_PRO))
}

#Yes#
Generate_4_12_IRT_parameters <- function(b2,reliability) {
  
  b1 <- b2 - 1
  b3 <- b2 + 1
  a1 <- 1
  
  ## Adjust a-parameter to acquire the desired reliability
  #a1 <- a1 * 1.33     # Create reliability of ~0.80
  # Increasing a1 increases reliability
  
  cf.simb <- data.frame(a1,b1,b2,b3)
  #print(cf.simb)
  
  # Transform b-parameters to d-parameters (mirt works with d-parameters)
  # difficulty (b) = easiness (d) / -a
  cf.sim <- cf.simb
  colnames(cf.sim) <- c("a1","d1","d2","d3")
  cf.sim$d1 <- -cf.simb$b1*cf.sim$a1
  cf.sim$d2 <- -cf.simb$b2*cf.sim$a1
  cf.sim$d3 <- -cf.simb$b3*cf.sim$a1
  
  return(cf.sim)
  
}


#### ####

#### Weekly ClinRO Data Gen ####
#Yes#
Simulate_5_7_ClinRo_IRT_Data_return_latents <- function(thetas_bar,per_filt_sd,n_subj,b0,reliability) {
  n_responses <- 5
  n_items <- 7
  
  filtered_latents<-Apply_filter(thetas_bar,per_filt_sd, n_subj)
  
  df.sim<-Generate_5_7_ClinRO_IRT_parameters(b0,reliability)
  
  Full_COA_data<-Simulate_IRT_data(df.sim,n_subj,n_responses,n_items,filtered_latents)
  
  Weekly_ClinRO <- rowSums(Full_COA_data)*(100/((n_responses-1)*n_items))
  
  return(data.frame(filtered_latents,Full_COA_data,Weekly_ClinRO))
}

#Yes#
Generate_5_7_ClinRO_IRT_parameters <- function(b0, reliability) {
  
  # 5-response, 7-items
  # The b2 parameters (the ‘middle’ location parameters) consist of a series
  # of numbers between -1 and +1. 
  # The b1 parameters (the first location parameters) are based on b2 minus 1.
  # The b3 parameters (the third location parameters)  are based on b2 plus 1.
  # Input vec b0 is used only to calc the actual b params used in the IRT model
  # The a parameters are initially set at 1.
  # Secondarily, the a parameters are multiplied by a factor to obtain a
  # dataset with the required reliability.

  #b0 = c(-1.0, -0.6, -0.2, 0, 0.2, 0.6, 1.0)  #This is the vector used in the sim, pasted here as a comment for convenience
  b1 <- b0 - 1#3#1
  b2 <- b0 - 1/3#1#1/3
  b3 <- b0 + 1/3#1#1/3
  b4 <- b0 + 1#3#1
  a1 <- 1
  
  ## Adjust a-parameter to acquire the desired reliability
  a1 <- a1*.8 #2     #Create reliability of ~0.70
  # Increasing a1 increases reliability
  
  df.simb <- data.frame(a1,b1,b2,b3,b4)
  #print(df.simb)
  
  # Transform b-parameters to d-parameters (mirt works with d-parameters)
  # difficulty (b) = easiness (d) / -a
  df.sim <- df.simb
  colnames(df.sim) <- c("a1","d1","d2","d3","d4")
  df.sim$d1 <- -df.simb$b1*df.sim$a1
  df.sim$d2 <- -df.simb$b2*df.sim$a1
  df.sim$d3 <- -df.simb$b3*df.sim$a1
  df.sim$d4 <- -df.simb$b4*df.sim$a1
  
  return(df.sim)
  
}

#### ####


#### Daily PRO Data Gen ####

#Yes#
Sim_Daily_Single_Item_return_latents <- function(thetas, n_subj,thresholds, filter_sd) {
  
  filtered_thetas<-Apply_filter(thetas,filter_sd,n_subj)
  
  Daily_Single_Items<-sapply(filtered_thetas,Check_Thresholds,thrds=thresholds,N=n_subj)
  
  colnames(Daily_Single_Items) <- c('PRO_day_1','PRO_day_2','PRO_day_3','PRO_day_4','PRO_day_5','PRO_day_6','PRO_day_7')
  
  return(data.frame(filtered_thetas, Daily_Single_Items))
  
}

#Yes#
Check_Thresholds <- function(latent_thetas, thrds, N) {
  
  trt <- numeric(N)
  trt[latent_thetas > thrds$thrd1] <- 1   
  trt[latent_thetas > thrds$thrd2] <- 2   
  trt[latent_thetas > thrds$thrd3] <- 3   
  trt[latent_thetas > thrds$thrd4] <- 4 
  
  return(trt)
  
}

