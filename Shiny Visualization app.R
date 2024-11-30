library(shiny)
library(ggplot2)
library(dplyr)
library(SimDesign)
library(mirt)
library(lavaan)
library(tidySEM)
library(missMethods)
library(parallel)
library(rsimsum)
library(tidySEM)


#### Load the necessary functions for data generations.
#### These are identical to the contents of the Toolkit file: 
#### "Functions for Simulation Data Gen.R"


#### Generate the values of the latent trait for each individual on each day of the study ####
### n_subj: number of study participants
### fluct_sd: Daily fluctuation in an individual's latent trait. 
###   Fluctuations are Gaussian noise with this value for the SD, and mean 0.
### Returns a data frame of the latent trait values for each individal over 7 consecutive days
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


### Apply a method or perception filter to each indivdual's latent trait values.
### The filter control's the ability of the digital measure to observe the latent trait,
### or the imperfect ability of an individual to perceive their mean latent trait value.
### thetas: Each individual latent trait values over the 7 days of the study.
### filter_sd: the standard deviation of the Gaussian noise inducing the method/perception filter.
### n_subj: number of study participants
### Returns: a data frame of method/perception filtered latent traits.
Apply_filter <- function(thetas,filter_sd, n_subj) {
  #usable for DHT and RM data.
  #Usable for the theta_bars as well as the daily thetas. 
  
  return(thetas + rnorm(n_subj,0,filter_sd))
  
}


### Digital measure data generation
### thetas: Each individual latent trait values over the 7 days of the study.
### n_subj: number of study participants
### fluct_sd: Daily fluctuation in an individual's latent trait. 
### method_filter_sd: he standard deviation of the Gaussian noise inducing the digital measure method filter.
### base_rate: the hypothesized mean of the digital measure for an individual from the population with mean physical ability
### latent_effect: the proportional effect of an individual's latent physical ability on their expected digital measure count
### meas_err_sd_ratio: Magnitude of the measurement error in the digital measure (as a multiple of latent_effect)
### Returns: a data frame of the digital measure data for each individual, along with the method-filtered latent traits.
Generate_DHT_Data_Full <- function(thetas,n_subj,fluct_sd,method_filter_sd,base_rate,latent_effect,meas_err_sd_ratio) {
  
  #Apply method filter
  filtered_thetas <- Apply_filter(thetas,method_filter_sd,n_subj)
  
  #Generate Poisson means
  poisson_means <- Generate_Poisson_Means(filtered_thetas,n_subj,base_rate,latent_effect,meas_err_sd_ratio)
  
  #Generate digital measure data
  Full_data <- data.frame(filtered_thetas,Generate_DHT_Data(poisson_means))
  
  
  return(Full_data)
  
}


### Use the method-filtered latent traits to generate a Poisson mean for each participant on each study day.
### An indiviudal's digital measure data are drawn randomly from a Poisson distribution with corresponding Poisson mean generated here.
### latent_thetas:  method filtered latent traits for each individual.
### n_subj: number of study participants
### base_rate: the hypothesized mean of the digital measure for an individual from the population with mean physical ability
### latent_effect: the proportional effect of an individual's latent physical ability on their expected digital measure count
### meas_err_sd_ratio: Magnitude of the measurement error in the digital measure (as a multiple of latent_effect)
### Returns: a data frame of Poisson means for each individual on each study day.
Generate_Poisson_Means <- function(latent_thetas,n_subj,base_rate,latent_effect,meas_err_sd_ratio) {
  
  #Calculate the method filter as a multiple of latent_effect
  meas_err_sd <- meas_err_sd_ratio * latent_effect
  
  
  #Generate daily Poisson means and store in data-frame
  poisson_means <- data.frame("day_1" = base_rate + latent_effect*latent_thetas$day_1 + rnorm(n_subj,0,meas_err_sd),
                              "day_2" = base_rate + latent_effect*latent_thetas$day_2 + rnorm(n_subj,0,meas_err_sd),
                              "day_3" = base_rate + latent_effect*latent_thetas$day_3 + rnorm(n_subj,0,meas_err_sd), 
                              "day_4" = base_rate + latent_effect*latent_thetas$day_4 + rnorm(n_subj,0,meas_err_sd),
                              "day_5" = base_rate + latent_effect*latent_thetas$day_5 + rnorm(n_subj,0,meas_err_sd), 
                              "day_6" = base_rate + latent_effect*latent_thetas$day_6 + rnorm(n_subj,0,meas_err_sd),
                              "day_7" = base_rate + latent_effect*latent_thetas$day_7 + rnorm(n_subj,0,meas_err_sd))
  
  #If measurement error induces a small or negative lambda, fix it to a notional small value
  poisson_means[poisson_means<100]=100
  
  return(poisson_means)
  
}


### Generate the digital measure data for each individual on each day, by drawing a random value from a Poisson
### distribution with corresponding Poisson mean.
### poisson_means: the data frame of Poisson means for all indiviudals on all study days.
### Returns: The digital measure data.
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


### Apply MCAR data missing mechanism to the digital measure data.
### For more details on data missingness, see the Toolkit for the manuscript supplementary materials.
### DHT_raw: the digital measure data for each individual, along with the method-filtered latent traits.
### missing_rate: Proprtion of data missingness in the digital measure data
### Returns: the digital measure data with data missingness applied. 
###     (method-filtered latent traits remain in the data frame, unaffected.)
Apply_MCAR <- function(DHT_raw,missing_rate) {
  
  #Set the correct exponential distribution parameter based on missingness proprtion
  if (missing_rate==0.2) {
    
    exp_rate <- 1/20  
    
  } else if (missing_rate == 0.5) {
    
    exp_rate <- 0.18   
    
  } else if (missing_rate==0.1) {
    
    exp_rate<-1/42
    
  } else if (missing_rate==0.25) {
    
    exp_rate<-1/15
    
  } else if (missing_rate==0.4) {
    
    exp_rate<-1/8  
    
  }
  
  #Simulate study drop-out due to low adherence or device failure, using an exponential distribution.
  #The generated value from the exponential distribution for a given individual 
  #indicated after which day in the trial their data should be deleted. 
  exp_vals<-rexp(nrow(DHT_raw),rate=exp_rate)
  
  for (i in 1:nrow(DHT_raw)) {
    if (round(exp_vals[i],0)<=7 & round(exp_vals[i],0)>=1){
      DHT_raw[i,(7+c(round(exp_vals[i],0))):14]=NA
    }
  }
  
  #Simulate Day 1 logistical issues for the study, e.g. postal issues
  #Delete day 1 values completely at random using a binomial distribution.
  binom_vals<- rbinom(nrow(DHT_raw),1,missing_rate)
  
  
  for (i in 1:nrow(DHT_raw)) {
    if (binom_vals[i]==1) {
      DHT_raw[i,8] = NA
    }
  }
  
  return(DHT_raw)
  
}

#### ####

### Take the mean of the latent traits for each individual over the study period.
### thetas: Each individual latent trait values over the 7 days of the study.
### day_to_include: the number of assessment days included in the study.
### Returns: means  of the latent traits for each individual over the study period.
Generate_avg_thetas <-function(thetas, days_to_include) {
  
  return(rowMeans(thetas[,days_to_include]))
  
}


#### ####




#### Weekly PRO Data Gen ####
###Simulate the primary reference measure data - the weekly 4-response, 12-item PRO
###This function uses Item Response Theory.
### thetas_bar: Weekly mean of the latent traits for each individual
### per_filt_sd: SD of the perception filter, the imperfect ability of an individual to perceive their mean latent trait value.
### n_subj: number of study participants
### b2: Difficulty thresholds for the PRO
### reliability: Reliability of the PRO
Simulate_4_12_IRT_data_return_latents <- function(thetas_bar,per_filt_sd, n_subj,b2,reliability) {
  n_responses <- 4
  n_items <- 12
  
  #Apply perception filter to latent traits
  filtered_latents<-Apply_filter(thetas_bar,per_filt_sd,n_subj)
  
  #Generate the IRT parameters for this PRO
  cf.sim<-Generate_4_12_IRT_parameters(b2,reliability)
  
  #Generate the PRO data
  Full_COA_data<-Simulate_IRT_data(cf.sim,n_subj,n_responses,n_items,filtered_latents)
  
  #Scale the total score for each individual to a 0-100 scale
  Weekly_PRO <- rowSums(Full_COA_data)*(100/((n_responses-1)*n_items)) #(100/36)
  
  return(data.frame(filtered_latents,Full_COA_data,Weekly_PRO))
}

### Generate the IRT parameters for the primary reference measure
Generate_4_12_IRT_parameters <- function(b2,reliability) {
  
  # 4-response, 12-items
  # The b2 parameters (the ‘middle’ location parameters) consist of a series
  # of numbers between -1 and +1. 
  # The b1 parameters (the first location parameters) are based on b2 minus 1.
  # The b3 parameters (the third location parameters)  are based on b2 plus 1.
  # The a parameters are initially set at 1., then scaled to obtain a
  # dataset with the required reliability.
  
  b1 <- b2 - 1
  b3 <- b2 + 1
  a1 <- 1 #a-parameter = 1 ensures the desired reliabiilty.
  
  ## Use this if a different realibilty is desired.
  ## Adjust a-parameter to acquire the desired reliability
  #a1 <- a1 * 1.33     # Create reliability of ~0.80
  # Increasing a1 increases reliability
  
  cf.simb <- data.frame(a1,b1,b2,b3)
  
  # Transform b-parameters to d-parameters (the mirt package used in the simulation needs d-parameters)
  # difficulty (b) = easiness (d) / -a
  cf.sim <- cf.simb
  colnames(cf.sim) <- c("a1","d1","d2","d3")
  cf.sim$d1 <- -cf.simb$b1*cf.sim$a1
  cf.sim$d2 <- -cf.simb$b2*cf.sim$a1
  cf.sim$d3 <- -cf.simb$b3*cf.sim$a1
  
  return(cf.sim)
  
}


#### ####

#### General Item Response Theory (IRT) functions ####

###Generate COA data using IRT and a supplied set of parameters.
### cf.sim: difficulty thresholds/d-parameters
### n_subj: number of study participants
### n_responses: number of response option in the COA
### n_items: number of items in the COA
### latent_COA: Perception-filtered latent traits for each individual
### Returns: Each individual's responses to the COA as a data frame.
Simulate_IRT_data <- function(cf.sim,n_subj,n_responses,n_items,latent_COA) {
  
  a1 <- as.matrix(cf.sim[ , 1])
  d1 <- as.matrix(cf.sim[ , -1])
  
  #generate data using a graded response model
  dat <- simdata(a1, d1, n_subj, itemtype="graded", Theta=latent_COA)
  dat <- as.data.frame(dat)
  
  return(dat)
  
}

#### ####

#### Weekly ClinRO Data Gen ####
###Simulate the secondary reference measure data - the weekly 5-response, 7-item ClinRO
###This function uses Item Response Theory.
### thetas_bar: Weekly mean of the latent traits for each individual
### per_filt_sd: SD of the perception filter, the imperfect ability of an individual to perceive their mean latent trait value.
### n_subj: number of study participants
### b0: Difficulty thresholds for the PRO
### reliability: Reliability of the PRO
Simulate_5_7_ClinRo_IRT_Data_return_latents <- function(thetas_bar,per_filt_sd,n_subj,b0,reliability) {
  n_responses <- 5
  n_items <- 7
  
  #Apply perception filter to latent traits
  filtered_latents<-Apply_filter(thetas_bar,per_filt_sd, n_subj)
  
  #Generate the IRT parameters for this PRO
  df.sim<-Generate_5_7_ClinRO_IRT_parameters(b0,reliability)
  
  #Generate the PRO data
  Full_COA_data<-Simulate_IRT_data(df.sim,n_subj,n_responses,n_items,filtered_latents)
  
  #Scale the total score for each individual to a 0-100 scale
  Weekly_ClinRO <- rowSums(Full_COA_data)*(100/((n_responses-1)*n_items))
  
  return(data.frame(filtered_latents,Full_COA_data,Weekly_ClinRO))
}


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
  
  #b0 = c(-1.0, -0.6, -0.2, 0, 0.2, 0.6, 1.0)  #This is the vector used in the sim, repeated here for convenience
  b1 <- b0 - 1#3#1
  b2 <- b0 - 1/3#1#1/3
  b3 <- b0 + 1/3#1#1/3
  b4 <- b0 + 1#3#1
  a1 <- 1
  
  ## Adjust a-parameter to acquire the desired reliability
  a1 <- a1*.8     #Create reliability of ~0.70
  # Increasing a1 increases reliability
  # Change this if a different realibilty is desired.
  
  df.simb <- data.frame(a1,b1,b2,b3,b4)
  
  # Transform b-parameters to d-parameters (the mirt package used in the simulation needs d-parameters)
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
### Simulate the secondary reference measure data - the daily single item Patient Global Impression of Severity
### This function uses ideas in Griffiths, Pip et al. “A confirmatory factor analysis approach was found to 
### accurately estimate the reliability of transition ratings.” 
### Journal of clinical epidemiology vol. 141 (2022): 36-45. doi:10.1016/j.jclinepi.2021.08.029
### thetas: Each individual latent trait values over the 7 days of the study.
### n_subj: number of study participants
### thresholds: Thresholds at which an individual's latent trait transitions them between the response categories.
### filter_sd: SD of the perception filter for this PRO.
Sim_Daily_Single_Item_return_latents <- function(thetas, n_subj,thresholds, filter_sd) {
  
  #Apply perception filter to latent traits
  filtered_thetas<-Apply_filter(thetas,filter_sd,n_subj)
  
  #Generate the PRO data
  Daily_Single_Items<-sapply(filtered_thetas,Check_Thresholds,thrds=thresholds,N=n_subj)
  
  #Rename columns
  colnames(Daily_Single_Items) <- c('PRO_day_1','PRO_day_2','PRO_day_3','PRO_day_4','PRO_day_5','PRO_day_6','PRO_day_7')
  
  #Return PRO responses, alongside the percepion filtered latent traits for each individual.
  return(data.frame(filtered_thetas, Daily_Single_Items))
  
}

### Use the response threshold transitions to generate the PRO data for each individual.
### latent_thetas: Perception-filtered latent traits
### thrds: transition thresholds between the PRO responses
### N: number of study participants
### Returns: An individual's response to PRO, based on the where their perception-filtered latent trait fits in the threshold sequence.
Check_Thresholds <- function(latent_thetas, thrds, N) {
  
  trt <- numeric(N)
  trt[latent_thetas > thrds$thrd1] <- 1   
  trt[latent_thetas > thrds$thrd2] <- 2   
  trt[latent_thetas > thrds$thrd3] <- 3   
  trt[latent_thetas > thrds$thrd4] <- 4 
  
  return(trt)
  
}

### Generic  rounding helper function
###df: A data frame to be rounded
### digits: number of digits to round to
### Returns: a rounded data frame
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}


#### App code begins here ####


# Define UI
ui <- fluidPage(
  titlePanel("Simulating an analytical validation study"),
  sidebarLayout(
    sidebarPanel(
      numericInput("N", "Sample Size (N):", value = 100, min = 1),
      numericInput("meas_error_mag", "Magnitude of Measurement Error:", value = 0.5, min = 0, step = 0.1),
      #numericInput("n_assess", "Number of Repeated Assessments:", value = 3, min = 1),
      #selectInput("missing_method", "Data Missingness Method:", choices = c("None", "MNAR", "MCAR")),
      selectInput("missing_rate", 
                  "Data Missingness Rate:", 
                  choices = c(0, 0.1, 0.25, 0.4), multiple = TRUE),
      numericInput("replications", "Number of Replications:", value = 10, min = 1),
      actionButton("run_sim", "Run Simulation"),
      # Add selectInput for grouping
      selectInput(
        "grouping_option",
        "Group By:",
        choices = c("Missingness Mechanism", "Repeated Assessments", "Both"),
        selected = "Missingness Mechanism"
      )
    ),
    mainPanel(
      h3("Summary Output"),
      tableOutput("summary_table"),  # For `results` dataset
      # Optional: Bias plot for bias estimates
      h3("Bias Estimates Plot"),
      plotOutput("bias_plot")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Reactive expression to run the simulation based on user inputs only when button is clicked
  simulation_results <- eventReactive(input$run_sim, {
    # Display a progress bar during the simulation run
    withProgress(message = "Running simulation...", {
      
      N <- input$N
      meas_error_mag <- input$meas_error_mag
      n_assess <- input$n_assess
      missing_method <- input$missing_method
      missing_rate <- as.numeric(input$missing_rate)
      replications <- input$replications
      
      #### Design and simulation code here. This code is identical to that found in the Toolkit file:
      #### "Final Simulation"
      
      ##### Create simulation design condition matrix. Conditions are:
      #N = Sample size
      #meas_error_mag = Magnitude of the measurement error in the digital measure (as a multiple of latent_effect, defined later)
      #n_assess = number of repeated assessments being included in the analysis
      #missing_method = the digital measure data missingness method (For more details on the data missingness methods, see the manuscript supplementary materials).
      #missing_rate = Proprtion of data missingness in the digital measure data ####
      Design <- createDesign(N = N, 
                             meas_error_mag = meas_error_mag,
                             n_assess = c(1,3,5,7),
                             missing_method = c('None', 'MNAR','MCAR'), 
                             missing_rate = missing_rate,
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
        
        #Generate the digital measure data for each individual, which depend on the method filter, the ability of the digital measure to observe an individual's latent trait
        DHT_raw <- Generate_DHT_Data_Full(thetas,N,fluct_sd,meth_filt_sd,base_rate,latent_effect,meas_error_mag)
        
        #Generate data missingness in the digital measure data
        if (missing_method == "MCAR") {
          
          DHT_data <- Apply_MCAR(DHT_raw,missing_rate)
          
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
        
        # Select the method-filtered latent traits means based on which assessment days
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
            ## Mean and Standard Error of PCC (observed data and filtered latent data) ##
            "mean_cor"=mean(results$PCC,na.rm=TRUE),
            "Emp_SE_cor"=sqrt(var(results$PCC,na.rm=TRUE)),
            
            "mean_filtered_latents_cor" = mean(results$PCC_filtered_latents,na.rm=TRUE),
            "EmpSE_filtered_latents_cor" = sqrt(var(results$PCC_filtered_latents,na.rm=TRUE)),
            
            ##  ## Mean empirical bias of PCC. ##
            "Mean_Emp_bias_cor" = mean(results$PCC - results$PCC_filtered_latents,na.rm=TRUE),
            
            ## Mean and Standard Error of the SLR model (observed data and filtered latent data) ##
            "wPRO_regr_R2_Mean"=mean(results$wPRO_regr_R2,na.rm=TRUE),
            "wPRO_regr_R2_SE"=sqrt(var(results$wPRO_regr_R2,na.rm=TRUE)),
            
            "wPRO_latent_regr_R2_Mean"=mean(results$wPRO_latent_regr_R2,na.rm=TRUE),
            "wPRO_latent_regr_R2_SE"=sqrt(var(results$wPRO_latent_regr_R2,na.rm=TRUE)),
            
            ## Mean empirical bias of the SLR model. ##
            "wPRO_regr_R2_Mean_Emp_Bias"= mean(results$wPRO_regr_R2 - results$wPRO_latent_regr_R2,na.rm=TRUE),
            
            
            ## Mean and Standard Error for each MLR model (observed data and filtered latent data) ##
            
            "all_adjR2_Mean"=mean(results$all_adjR2,na.rm=TRUE),
            "all_adjR2_SE"=sqrt(var(results$all_adjR2,na.rm=TRUE)),
            
            "all_filtered_latents_adjR2_Mean"=mean(results$all_filtered_latents_adjR2,na.rm=TRUE),
            "all_filtered_latents_adjR2_SE"=sqrt(var(results$all_filtered_latents_adjR2,na.rm=TRUE)),
            
            
            "all_daybyday_adjR2_Mean"=mean(results$all_daybyday_adjR2,na.rm=TRUE),
            "all_daybyday_adjR2_SE"=sqrt(var(results$all_daybyday_adjR2,na.rm=TRUE)),
            
            "all_filtered_latents_daybyday_adjR2_Mean"=mean(results$all_filtered_latents_daybyday_adjR2,na.rm=TRUE),
            "all_filtered_latents_daybyday_adjR2_SE"=sqrt(var(results$all_filtered_latents_daybyday_adjR2,na.rm=TRUE)),
            
            
            "MLR_weekly_COAs_adjR2_Mean"=mean(results$MLR_weekly_COAs_adjR2,na.rm=TRUE),
            "MLR_weekly_COAs_adjR2_SE"=sqrt(var(results$MLR_weekly_COAs_adjR2,na.rm=TRUE)),
            
            "MLR_weekly_COAs_filtered_latents_adjR2_Mean"=mean(results$MLR_weekly_COAs_filtered_latents_adjR2,na.rm=TRUE),
            "MLR_weekly_COAs_filtered_latents_adjR2_SE"=sqrt(var(results$MLR_weekly_COAs_filtered_latents_adjR2,na.rm=TRUE)),
            
            
            ## Mean empirical bias of each MLR model. ##
            "all_adjR2_Mean_Emp_Bias"= mean(results$all_adjR2 - results$all_filtered_latents_adjR2,na.rm=TRUE),
            
            "all_adjR2_Mean_Emp_Bias_targeting_daybyday"= mean(results$all_adjR2 - results$all_filtered_latents_daybyday_adjR2,na.rm=TRUE),   
            
            "all_daybyday_adjR2_Mean_Emp_Bias" = mean(results$all_daybyday_adjR2 - results$all_filtered_latents_daybyday_adjR2,na.rm=TRUE),
            
            "MLR_weekly_COAs_adjR2_Mean_Emp_Bias" = mean(results$MLR_weekly_COAs_adjR2 - results$MLR_weekly_COAs_filtered_latents_adjR2,na.rm=TRUE),
            
            
            ## Mean, standard error, and mean empirical bias for CFA model
            "cfa1_cor_Mean" = mean(results$Aim_1_cor,na.rm=TRUE),
            "cfa1_cor_SE" = sqrt(var(results$Aim_1_cor[abs(results$Aim_1_cor)<=2],na.rm=TRUE)),
            "cfa1_Mean_Emp_bias" = mean(results$Aim_1_cor - results$PCC_filtered_latents,na.rm=TRUE),
            
            ## Mean, standard error, and rate of acceptable fit for CFI ###
            "cfa1_cfi_Mean" = mean(results$Aim_1_cfi,na.rm=TRUE),
            "cfa1_cfi_SE" = sqrt(var(results$Aim_1_cfi,na.rm=TRUE)),
            "cfa1_cfi_acceptable_rate" = sum(results$Aim_1_cfi>=0.9, na.rm = TRUE)/sum(!is.na(results$Aim_1_cfi)),   #exclude >1 cor conditions 
            
            ## Mean, standard eroor, and rate of acceptable fit for TLI ###
            "cfa1_tli_Mean" = mean(results$Aim_1_tli,na.rm=TRUE),
            "cfa1_tli_SE" = sqrt(var(results$Aim_1_tli,na.rm=TRUE)),
            "cfa1_tli_acceptable_rate" = sum(results$Aim_1_tli>=0.9, na.rm = TRUE)/sum(!is.na(results$Aim_1_tli)),
            
            ## Mean, standard eroor, and rate of acceptable fit for RMSEA ###
            "cfa1_rmsea_Mean" = mean(results$Aim_1_rmsea,na.rm=TRUE),
            "cfa1_rmsea_SE" = sqrt(var(results$Aim_1_rmsea,na.rm=TRUE)),
            "cfa1_rmsea_acceptable_rate" = sum(results$Aim_1_rmsea<0.08, na.rm=TRUE)/sum(!is.na(results$Aim_1_rmsea)),
            
            ## Mean, standard eroor, and rate of acceptable fit for SRMR ###
            "cfa1_srmr_Mean" = mean(results$Aim_1_srmr,na.rm=TRUE),
            "cfa1_srmr_SE" = sqrt(var(results$Aim_1_srmr,na.rm=TRUE)),
            "cfa1_srmr_acceptable_rate" = sum(results$Aim_1_srmr<0.08, na.rm=TRUE)/sum(!is.na(results$Aim_1_srmr)),
            
            #relative % increase in precision (CFA vs PCC)
            "rel_precision_CFA1_vs_Pearson" = 100*((sqrt(var(results$PCC,na.rm=TRUE))/sqrt(var(results$Aim_1_cor,na.rm=TRUE)))^2-1)
            
          )
        
        
        return(ret)
        
      }
      
      #####################################################################################################################
      
      #### Run the simulation. ####
      #### Refer to the simDesign package for more details on the runSimulation function.
      #### 10 replications to quickly test the principle. For full simulation, use replications = 500. 
      #### Note that 500 replications will take several hours at least, depending on the processing power of your machine.
      #### Final_Results contains the performance measures for each simulation condition from the design matrix Design
      #### along with the condition's parameters and other information.
      results <- runSimulation(design=Design, replications=10, 
                                     generate=Generate, analyse=Analyse, summarise=Summarise,
                                     seed=c(86532:(86532+nrow(Design)-1)),
                                     save_results = TRUE,parallel = TRUE,packages = c('mirt','missMethods','lavaan'),beep=TRUE,
                                     control = list(print_RAM=FALSE), progress=FALSE
      )
      
      
      ####Collate the estimates from each simmulation into one grand data frame, 
      ####and write the estimates data frame and summary results to CSV.
      
      data_set <- results
      
      #Import estimates from first condition. 
      #Estimates from each condition are saved locally based on your workspace settings.
      final_sim_res_1<- SimResults(data_set,1)
      
      #Append simulation condition parameters to estimates
      estimates_data_1<- with(final_sim_res_1,
                              data.frame(N=rep(condition$N,nrow(results)),
                                         meas_error_mag = rep(condition$meas_error_mag,nrow(results)),
                                         n_assess=rep(condition$n_assess,nrow(results)),
                                         missing_method=rep(condition$missing_method,nrow(results)),
                                         missing_rate = rep(condition$missing_rate,nrow(results)),results)
                              
      )
      
      #Add estimates and condition info from condition 1 to the grand data.frame of estimates
      collated_estimates <- estimates_data_1
      
      #loop through remaining simulation design conditions, adding estimates 
      # and condition information to the grand data.frame of estimates.                          
      for (i in 2:nrow(data_set)) {
        
        nam <- paste("final_sim_res_",i,sep="")
        
        assign(nam,SimResults(data_set,i))
        
        nam <- paste("estimates_data_",i,sep="")
        
        assign(nam,with(eval(parse(text=paste("final_sim_res_",i,sep=""))),
                        data.frame(N=rep(condition$N,nrow(results)),
                                   meas_error_mag = rep(condition$meas_error_mag,nrow(results)),
                                   n_assess=rep(condition$n_assess,nrow(results)),
                                   missing_method=rep(condition$missing_method,nrow(results)),
                                   missing_rate = rep(condition$missing_rate,nrow(results)),results))
        )
        
        collated_estimates <- rbind(collated_estimates, eval(parse(text=nam)))
        
        #Keep the workspace tidy while looping by removing unneeded variables
        rm(list = ls(pattern = paste("final_sim_res_",i,sep="")))  
        rm(list = ls(pattern = nam))
        
      }
      
    })
    # Make sure to return stuff we want to display
    return(list(results = results, collated_estimates = collated_estimates)) # return both results and collated_estimates
  })
  
  # Output collated estimates table
  # Output customized `results` table
  output$summary_table <- renderTable({
    req(simulation_results())
    simulation_results()$results %>%
      dplyr::select(
        N,
        meas_error_mag,
        n_assess,
        missing_method,
        missing_rate,
        Mean_Emp_bias_cor,
        cfa1_Mean_Emp_bias,
        all_daybyday_adjR2_Mean_Emp_Bias,
        MLR_weekly_COAs_adjR2_Mean_Emp_Bias
      ) %>%
      dplyr::rename(
        `Sample Size Condition` = N,
        `Error Condition` = meas_error_mag,
        `Repeated Assessments Condition` = n_assess,
        `Missingness Mechanism Condition` = missing_method,
        `Missingness Rate Condition` = missing_rate,
        `Pearson-based Mean Correlation Bias` = Mean_Emp_bias_cor,
        `CFA-based Mean Correlation Bias` = cfa1_Mean_Emp_bias,
        `Simple Linear Regression-based Mean Bias` = all_daybyday_adjR2_Mean_Emp_Bias,
        `Multiple Linear Regression-based Mean Bias` = MLR_weekly_COAs_adjR2_Mean_Emp_Bias
      )
  })
  
  output$bias_plot <- renderPlot({
    req(simulation_results())
    
    # Extract `collated_estimates` from `simulation_results`
    collated_data <- simulation_results()$collated_estimates
    
    # Calculate the biases
    collated_data <- collated_data %>%
      dplyr::mutate(
        `Simple Pearson-based Mean Correlation` = PCC - PCC_filtered_latents,
        `CFA-based Correlation` = Aim_1_cor - PCC_filtered_latents,
        `Simple Linear Regression` = all_daybyday_adjR2 - all_filtered_latents_daybyday_adjR2,
        `Multiple Linear Regression` = MLR_weekly_COAs_adjR2 - MLR_weekly_COAs_filtered_latents_adjR2
      )
    
    observe({
      print(head(collated_data))
    })
    
    # Gather data for box plot
    bias_data <- collated_data %>%
      tidyr::pivot_longer(
        cols = c(`Simple Pearson-based Mean Correlation`,
                 `CFA-based Correlation`,
                 `Simple Linear Regression`,
                 `Multiple Linear Regression`),
        names_to = "Bias Type",
        values_to = "Bias Estimate"
      ) %>%
      select(n_assess, missing_method, `Bias Type`, `Bias Estimate`)
    
    
    observe({
      print(head(bias_data))
    })
    
    
    # Determine grouping based on user input
    if (input$grouping_option == "Missingness Mechanism") {
      bias_data <- bias_data %>%
        mutate(Grouping = bias_data$missing_method)
    } else if (input$grouping_option == "Repeated Assessments") {
      bias_data <- bias_data %>%
        mutate(Grouping = as.factor(bias_data$n_assess))
    } else {
      bias_data <- bias_data %>%
        mutate(Grouping = interaction(bias_data$missing_method, bias_data$n_assess))
    }
    
    # Plot with dynamic grouping
    ggplot(bias_data, aes(x = `Bias Type`, y = `Bias Estimate`, fill = Grouping)) +
      geom_boxplot() +
      labs(
        x = "Bias Type",
        y = "Bias Estimate",
        title = "Bias Estimates Across Different Models"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    
  })
  
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
