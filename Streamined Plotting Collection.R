library(bigrquery)
library(ggplot2)
library(ggpattern)
library(dplyr)
library(data.table)
library(SimDesign)

options(repr.plot.width=12, repr.plot.height=8)

#Function to Resummarize the results, in case they are not stored in your local environment
Resumm_exclude_CFA_cor_out_of_range <- function(condition, results, fixed_objects = NULL) {
  
  results_cfa <- results %>%
    subset(abs(results$Aim_1_cor)<=1)
  
  ret<-
    c(
      "mean_cor"=mean(results$Avg,na.rm=TRUE),
      "Emp_SE_cor"=sqrt(var(results$Avg,na.rm=TRUE)),
      "mean_filtered_latents_cor" = mean(results$Avg_filtered_latents,na.rm=TRUE),
      "EmpSE_filtered_latents_cor" = sqrt(var(results$Avg_filtered_latents,na.rm=TRUE)),
      
      "wPRO_regr_R2_Mean"=mean(results$wPRO_regr_R2,na.rm=TRUE),
      "wPRO_regr_R2_SE"=sqrt(var(results$wPRO_regr_R2,na.rm=TRUE)),
      
      "all_adjR2_Mean"=mean(results$all_adjR2,na.rm=TRUE),
      "all_adjR2_SE"=sqrt(var(results$all_adjR2,na.rm=TRUE)),
      
      "all_daybyday_adjR2_Mean"=mean(results$all_daybyday_adjR2,na.rm=TRUE),
      "all_daybyday_adjR2_SE"=sqrt(var(results$all_daybyday_adjR2,na.rm=TRUE)),
      
      "MLR_weekly_COAs_adjR2_Mean"=mean(results$MLR_weekly_COAs_adjR2,na.rm=TRUE),
      "MLR_weekly_COAs_adjR2_SE"=sqrt(var(results$MLR_weekly_COAs_adjR2,na.rm=TRUE)),
      
      "cfa1_cor_Mean" = mean(results_cfa$Aim_1_cor,na.rm=TRUE),
      "cfa1_cor_SE" = sqrt(var(results_cfa$Aim_1_cor,na.rm=TRUE)),
      #"cfa1_cor_SE" = sqrt(var(results$Aim_1_cor[abs(results$Aim_1_cor)<=1],na.rm=TRUE)),
      
      "cfa1_cfi_Mean" = mean(results_cfa$Aim_1_cfi,na.rm=TRUE),
      "cfa1_cfi_SE" = sqrt(var(results_cfa$Aim_1_cfi,na.rm=TRUE)),
      "cfa1_cfi_acceptable_rate" = sum(results_cfa$Aim_1_cfi>=0.9, na.rm = TRUE)/sum(!is.na(results_cfa$Aim_1_cfi)),   #exclude >1 cor conditions
      
      "cfa1_tli_Mean" = mean(results_cfa$Aim_1_tli,na.rm=TRUE),
      "cfa1_tli_SE" = sqrt(var(results_cfa$Aim_1_tli,na.rm=TRUE)),
      "cfa1_tli_acceptable_rate" = sum(results_cfa$Aim_1_tli>=0.9, na.rm = TRUE)/sum(!is.na(results_cfa$Aim_1_tli)),
      
      "cfa1_rmsea_Mean" = mean(results_cfa$Aim_1_rmsea,na.rm=TRUE),
      "cfa1_rmsea_SE" = sqrt(var(results_cfa$Aim_1_rmsea,na.rm=TRUE)),
      "cfa1_rmsea_acceptable_rate" = sum(results_cfa$Aim_1_rmsea<0.08, na.rm=TRUE)/sum(!is.na(results_cfa$Aim_1_rmsea)),
      
      "cfa1_srmr_Mean" = mean(results_cfa$Aim_1_srmr,na.rm=TRUE),
      "cfa1_srmr_SE" = sqrt(var(results_cfa$Aim_1_srmr,na.rm=TRUE)),
      "cfa1_srmr_acceptable_rate" = sum(results_cfa$Aim_1_srmr<0.08, na.rm=TRUE)/sum(!is.na(results_cfa$Aim_1_srmr)),
      
      "cfa1_cor_OOB" = nrow(results) - nrow(results_cfa)
      
    )
  
  
  return(ret)
  
}


## Estimates ##

#### Collect the estimates ####
#This will prompt you to choose a file. Choose the file "Collated Estimates.csv"
#saved as the output from "Final Simulation for Manuscript.R"
est <- read.csv(file.choose())[2:21] 

#Rename columns
names(est)[1]<-"N"
names(est)[2]<-"MEM"
names(est)[3]<-'RA'

#Treat variables as factors
est$N <- as.factor(est$N)
est$MEM <- as.factor(est$MEM)
est$RA <- as.factor(est$RA)
est$missing_rate <- as.factor(est$missing_rate)

#remove conditions where CFA correlation is outside [-1,1]
est[abs(est$Aim_1_cor)>1,]$Aim_1_cor <-NA
est[is.na(est$Aim_1_cor),]$Aim_1_cfi <-NA
est[is.na(est$Aim_1_cor),]$Aim_1_tli <-NA
est[is.na(est$Aim_1_cor),]$Aim_1_rmsea <-NA
est[is.na(est$Aim_1_cor),]$Aim_1_srmr <-NA

#Remove CFA conditions for RA = 1
est[est$RA==1,]$Aim_1_cor <-NA
est[est$RA==1,]$Aim_1_cfi <-NA
est[est$RA==1,]$Aim_1_tli <-NA
est[est$RA==1,]$Aim_1_rmsea <-NA
est[est$RA==1,]$Aim_1_srmr <-NA

#Some quality of life operations
est$PCC <-  est$Avg - est$Avg_filtered_latents
est$CFA <- est$Aim_1_cor -  est$Avg_filtered_latents 
method_cols <- c('PCC', 'CFA')

#### ####


#### Collect the summary data ####
#Here we choose the directory containing all the simulation results.
#The prompt will ask you to choose a file. Choose any file in the folder
#containing the simulation results. For example, "results-row-1.rds"
#The directory containing the sim results will then be isolated.
dir_name <- file.choose()
dir_name<-dirname(dir_name)

#Resummarize the results, in case they are not stored in your local environment
summ<- reSummarise(Resumm_exclude_CFA_cor_out_of_range,dir_name)

#Relabel columns
names(summ)[2]<-'MEM'
names(summ)[3] <- 'RA'
names(summ)[7]<-'PCC'
names(summ)[19]<-'CFA'
names(summ)[28]<-'RMSEA'
names(summ)[31]<-'SRMR'

#Treat variables as factors
summ$N <- as.factor(summ$N)
summ$MEM <- as.factor(summ$MEM)
summ$RA <- as.factor(summ$RA)
summ$missing_rate <- as.factor(summ$missing_rate)

#Throw out single-assessment CFA summary data
summ[summ$RA==1,grep("cfa", names(summ))]<-NA

#Quality-of-life operation
SE_method_cols <- c('PCC', 'CFA')
#### ####




#FIGURE 4 in the manuscript
#Empirical bias estimates for 7-assessment cases and with no missing data

assessments <- 7  
wide_df <- est %>%
  subset(RA==assessments) %>%
  subset(missing_rate==0) %>%
  select(N, MEM, RA, missing_method, missing_rate, all_of(method_cols))
long_df <- melt(setDT(wide_df), 
                measure = method_cols,
                variable.name = 'method',
                value.name = 'Bias_value')



fig <- ggplot(long_df, aes(x = N, y = Bias_value, fill = method, pattern = MEM)) +
  geom_boxplot(aes(fill = method)) + 
  geom_hline(aes(yintercept=0))+
  scale_fill_manual(name = "method", values = c("#FDECCD", "#BAE4B3", "#329A55", "#075507")) +
  geom_boxplot_pattern(
    position = position_dodge(preserve = "single"),
    color = "black",
    pattern_fill = "white",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.025,
    pattern_key_scale_factor = 0.8
  ) +
  scale_pattern_manual(values=c('none', 'stripe', 'crosshatch', 'circle')) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), 
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  ylim(-0.5,0.5) +
  labs(
    y = 'Empirical Bias',
    title = paste0('Empirical bias estimates for 7-assessment cases and with no missing data'))

print(fig)



#### FIGURE 5 in the manuscript ####
#Estimates for empirical standard error by method

wide_df <- summ %>%
  select(N, MEM, RA, missing_method, missing_rate, all_of(SE_method_cols))
long_df <- melt(setDT(wide_df), 
                measure = SE_method_cols,
                variable.name = 'method',
                value.name = 'EmpSE_value')


fig <- ggplot(wide_df, aes(x=PCC,y=CFA,colour=N)) +
  geom_point() + 
  geom_abline(slope=1,intercept=0) +
  xlim(0,0.4) + 
  ylim(0,0.4) +
  labs(
    x = 'PCC', 
    y = 'CFA', 
    title = paste0('Estimates for empirical standard error by method'))
fig


#### ####


#### FIGURE 4A ####
#Distribution of empirical standard error when missing data rate is 0.40

missing_pattern <- 0.40
wide_df <- summ %>%
  subset(missing_rate==missing_pattern) %>%
  select(N, MEM, RA, missing_method, missing_rate, all_of(SE_method_cols))
long_df <- melt(setDT(wide_df), 
                measure = SE_method_cols,
                variable.name = 'method',
                value.name = 'SE_value')

fig <- ggplot(long_df, aes(x = N, y = SE_value, fill = method, pattern = MEM)) +
  geom_boxplot(aes(fill = method)) + 
  scale_fill_manual(name = "method", values = c("#FDECCD", "#BAE4B3")) + 
  geom_boxplot_pattern(
    position = position_dodge(preserve = "single"), 
    color = "black", 
    pattern_fill = "white", 
    pattern_angle = 45, 
    pattern_density = 0.1, 
    pattern_spacing = 0.025, 
    pattern_key_scale_factor = 0.8) +
  scale_pattern_manual(values=c('none', 'stripe', 'crosshatch', 'circle')) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), 
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  labs(
    y = 'Empirical Standard Error',
    title = paste0('Distribution of empirical standard error when missing data rate is 0.40')) 


print(fig)


#Fix missing rate, display by method, sample size and measurement error magnitude 
#### FIGURE 4B in the manuscript ####
#Distribution of empirical bias when missing rate is 0.40
missing_pattern <- 0.4
wide_df <- est %>%
  subset(missing_rate==missing_pattern) %>%
  select(N, MEM, RA, missing_method, missing_rate, all_of(method_cols))
long_df <- melt(setDT(wide_df), 
                measure = method_cols,
                variable.name = 'method',
                value.name = 'Bias_value')

fig <- ggplot(long_df, aes(x = N, y = Bias_value, fill = method, pattern = MEM)) +
  geom_boxplot(aes(fill = method)) + 
  geom_hline(yintercept = 0) +
  scale_fill_manual(name = "method", values = c("#FDECCD", "#BAE4B3")) + 
  geom_boxplot_pattern(
    position = position_dodge(preserve = "single"), 
    color = "black", 
    pattern_fill = "white", 
    pattern_angle = 45, 
    pattern_density = 0.1, 
    pattern_spacing = 0.025, 
    pattern_key_scale_factor = 0.8) +
  scale_pattern_manual(values=c('none', 'stripe', 'crosshatch', 'circle')) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), 
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  labs(
    y = 'Empirical Bias',
    title = paste0('Distribution of empirical bias when missing rate is 0.40'))

print(fig)




#### FIGURE 5a in the manuscript####
#Fix sample size, display by method, measurement error magnitude
#and number of assessments (exclude single assessment) 
#Distribution of empirical standard error when sample size is 35

samples <- 35
wide_df <- summ %>%
  subset(RA!=1) %>%
  subset(N==samples) %>%
  select(N, MEM, RA, missing_method, missing_rate, all_of(SE_method_cols))
long_df <- melt(setDT(wide_df), 
                measure = SE_method_cols,
                variable.name = 'method',
                value.name = 'SE_value')


fig <- ggplot(long_df, aes(x = MEM, y = SE_value, fill = method, pattern = RA)) +
  geom_boxplot(aes(fill = method)) + 
  scale_fill_manual(name = "method", values = c("#FDECCD", "#BAE4B3")) + 
  geom_boxplot_pattern(
    position = position_dodge(preserve = "single"), 
    color = "black", 
    pattern_fill = "white", 
    pattern_angle = 45, 
    pattern_density = 0.1, 
    pattern_spacing = 0.025, 
    pattern_key_scale_factor = 0.8) +
  scale_pattern_manual(values=c('none', 'stripe', 'crosshatch', 'circle')) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), 
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  labs(
    y = 'Empirical Standard Error',
    title = paste0('Distribution of empirical standard error when sample size is 35')) 

print(fig)



#Distribution of empirical bias when sample size is 35 or 100
for (i in c(35,100)) {
  
  samples <- i 
  wide_df <- est %>%
    subset(RA!=1) %>%
    subset(N==samples) %>%
    select(N, MEM, RA, missing_method, missing_rate, all_of(method_cols))
  long_df <- melt(setDT(wide_df), 
                  measure = method_cols,
                  variable.name = 'method',
                  value.name = 'Bias_value')
  
  
  fig <- ggplot(long_df, aes(x = MEM, y = Bias_value, fill = method, pattern = RA)) +
    geom_boxplot(aes(fill = method)) + 
    geom_hline(yintercept = 0) +
    scale_fill_manual(name = "method", values = c("#FDECCD", "#BAE4B3")) + 
    geom_boxplot_pattern(
      position = position_dodge(preserve = "single"), 
      color = "black", 
      pattern_fill = "white", 
      pattern_angle = 45, 
      pattern_density = 0.1, 
      pattern_spacing = 0.025, 
      pattern_key_scale_factor = 0.8) +
    scale_pattern_manual(values=c('none', 'stripe', 'crosshatch', 'circle')) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")), 
           fill = guide_legend(override.aes = list(pattern = "none"))) +
    labs(
      y = 'Empirical Bias',
      title = paste0('Distribution of empirical bias when sample size is ', samples)) 
 
  print(fig)
  
}


#### FIGURE 7 ####
#Rate of acceptable fit when N=100
#By measurement error

method_cols <- c('RMSEA','SRMR')
wide_df <- summ %>%
  subset(N==100) %>%
  select(N, MEM, RA, missing_method, missing_rate, all_of(method_cols))
long_df <- melt(setDT(wide_df), 
                measure = method_cols,
                variable.name = 'method',
                value.name = 'CFA_Fit_value')

fig <- ggplot(long_df, aes(x = MEM, y = CFA_Fit_value, fill = method)) +
  geom_boxplot(aes(fill = method)) +
  #geom_violin(aes(fill = method)) +
  scale_fill_manual(name = "method", values = c("#BAE4B3","#329A55")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), 
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  labs(
    y = 'Rate of acceptable fit',
    title = paste0('Rate of acceptable fit when N=100')) #, by measurement error magnitude, '))

print(fig)


#### FIGURE Supp2a ####
#Fix sample size, display by method, measurement error magnitude
#and number of assessments (exclude single assessment) 
#Distribution of empirical standard error when sample size is 100
samples <- 100
wide_df <- summ %>%
  subset(RA!=1) %>%
  subset(N==samples) %>%
  select(N, MEM, RA, missing_method, missing_rate, all_of(SE_method_cols))
long_df <- melt(setDT(wide_df), 
                measure = SE_method_cols,
                variable.name = 'method',
                value.name = 'SE_value')


fig <- ggplot(long_df, aes(x = MEM, y = SE_value, fill = method, pattern = RA)) +
  geom_boxplot(aes(fill = method)) + 
  scale_fill_manual(name = "method", values = c("#FDECCD", "#BAE4B3")) + 
  geom_boxplot_pattern(
    position = position_dodge(preserve = "single"), 
    color = "black", 
    pattern_fill = "white", 
    pattern_angle = 45, 
    pattern_density = 0.1, 
    pattern_spacing = 0.025, 
    pattern_key_scale_factor = 0.8) +
  scale_pattern_manual(values=c('none', 'stripe', 'crosshatch', 'circle')) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), 
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  labs(
    y = 'Empirical Standard Error',
    title = paste0('Distribution of empirical standard error when sample size is 100')) #, samples))


print(fig)


#### Supp 1a ####
#Grouped by number of repeated assessments and N, exclude single assessment case
#Empirical SE grouped by number of repeated assessments and sample size
wide_df <- summ %>%
  subset(RA!=1) %>%
  select(N, MEM, RA, missing_method, missing_rate, all_of(SE_method_cols))
long_df <- melt(setDT(wide_df), 
                measure = SE_method_cols,
                variable.name = 'method',
                value.name = 'EmpSE_value')


fig <- ggplot(long_df, aes(x = RA, y = EmpSE_value, fill = method, pattern = N)) +
  geom_boxplot(aes(fill = method)) + 
  scale_fill_manual(name = "method", values = c("#FDECCD", "#BAE4B3")) + 
  geom_boxplot_pattern(
    position = position_dodge(preserve = "single"), 
    color = "black", 
    pattern_fill = "white", 
    pattern_angle = 45, 
    pattern_density = 0.1, 
    pattern_spacing = 0.025, 
    pattern_key_scale_factor = 0.8) +
  scale_pattern_manual(values=c('none', 'stripe', 'crosshatch', 'circle')) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), 
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  labs(
    y = 'Empirical Standard Error',
    title = paste0('Empirical SE grouped by number of repeated assessments and sample size')) 


print(fig)


#### Supp 1b ####
#Grouped by repeated assessments and measurement error magnitude,  exclude single assessment case
#Empirical SE grouped by number of repeated assessments and MEM')) 

wide_df <- summ %>%
  subset(RA!=1) %>%
  select(N, MEM, RA, missing_method, missing_rate, all_of(SE_method_cols))
long_df <- melt(setDT(wide_df), 
                measure = SE_method_cols,
                variable.name = 'method',
                value.name = 'EmpSE_value')


fig <- ggplot(long_df, aes(x = RA, y = EmpSE_value, fill = method, pattern = MEM)) +
  geom_boxplot(aes(fill = method)) + 
  scale_fill_manual(name = "method", values = c("#FDECCD", "#BAE4B3")) + 
  geom_boxplot_pattern(
    position = position_dodge(preserve = "single"), 
    color = "black", 
    pattern_fill = "white", 
    pattern_angle = 45, 
    pattern_density = 0.1, 
    pattern_spacing = 0.025, 
    pattern_key_scale_factor = 0.8) +
  scale_pattern_manual(values=c('none', 'stripe', 'crosshatch', 'circle')) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")), 
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  labs(
    y = 'Empirical Standard Error',
    title = paste0('Empirical SE grouped by number of repeated assessments and MEM')) 


print(fig)

