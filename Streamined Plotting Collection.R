library(bigrquery)
library(ggplot2)
library(ggpattern)
library(dplyr)
library(data.table)

options(repr.plot.width=12, repr.plot.height=8)


## Estimates ##

#### Collect the estimates ####
est <- read.csv(file.choose())[2:21] #collated_estimates

names(est)[1]<-"N"
names(est)[2]<-"MEM"
names(est)[3]<-'RA'
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

est$PCC <-  est$Avg - est$Avg_filtered_latents
est$CFA <- est$Aim_1_cor -  est$Avg_filtered_latents 
method_cols <- c('PCC', 'CFA')

#### ####


#### Collect the summary data ####
#summ <- read.csv(file.choose())[2:43]
dir_name <- file.choose()
dir_name<-dirname(dir_name)
summ<- reSummarise(Resumm_exclude_CFA_cor_out_of_range,dir_name)
#Final_Results_Run_5
#reSummarise(Summarise_with_latents_ignore_large_CFA_vals, dir_name,Design = Design4)
#Absolutely_Final_wCRO_alpha_07
#Absolutely_Final_no_CFA_2

names(summ)[2]<-'MEM'
names(summ)[3] <- 'RA'
names(summ)[7]<-'PCC'
names(summ)[19]<-'CFA'


summ$N <- as.factor(summ$N)
summ$MEM <- as.factor(summ$MEM)
summ$RA <- as.factor(summ$RA)
summ$missing_rate <- as.factor(summ$missing_rate)

#Throw out single-assessment CFA summary data
summ[summ$RA==1,grep("cfa", names(summ))]<-NA


SE_method_cols <- c('PCC', 'CFA')
#### ####




#FIGURE 4 in the manuscript
#Split by MEM and repeated assessments, missing rate fixed at 0
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
    title = paste0('Empirical bias estimates for 7-assessment cases and with no missing data'))#, split by sample size and measurement error')) #+



print(fig)

jpeg("Figure 4.jpeg", units="in", width=10, height=8, res=300) 
fig
dev.off()


for (samp in c(35,100,200,1000)) {  
  print(samp)  
  for (MEM in c(0.5,1,1.5,2)) {  
    wide_df_sub<-est %>%
      subset(N==samp) %>%
      subset(MEM==MEM) %>%
      subset(missing_rate==0)%>%
      subset(RA==1) %>%
      select(N, MEM, RA, missing_method, missing_rate, all_of(method_cols))  
    
    print(MEM)
    print(median(wide_df_sub$PCC))
    print(median(wide_df_sub$CFA))
    
  }
}

#### FIGURE 5 in the manuscript ####
#scatter plot
fig <- ggplot(wide_df, aes(x=Emp_SE_cor,y=cfa1_cor_SE,colour=N)) +
  geom_point() + 
  geom_abline(slope=1,intercept=0) +
  xlim(0,0.4) + 
  ylim(0,0.4) +
  labs(
    x = 'PCC', #'Pearson correlation',
    y = 'CFA', #'CFA factor correlation',
    title = paste0('Estimates for empSE by method')) #+
fig

jpeg("Figure 5.jpeg", units="in", width=10, height=8, res=300) 
fig
dev.off()

#### ####

