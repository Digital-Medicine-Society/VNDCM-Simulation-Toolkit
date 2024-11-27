library(ggplot2)
library(tidySEM)

#Choose data set
data_set <- Final_Results_Run_5
  #Absolutely_Final_no_CFA_2
  #Absolutely_Final_wCRO_alpha_07

# dir_name <- file.choose()
# dir_name<-dirname(dir_name)

#Import estimates from first condition
final_sim_res_1<- SimResults(data_set,1)
  #readRDS(paste(dir_name,"/results-row-",1,".rds",sep=""))                


estimates_data_1<- with(final_sim_res_1,
                       data.frame(N=rep(condition$N,nrow(results)),
                                  meas_error_mag = rep(condition$meas_error_mag,nrow(results)),
                                  n_assess=rep(condition$n_assess,nrow(results)),
                                  missing_method=rep(condition$missing_method,nrow(results)),
                                  missing_rate = rep(condition$missing_rate,nrow(results)),results)
                       
)

#grand estimate data.frame
collated_estimates <- estimates_data_1

#loop through remaining conditions, adding estimates to the grand estimate data.frame                          
for (i in 2:nrow(data_set)) {
  
  nam <- paste("final_sim_res_",i,sep="")
  
  #assign(nam,readRDS(paste(dir_name,"/results-row-",i,".rds",sep="")))
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
  
  #Keep the workspace tidy
  rm(list = ls(pattern = paste("final_sim_res_",i,sep="")))  
  rm(list = ls(pattern = nam))
  
}

dir_name <- file.choose()
dir_name<-dirname(dir_name)

write.csv(collated_estimates, paste(dir_name,"/Collated Estimates.csv",sep=""))

write.csv(data_set,paste(dir_name,"/Summarised results.csv",sep=""))
