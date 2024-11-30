library(ggplot2)
library(tidySEM)

####Collate the estimates from each simmulation into one grand data frame, 
####and write the estimates data frame and summary results to CSV.

data_set <- Final_Results

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

#User is prompted to choose a directory to save the outputs to.
dir_name <- file.choose()
dir_name<-dirname(dir_name)

#Write the grand data frame of estimates to CSV
write.csv(collated_estimates, paste(dir_name,"/Collated Estimates.csv",sep=""))

#Write the simulation performance measures summary to CSV.
write.csv(data_set,paste(dir_name,"/Summarised results.csv",sep=""))
