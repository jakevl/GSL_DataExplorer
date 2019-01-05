
#Task scheduler
#library(taskscheduleR)
myscript="F:\\Shiny\\GSL_data_explorer\\retrieveAndDeploy.R"
start.date="12/20/2017"
taskscheduler_create(taskname = "retrieveAndDeploy_GSL_trophic_explorer", rscript = myscript,schedule = "MONTHLY", starttime = "08:30", startdate=start.date,days='1')
#taskscheduler_delete("retrieveAndDeploy_GSL_trophic_explorer")



