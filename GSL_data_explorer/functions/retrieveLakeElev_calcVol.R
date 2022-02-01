
#Pull lake elevation data and calculate UBL & DBL volumes

options(repos = c(CRAN = "https://cran.rstudio.com"))
setwd("F:\\Shiny\\GSL_data_explorer\\GSL_data_explorer")
library(waterData)

#Query south arm elevation. Calculate estimated BL volumes.
#Input codes and matching names
par_codes="62614"
par_names="LakeElevFt"

####Run code:
parameters=cbind(par_codes,par_names)
row.names(parameters)=parameters[,"par_codes"]

elev_data=importDVs(staid="10010000",code="62614",stat="00003",sdate="2010-01-01")
colnames(elev_data)=c("staid","LakeElevFt","Date","qualcode")

elev_data$LakeElev_m=elev_data$LakeElevFt/3.28084
elev_data$DBL_elev_N=-46.51008+1.03163*elev_data$LakeElev_m
elev_data$DBL_elev_N[elev_data$DBL_elev_N<1269.43714800001]=1269.43714800001 #Setting minimum possible N DBL elevation at inflection point of formula
elev_data$DBL_elev_S=-104.967+1.07643*elev_data$LakeElev_m
elev_data$DBL_elev_S[elev_data$DBL_elev_S<1270.92582600001]=1270.92582600001 #Setting minimum possible S DBL elevation at inflection point of formula
elev_data$DBL_vol_N_m3=59633174862208.5-93949911022.036*elev_data$DBL_elev_N+37003662.2562751*elev_data$DBL_elev_N^2
elev_data$DBL_vol_S_m3=120180569424473-189400877631.774*elev_data$DBL_elev_S+74622375.23038100*elev_data$DBL_elev_S^2
elev_data$DBL_vol_total_m3=elev_data$DBL_vol_S_m3+elev_data$DBL_vol_N_m3
elev_data$UBL_vol_m3=(137054926141519-216113402592.612*elev_data$LakeElev_m+85193657.0902014*elev_data$LakeElev_m^2)-elev_data$DBL_vol_total_m3
elev_data$GB_vol_m3=137054926141519-216113402592.612*elev_data$LakeElev_m+85193657.0902014*elev_data$LakeElev_m^2
elev_data$DBL_vol_acreft=elev_data$DBL_vol_total_m3*0.000810714
elev_data$UBL_vol_acreft=elev_data$UBL_vol_m3*0.000810714
elev_data$GB_vol_acreft=elev_data$GB_vol_m3*0.000810714

write.csv(file="data/UBL_DBL_vol.csv",elev_data)




######Long term GSL elevation data and volume estimates data pull:

#Query south arm elevation. Calculate estimated BL volumes.
#Input codes and matching names
par_codes="62614"
par_names="LakeElevFt"

####Run code:
parameters=cbind(par_codes,par_names)
row.names(parameters)=parameters[,"par_codes"]

elev_data=importDVs(staid="10010000",code="62614",stat="00003")
colnames(elev_data)=c("staid","LakeElevFt","Date","qualcode")

elev_data$LakeElev_m=elev_data$LakeElevFt/3.28084
elev_data$DBL_elev_N=-46.51008+1.03163*elev_data$LakeElev_m
elev_data$DBL_elev_N[elev_data$DBL_elev_N<1269.43714800001]=1269.43714800001 #Setting minimum possible N DBL elevation at inflection point of formula
elev_data$DBL_elev_S=-104.967+1.07643*elev_data$LakeElev_m
elev_data$DBL_elev_S[elev_data$DBL_elev_S<1270.92582600001]=1270.92582600001 #Setting minimum possible S DBL elevation at inflection point of formula
elev_data$DBL_vol_N_m3=59633174862208.5-93949911022.036*elev_data$DBL_elev_N+37003662.2562751*elev_data$DBL_elev_N^2
elev_data$DBL_vol_S_m3=120180569424473-189400877631.774*elev_data$DBL_elev_S+74622375.23038100*elev_data$DBL_elev_S^2
elev_data$DBL_vol_total_m3=elev_data$DBL_vol_S_m3+elev_data$DBL_vol_N_m3
elev_data$UBL_vol_m3=(137054926141519-216113402592.612*elev_data$LakeElev_m+85193657.0902014*elev_data$LakeElev_m^2)-elev_data$DBL_vol_total_m3
elev_data$GB_vol_m3=137054926141519-216113402592.612*elev_data$LakeElev_m+85193657.0902014*elev_data$LakeElev_m^2
elev_data$DBL_vol_acreft=elev_data$DBL_vol_total_m3*0.000810714
elev_data$UBL_vol_acreft=elev_data$UBL_vol_m3*0.000810714
elev_data$GB_vol_acreft=elev_data$GB_vol_m3*0.000810714


###Query north arm elevation. Calculate estimated NA volume.

#Input codes and matching names
par_codes="62614"
par_names="LakeElevFt"

####Run code:
parameters=cbind(par_codes,par_names)
row.names(parameters)=parameters[,"par_codes"]

NA_elev_data=importDVs(staid="10010100",code="62614",stat="00003")
colnames(NA_elev_data)=c("staid","NA_LakeElevFt","Date","qualcode")

NA_elev_data$NA_vol_acreft=5013.32379115*NA_elev_data$NA_LakeElevFt^2-41775852.41544460*NA_elev_data$NA_LakeElevFt+87029164222.46370000
NA_elev_data=NA_elev_data[,c("NA_LakeElevFt","Date","NA_vol_acreft")]

#Append to SA calculations in new data frame (first manipulate SA columns)
SA_elev_data=elev_data[,c("LakeElevFt","Date","GB_vol_acreft")]
names(SA_elev_data)=c("SA_LakeElevFt","Date","SA_vol_acreft")

combined_elev_data=merge(SA_elev_data,NA_elev_data)
combined_elev_data$SmNhead=combined_elev_data$SA_LakeElevFt-combined_elev_data$NA_LakeElevFt

#Flatten
library(reshape2)
combined_elev_data_flat=melt(combined_elev_data,id.vars="Date")

#Write to data folder
write.csv(file="data/combined_elev.csv",row.names=FALSE,combined_elev_data_flat)





