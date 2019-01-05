###GSL water quality retrieval
options(repos = c(CRAN = "https://cran.rstudio.com"))
setwd("F:\\Shiny\\GSL_data_explorer")
library(dataRetrieval)

#select desired sites
sites=c(410224112095101,
		410401112134801,
		404607112193801,
		405356112205601,
		410422112200001,
		410228112245301,
		410323112301901,
		410637112270401,
		410644112382601,
		411116112244401,
		411403112200801,
		411218112414401,
		410153112082301)

#Input codes and matching names
par_codes=c("00010",
			"00078",
			"00098",			
			"00600",
			"00665",
			"00666",
			#"00680",
			#"00681",
			"32209",
			"62854",
			"62855",
			"70305",
			"70953",
			"72263")
			
par_names=c("Water_temp_C",
			"Secchi_depth_m",
			"Sampling_depth_m",
			"Total_nitrogen_mgL",
			"Total_phosphorus_mgL",
			"Filtered_phosphorus_mgL",
			#"OC_mgL",
			#"DOC_mgL",
			"Chlorophyll_a_ugL",
			"Filtered_nitrogen_mgL",
			"Total_nitrogen_mgL",
			"Salinity_gL",
			"Chlorophyll_a_ugL",
			"Density")

####Run code:
parameters=cbind(par_codes,par_names)
row.names(parameters)=parameters[,"par_codes"]

data_raw=readNWISqw(sites,parameters[,1])




library(reshape)
data_mat=cast(data_raw,site_no+sample_dt+sample_tm~parm_cd,value="result_va")
names=colnames(data_mat)[4:(dim(parameters)[1]+3)]
matNames=as.vector(c("site_no","sample_dt","sample_tm",c(parameters[,"par_codes"])))
colnames(data_mat)=matNames
data_mat=data.frame(data_mat)

#Merge salinity cols, chla, and TN cols.
for(n in 1:dim(data_mat)[1]){
	#if(is.na(data_mat$X70305[n])){data_mat$X70305[n]=data_mat$X00480[n]}
	if(is.na(data_mat$X32209[n])){data_mat$X32209[n]=data_mat$X70953[n]}
	if(is.na(data_mat$X62855[n])){data_mat$X62855[n]=data_mat$X00600[n]}
	}
	
#Delete old salinity, chla, and TN cols
keep=!colnames(data_mat) %in% c("X00480","X70953","X00600")
data_mat=data_mat[,keep]

#Convert water temp > 100 C to NA
data_mat$X00010[data_mat$X00010>=100]=NA

#drop NA depth rows
data_mat=data_mat[!is.na(data_mat$X00098),]

#Rename columns
row.names(parameters)=paste0("X",row.names(parameters))
parameters=parameters[colnames(data_mat)[4:dim(data_mat)[2]],]
names=c(colnames(data_mat)[1:3],parameters[,"par_names"])
colnames(data_mat)=names
#head(data_mat)

#Add site names
site_trans=read.csv(file="data/GSL_USGS_SiteNoSiteName.csv")
site_trans=site_trans[,c("site_no","Baseline_site_name")]
data_mat=merge(data_mat,site_trans,by="site_no")
names=c(colnames(data_mat)[1:3],"Baseline_site_name",parameters[,"par_names"])
data_mat=data_mat[,names]
colnames(data_mat)=c(colnames(data_mat)[1:3],"site_name",parameters[,"par_names"])

#Divide salinity (g/L) by (density(g/cm3)*10) for % salinity, convert salinity >=80% to NA, convert salinity = 0% to 0.1% for log plots
data_mat$Salinity_pct=data_mat$Salinity_gL/(data_mat$Density*10)
data_mat$Salinity_pct[data_mat$Salinity_pct>=80]=NA
data_mat$Salinity_pct[data_mat$Salinity_pct==0]=0.1
data_mat$Salinity_gL[data_mat$Salinity_gL>=800]=NA
data_mat$Salinity_gL[data_mat$Salinity_gL==0]=1

#Removing TP value > 9 mg/L
data_mat$Total_phosphorus_mgL[data_mat$Total_phosphorus_mgL>=9]=NA

write.csv(file="data/GSL_WQ.csv",data_mat,row.names=FALSE)

source("GSL_Data_Explorer/functions/retrieveLakeElev_calcVol.R")


library(rsconnect)
deployApp("F:\\Shiny\\GSL_data_explorer",account="udwq")
#deployApp("F:\\Shiny\\GSL_data_explorer",account="jakevl")

library(shiny)
#runApp("F:\\Shiny\\GSL_data_explorer")

