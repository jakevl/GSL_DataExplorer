####Pre-processing for nutrient mass estimates

#setwd("F:\\Shiny\\GSL_data_explorer")
vol_data=read.csv(file="data/UBL_DBL_vol.csv")
wq_data=read.csv("data/GSL_WQ.csv")

vol_data$Date=as.Date(vol_data$Date)
vol_data$Year=year(vol_data$Date)
vol_data$Month=month(vol_data$Date)
wq_data$Date=as.Date(wq_data$sample_dt)
wq_data$Year=year(wq_data$Date)
wq_data$Month=month(wq_data$Date)
wq_data=wq_data[wq_data$site_name!="BRB11"&wq_data$site_name!="FB9"&wq_data$site_name!="FB10"&wq_data$site_name!="GB8",]

DBL_vol_YrMo=aggregate(DBL_vol_total_m3~Year+Month,FUN='mean',vol_data)
UBL_vol_YrMo=aggregate(UBL_vol_m3~Year+Month,FUN='mean',vol_data)
vol_YrMo=merge(UBL_vol_YrMo,DBL_vol_YrMo,by=c("Year","Month"))

wq_data$bl=vector(length=dim(wq_data)[1])
for(n in 1:dim(wq_data)[1]){
	if((wq_data[n,"site_name"]=="GB2"|wq_data[n,"site_name"]=="GB5")&wq_data[n,"Sampling_depth_m"]>=5){
		wq_data[n,"bl"]="DBL"}
	if(wq_data[n,"Sampling_depth_m"]<=2){wq_data[n,"bl"]="UBL"}
	}

wq_data$Salinity_mgL=wq_data$Salinity_gL*1000

TP_UBL=aggregate(Total_phosphorus_mgL~Year+Month,FUN='mean',wq_data[wq_data$bl=="UBL",])
FP_UBL=aggregate(Filtered_phosphorus_mgL~Year+Month,FUN='mean',wq_data[wq_data$bl=="UBL",])
TN_UBL=aggregate(Total_nitrogen_mgL~Year+Month,FUN='mean',wq_data[wq_data$bl=="UBL",])
FN_UBL=aggregate(Filtered_nitrogen_mgL~Year+Month,FUN='mean',wq_data[wq_data$bl=="UBL",])
salt_UBL=aggregate(Salinity_mgL~Year+Month,FUN='mean',wq_data[wq_data$bl=="UBL",])

TP_DBL=aggregate(Total_phosphorus_mgL~Year+Month,FUN='mean',wq_data[wq_data$bl=="DBL",])
FP_DBL=aggregate(Filtered_phosphorus_mgL~Year+Month,FUN='mean',wq_data[wq_data$bl=="DBL",])
TN_DBL=aggregate(Total_nitrogen_mgL~Year+Month,FUN='mean',wq_data[wq_data$bl=="DBL",])
FN_DBL=aggregate(Filtered_nitrogen_mgL~Year+Month,FUN='mean',wq_data[wq_data$bl=="DBL",])
salt_DBL=aggregate(Salinity_mgL~Year+Month,FUN='mean',wq_data[wq_data$bl=="DBL",])

colnames(TP_UBL)=c("Year","Month","TP_UBL")
colnames(FP_UBL)=c("Year","Month","FP_UBL")
colnames(TN_UBL)=c("Year","Month","TN_UBL")
colnames(FN_UBL)=c("Year","Month","FN_UBL")
colnames(salt_UBL)=c("Year","Month","salt_UBL")
colnames(TP_DBL)=c("Year","Month","TP_DBL")
colnames(FP_DBL)=c("Year","Month","FP_DBL")
colnames(TN_DBL)=c("Year","Month","TN_DBL")
colnames(FN_DBL)=c("Year","Month","FN_DBL")
colnames(salt_DBL)=c("Year","Month","salt_DBL")

mass_data=merge(vol_YrMo,TP_UBL,by=c("Year","Month"),all.x=T)
mass_data=merge(mass_data,FP_UBL,by=c("Year","Month"),all.x=T)
mass_data=merge(mass_data,TN_UBL,by=c("Year","Month"),all.x=T)
mass_data=merge(mass_data,FN_UBL,by=c("Year","Month"),all.x=T)
mass_data=merge(mass_data,salt_UBL,by=c("Year","Month"),all.x=T)
mass_data=merge(mass_data,TP_DBL,by=c("Year","Month"),all.x=T)
mass_data=merge(mass_data,FP_DBL,by=c("Year","Month"),all.x=T)
mass_data=merge(mass_data,TN_DBL,by=c("Year","Month"),all.x=T)
mass_data=merge(mass_data,FN_DBL,by=c("Year","Month"),all.x=T)
mass_data=merge(mass_data,salt_DBL,by=c("Year","Month"),all.x=T)
mass_data$YrMo=as.Date(paste0(mass_data$Year,"-",mass_data$Month,"-01"),format='%Y-%m-%d')
mass_data=mass_data[order(mass_data$YrMo),]

#ublnames=c("TP_UBL","FP_UBL","TN_UBL","FN_UBL")
#dblnames=c("TP_DBL","FP_DBL","TN_DBL","FN_DBL")
#
#pdf(file="GB_nut_mass.pdf",height=11,width=8.5)
#par(mfrow=c(4,1),mar=c(3,4.1,3,4))
#for(n in 1:length(ublnames)){
#	ublmass_n=mass_data[,ublnames[n]]*mass_data[,"UBL_vol_m3"]*0.001
#	dblmass_n=mass_data[,dblnames[n]]*mass_data[,"DBL_vol_total_m3"]*0.001
#	totmass_n=ublmass_n+dblmass_n
#	label=paste0(substr(ublnames[n], start = 1, stop = 2)," (kg)")
#	ymin_n=min(ublmass_n,dblmass_n,totmass_n)
#	ymax_n=max(ublmass_n,dblmass_n,totmass_n)
#	plot(totmass_n~mass_data$YrMo,log='y',ylim=c(ymin_n,ymax_n),type='b',col="blue",pch=21,ylab=label,xlab="")
#	points(ublmass_n~mass_data$YrMo,type='b',col="orange",pch=19)
#	points(dblmass_n~mass_data$YrMo,type='b',col="green",pch=17)
#	legend("bottomleft",bty='n',pch=c(21,19,17),col=c("blue","orange","green"),lty=1,legend=c("Total","UBL","DBL"))
#	plotmeans(totmass_n~mass_data$Year,ylab=label,main="Gilbert Bay total pool",xlab="")
#	}
#dev.off()

massPlot=function(param,mass_data,label){
	ublmass=mass_data[,paste0(param,"_UBL")]*mass_data[,"UBL_vol_m3"]*0.001
	dblmass=mass_data[,paste0(param,"_DBL")]*mass_data[,"DBL_vol_total_m3"]*0.001
	totmass=ublmass+dblmass
	ymin=min(ublmass,dblmass,totmass,na.rm=T)
	ymax=max(ublmass,dblmass,totmass,na.rm=T)
	plot(totmass~mass_data$YrMo,log='y',ylim=c(ymin,ymax),pch=NA,ylab=label,xlab="",cex=2,cex.lab=1.5,cex.axis=1.5)
	data1=na.omit(data.frame(totmass,as.Date(mass_data$YrMo)))
	colnames(data1)=c("totmass","YrMo")
	points(totmass~YrMo,data=data1,ylim=c(ymin,ymax),type='b',col="blue",pch=21,cex=2)
	data2=na.omit(data.frame(ublmass,as.Date(mass_data$YrMo)))
	colnames(data2)=c("ublmass","YrMo")
	points(ublmass~YrMo,data=data2,type='b',col="orange",pch=19,cex=2)
	data3=na.omit(data.frame(dblmass,as.Date(mass_data$YrMo)))
	colnames(data3)=c("dblmass","YrMo")
	points(dblmass~YrMo,data=data3,type='b',col="green",pch=17,cex=2)
	legend("bottomleft",bty='n',pch=c(21,19,17),col=c("blue","orange","green"),lty=1,legend=c("Total","UBL","DBL"),cex=2)
	massTable=data.frame(ublmass,dblmass,totmass,mass_data$YrMo)
	#return(massTable)
}





