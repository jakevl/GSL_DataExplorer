
makeMetalsSummaryTable=function(data,years_metals,surf_deep_metals,matrix,parameter,aggby){
	#Subset by year input
	metalsData_table=data[data$Year>=years_metals[1]&data$Year<=years_metals[2],]
	metalsData_table=data[data$Year>=years_metals[1]&data$Year<=years_metals[2],]
	if(matrix==1){
		#Subset by surface/deep input
		surf_deep_metals[surf_deep_metals=="Surface"]="S"
		surf_deep_metals[surf_deep_metals=="Deep"]="D"
		surf_deep_metals=append(surf_deep_metals,rep("NA",2-length(surf_deep_metals)))
		metalsData_table=metalsData_table[metalsData_table$WaterColLoc==surf_deep_metals[1]|metalsData_table$WaterColLoc==surf_deep_metals[2],]}
	if(matrix==2){
		metalsData_table=metalsData_table[metalsData_table$Matrix=="Brine shrimp",]}
	if(matrix==3){
		metalsData_table=metalsData_table[metalsData_table$Matrix=="Bird eggs",]}
	metalsData_table=metalsData_table[metalsData_table$Parameter==parameter,]
	metalsData_table=metalsData_table[,c("Result",aggby)]
	metalsData_table=na.omit(metalsData_table)
	if(dim(metalsData_table)[1]>0){
		min=aggregate(metalsData_table[,"Result"]~metalsData_table[,aggby],FUN='min')
		median=aggregate(metalsData_table[,"Result"]~metalsData_table[,aggby],FUN='median')
		mean=aggregate(metalsData_table[,"Result"]~metalsData_table[,aggby],FUN='mean')
		max=aggregate(metalsData_table[,"Result"]~metalsData_table[,aggby],FUN='max')
		n=aggregate(metalsData_table[,"Result"]~metalsData_table[,aggby],FUN='length')
		minT=min(metalsData_table[,"Result"])
		medianT=median(metalsData_table[,"Result"])
		meanT=mean(metalsData_table[,"Result"])
		maxT=max(metalsData_table[,"Result"])
		nT=length(metalsData_table[,"Result"])
		stats=cbind(min,median[,2],mean[,2],max[,2],n[,2])
		Total=c(minT,medianT,meanT,maxT,nT)
		columnNames=stats[,1]
		stats=t(stats[,2:6])
		colnames(stats)=columnNames
		stats=cbind(stats,Total)
		rownames(stats)=c("Min","Median","Mean","Max","Count")}
	else{stats=""}
	return(stats)
	}
