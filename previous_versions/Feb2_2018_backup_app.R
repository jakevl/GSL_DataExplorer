
#setwd("F:/Shiny/GSL_trophic_explorer")
gsl_data=read.csv(file="GSL_WQ.csv")
NLA_data=read.csv(file="NLA2012_TrophicData.csv")
NCCA_data=read.csv(file="assessed_ncca2010_waterchem.csv")

library(shiny)
library(lubridate)
library(reshape2)
colnames(gsl_data)[colnames(gsl_data)=="site_name"]="Site"
options(scipen=999)
gsl_data$Date=as.Date(gsl_data$sample_dt,format='%Y-%m-%d')
gsl_data$Year=year(gsl_data$Date)
gsl_data$Depth=rep(NA,length=dim(gsl_data)[1])
gsl_data$Depth[gsl_data$Sampling_depth_m<=1]="Surface"
gsl_data$Depth[gsl_data$Sampling_depth_m>=4]="Deep"
gsl_data$TNTP_mol=gsl_data$Total_nitrogen_mgL/gsl_data$Total_phosphorus_mgL*2.211323
compSP_choices=c("Total_nitrogen_mgL","Total_phosphorus_mgL","Chlorophyll_a_ugL","TNTP_mol")
names(compSP_choices)=c("Total nitrogen (mg/L)","Total phosphorus (mg/L)","Chlorophyll a (ug/L)","TN:TP (mol)")
param_names=c(
"Water_temp_C",
"Secchi_depth_m",
"Sampling_depth_m",
"Total_phosphorus_mgL",
"Filtered_phosphorus_mgL",
"Chlorophyll_a_ugL",
"Filtered_nitrogen_mgL",
"Total_nitrogen_mgL",
"Salinity_pct",
"TNTP_mol")
names(param_names)=c(
"Water temp (C)",
"Secchi depth (m)",
"Sampling depth (m)",
"Total phosphorus (mg/L)",
"Filtered phosphorus (mg/L)",
"Chlorophyll a (ug/L)",
"Filtered nitrogen (mg/L)",
"Total nitrogen (mg/L)",
"Salinity (%)",
"TN:TP (mol)")

mass_names=c("TP","FP","TN","FN","salt")
names(mass_names)=c("Total phosphorus","Filtered phosphorus","Total nitrogen","Filtered nitrogen","Salt")

NCCA_data_mat=dcast(NCCA_data,UID~PARAMETER_2,value.var="RESULT",)
NCCA_data_mat$Secchi_m=rep(NA,dim(NCCA_data_mat)[1])
NCCA_data_mat$TNTP_mol=NCCA_data_mat$Total_nitrogen_mgL/NCCA_data_mat$Total_phosphorus_mgL*2.211323

source("nut_mass_preprocess.R")

# Define UI for app ----
ui <- fluidPage(
   tags$style("
              body {
    -moz-transform: scale(0.85, 0.85); /* Moz-browsers */
    zoom: 0.85; /* Other non-webkit browsers */
    zoom: 85%; /* Webkit browsers */
	}
              "),
	
headerPanel(HTML('<img src="deq_dwq_logo.png" height="70" width="199.5"/>'),
	windowTitle="GSL trophic data explorer"),

  # App title ----
	titlePanel("GSL trophic data explorer",
		tags$head(tags$link(rel = "icon", type = "image/png", href = "dwq_logo_small.png"),
                        tags$title("GSL trophic data explorer"))
	),
	
  # Sidebar layout with input and output definitions ----
	sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: ----
    radioButtons(inputId = "analysisType",
                label = "Analysis:",
				choiceNames=c("GSL concentration","GSL mass","GSL-NARS comparison"),
				choiceValues=c(1,3,2)),
		conditionalPanel(
			condition="input.analysisType==1",
			radioButtons(inputId="plot_type","Plot type:",inline=TRUE,
				choiceNames=c("Time-series","Boxplots","Scatter plots"),
				choiceValues=c(2,1,3)),
			sliderInput(inputId="years","Year range:",min=min(gsl_data$Year),max=max(gsl_data$Year),value=c(2010,2017),step=1,sep=""),
			checkboxGroupInput("surf_deep","Include:",choiceNames=c("Surface (<=1 m)","Deep (>=4 m)"),choiceValues=c("Surface","Deep"),selected=c("Surface","Deep"),inline=TRUE),
			conditionalPanel(
				condition="input.plot_type==1",
				selectInput("aggby","Aggregate by:", choices=c("Site","Year","Depth"),selected="Year"),
				selectInput("param1","Parameter 1:",choices=param_names,selected="Total_phosphorus_mgL"),
				checkboxInput(inputId="logy1","log(y)",value=TRUE),
				selectInput("param2","Parameter 2:",choices=param_names,selected="Total_nitrogen_mgL"),
				checkboxInput(inputId="logy2","log(y)",value=TRUE),
				selectInput("param3","Parameter 3:",choices=param_names,selected="Salinity_pct"),
				checkboxInput(inputId="logy3","log(y)",value=FALSE)
			),
			conditionalPanel(	
				condition="input.plot_type==2",
				selectInput("param12","Parameter 1:",choices=param_names,selected="Total_phosphorus_mgL"),
				checkboxInput(inputId="logy12","log(y)",value=TRUE),
				selectInput("param22","Parameter 2:",choices=param_names,selected="Total_nitrogen_mgL"),
				checkboxInput(inputId="logy22","log(y)",value=TRUE),
				selectInput("param32","Parameter 3:",choices=param_names,selected="Salinity_pct"),
				checkboxInput(inputId="logy32","log(y)",value=FALSE)
			),
			conditionalPanel(	
				condition="input.plot_type==3",
				fluidRow(
					column(4,selectInput("sp_y1","Y variable 1:",choices=param_names,selected="Total_phosphorus_mgL")),
					column(2,checkboxInput("sp_logy1","log(y1)",value=TRUE)),
					column(4,selectInput("sp_x1","X variable 1:",choices=param_names,selected="Salinity_pct")),
					column(2,checkboxInput("sp_logx1","log(x1)",value=FALSE))
				),

				fluidRow(
					column(4,selectInput("sp_y2","Y variable 2:",choices=param_names,selected="Filtered_phosphorus_mgL")),
					column(2,checkboxInput("sp_logy2","log(y2)",value=TRUE)),
					column(4,selectInput("sp_x2","X variable 2:",choices=param_names,selected="Salinity_pct")),
					column(2,checkboxInput("sp_logx2","log(x2)",value=FALSE))
				),
				fluidRow(
					column(4,selectInput("sp_y3","Y variable 3:",choices=param_names,selected="Total_nitrogen_mgL")),
					column(2,checkboxInput("sp_logy3","log(y3)",value=TRUE)),
					column(4,selectInput("sp_x3","X variable 3:",choices=param_names,selected="Salinity_pct")),
					column(2,checkboxInput("sp_logx3","log(x3)",value=FALSE))
				),
				fluidRow(
					column(4,selectInput("sp_y4","Y variable 4:",choices=param_names,selected="Filtered_nitrogen_mgL")),
					column(2,checkboxInput("sp_logy4","log(y4)",value=TRUE)),
					column(4,selectInput("sp_x4","X variable 4:",choices=param_names,selected="Salinity_pct")),
					column(2,checkboxInput("sp_logx4","log(x4)",value=FALSE))
				),
				checkboxInput(inputId="sp_show_reg","Show regressions",value=TRUE)
			)
		),
		conditionalPanel(
			condition="input.analysisType==2",
			helpText("This portion of the tool allows comparison between data collected from GSL and two EPA national monitoring programs: the National Lakes Assessment (NLA, 2012) and the National Coastal Condition Assessment (NCCA, 2010)."),
			radioButtons("comp_plot_type","Plot type:",choiceNames=c("Scatter plot","Boxplot"),choiceValues=c(1,2)),
			checkboxGroupInput("comp_include","Include:",choices=c("GSL surface","GSL deep","NLA","NCCA"),selected=c("GSL surface","NLA")),
			selectInput("compSPy","Y variable:",choices=compSP_choices,selected="Chlorophyll_a_ugL"),
			checkboxInput(inputId="complogy","log(y)",value=TRUE),
			conditionalPanel(
				condition="input.comp_plot_type==1",
				selectInput("compSPx","X variable:",choices=compSP_choices,selected="Total_nitrogen_mgL"),
				checkboxInput(inputId="complogx","log(x)",value=TRUE)
			)	
		),
		
	
		conditionalPanel(
			condition="input.analysisType==3",
			helpText("This portion of the tool estimates DBL, UBL, and total nutrient masses in Gilbert Bay following Naftz (2017). DBL and UBL volumes (and therefore constituent masses) outside 2011-2014 are extrapolated from models presented in Naftz (2017). Bear River and Farmington Bays and inflows excluded from these analyses"),
			sliderInput(inputId="mass_years","Year range:",min=min(mass_data$Year),max=max(mass_data$Year),value=c(2011,max(mass_data$Year)),step=1,sep=""),
			selectInput("mass1","Parameter 1:",choices=mass_names,selected="TP"),
			selectInput("mass2","Parameter 2:",choices=mass_names,selected="TN"),
			selectInput("mass3","Parameter 3:",choices=mass_names,selected="salt")
		),
		

		helpText("For help with this tool, or to report a bug, please contact Jake Vander Laan, UDWQ, jvander@utah.gov, (801) 536-4350")
		
		
	),

    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Profile Plots ----
		plotOutput(outputId = "plot1",width="800px",height="800px")
	  )
  )
)



# Define server logic

server <- function(input, output) {

output$plot1=renderPlot({
	if(input$analysisType==1){
		gsl_data1=gsl_data[gsl_data$Depth==input$surf_deep[1]|gsl_data$Depth==input$surf_deep[2],]
		if(length(input$surf_deep)==1){gsl_data1=gsl_data[gsl_data$Depth==input$surf_deep[1],]}else{gsl_data1=gsl_data}
		gsl_data1=gsl_data1[gsl_data1$Year>=input$years[1]&gsl_data1$Year<=input$years[2],]
		gsl_data2=gsl_data[gsl_data$Year>=input$years[1]&gsl_data$Year<=input$years[2],]
		if(input$plot_type==1){
			par(mfrow=c(3,1),mar=c(4.1,6.1,2.1,4.1))
			if(input$logy1==TRUE){logy1="y"}else{logy1=""}
			if(input$logy2==TRUE){logy2="y"}else{logy2=""}
			if(input$logy3==TRUE){logy3="y"}else{logy3=""}
			if(length(input$surf_deep)==0){bpcol=NA}else{bpcol="black"}
			boxplot(gsl_data1[,input$param1]~factor(gsl_data1[,input$aggby]),log=logy1,ylab=names(param_names[param_names==input$param1]),cex.lab=1.5,cex.axis=1.5,border=bpcol)
			if(input$param1=="TNTP_mol"){abline(h=16,lty=2,lwd=2,col="purple")}
			boxplot(gsl_data1[,input$param2]~factor(gsl_data1[,input$aggby]),log=logy2,ylab=names(param_names[param_names==input$param2]),cex.lab=1.5,cex.axis=1.5,border=bpcol)
			if(input$param2=="TNTP_mol"){abline(h=16,lty=2,lwd=2,col="purple")}
			boxplot(gsl_data1[,input$param3]~factor(gsl_data1[,input$aggby]),log=logy3,ylab=names(param_names[param_names==input$param3]),cex.lab=1.5,cex.axis=1.5,border=bpcol)}
			if(input$param3=="TNTP_mol"){abline(h=16,lty=2,lwd=2,col="purple")}
		if(input$plot_type==2){
			par(mfrow=c(3,1),mar=c(2.1,6.1,2.1,4.1))
			if(input$logy12==TRUE){logy12="y"}else{logy12=""}
			if(input$logy22==TRUE){logy22="y"}else{logy22=""}
			if(input$logy32==TRUE){logy32="y"}else{logy32=""}
			plot(gsl_data2[,input$param12]~gsl_data2$Date,ylab=names(param_names[param_names==input$param12]),xlab="",log=logy12,pch=NA,cex.lab=1.5,cex.axis=1.5)
				if(("Deep"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Deep",input$param12]~gsl_data2[gsl_data2$Depth=="Deep","Date"],pch=17,col="orange",cex=2)}
				if(("Surface"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Surface",input$param12]~gsl_data2[gsl_data2$Depth=="Surface","Date"],pch=21,col="blue",cex=2)}
				if(input$param12=="TNTP_mol"){abline(h=16,lty=2,lwd=2,col="purple")}
				legend("topleft",bty='n',pch=c(21,17),col=c("blue","orange"),legend=c("Surface","Deep"),cex=2)
			plot(gsl_data2[,input$param22]~gsl_data2$Date,pch=NA,ylab=names(param_names[param_names==input$param22]),xlab="",log=logy22,cex.lab=1.5,cex.axis=1.5)
				if(("Deep"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Deep",input$param22]~gsl_data2[gsl_data2$Depth=="Deep","Date"],pch=17,col="orange",cex=2)}
				if(("Surface"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Surface",input$param22]~gsl_data2[gsl_data2$Depth=="Surface","Date"],pch=21,col="blue",cex=2)}
				if(input$param22=="TNTP_mol"){abline(h=16,lty=2,lwd=2,col="purple")}
				legend("topleft",bty='n',pch=c(21,17),col=c("blue","orange"),legend=c("Surface","Deep"),cex=2)
			plot(gsl_data2[,input$param32]~gsl_data2$Date,pch=NA,ylab=names(param_names[param_names==input$param32]),xlab="",log=logy32,cex.lab=1.5,cex.axis=1.5)
				if(("Deep"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Deep",input$param32]~gsl_data2[gsl_data2$Depth=="Deep","Date"],pch=17,col="orange",cex=2)}
				if(("Surface"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Surface",input$param32]~gsl_data2[gsl_data2$Depth=="Surface","Date"],pch=21,col="blue",cex=2)}
				if(input$param32=="TNTP_mol"){abline(h=16,lty=2,lwd=2,col="purple")}
				legend("topleft",bty='n',pch=c(21,17),col=c("blue","orange"),legend=c("Surface","Deep"),cex=2)
		}
		if(input$plot_type==3){
			par(mfrow=c(2,2),mar=c(6.1,6.1,2.1,2.1))
			if(input$sp_logx1==TRUE&input$sp_logy1==TRUE){log1='xy'}
			if(input$sp_logx1==TRUE&input$sp_logy1==FALSE){log1='x'}
			if(input$sp_logx1==FALSE&input$sp_logy1==TRUE){log1='y'}
			if(input$sp_logx1==FALSE&input$sp_logy1==FALSE){log1=''}
			if(input$sp_logx2==TRUE&input$sp_logy2==TRUE){log2='xy'}
			if(input$sp_logx2==TRUE&input$sp_logy2==FALSE){log2='x'}
			if(input$sp_logx2==FALSE&input$sp_logy2==TRUE){log2='y'}
			if(input$sp_logx2==FALSE&input$sp_logy2==FALSE){log2=''}
			if(input$sp_logx3==TRUE&input$sp_logy3==TRUE){log3='xy'}
			if(input$sp_logx3==TRUE&input$sp_logy3==FALSE){log3='x'}
			if(input$sp_logx3==FALSE&input$sp_logy3==TRUE){log3='y'}
			if(input$sp_logx3==FALSE&input$sp_logy3==FALSE){log3=''}
			if(input$sp_logx4==TRUE&input$sp_logy4==TRUE){log4='xy'}
			if(input$sp_logx4==TRUE&input$sp_logy4==FALSE){log4='x'}
			if(input$sp_logx4==FALSE&input$sp_logy4==TRUE){log4='y'}
			if(input$sp_logx4==FALSE&input$sp_logy4==FALSE){log4=''}
			lmdata=data.frame(matrix(nrow=0,ncol=dim(gsl_data2)[2]))
			colnames(lmdata)=colnames(gsl_data2)
			if(("Deep"%in%input$surf_deep)==TRUE){lmdata=rbind(lmdata,gsl_data2[gsl_data2$Depth=="Deep",])}
			if(("Surface"%in%input$surf_deep)==TRUE){lmdata=rbind(lmdata,gsl_data2[gsl_data2$Depth=="Surface",])}
			plot(gsl_data2[,input$sp_y1]~gsl_data2[,input$sp_x1],pch=NA,log=log1,ylab=names(param_names[param_names==input$sp_y1]),xlab=names(param_names[param_names==input$sp_x1]),cex.lab=1.5,cex.axis=1.5)
				if(("Deep"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Deep",input$sp_y1]~gsl_data2[gsl_data2$Depth=="Deep",input$sp_x1],pch=17,col="orange",cex=2)}
				if(("Surface"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Surface",input$sp_y1]~gsl_data2[gsl_data2$Depth=="Surface",input$sp_x1],pch=21,col="blue",cex=2)}
				if(input$sp_show_reg==TRUE){
					if(log1==""){lm1=lm(lmdata[,input$sp_y1]~lmdata[,input$sp_x1])}
					if(log1=="x"){lm1=lm(lmdata[,input$sp_y1]~log10(lmdata[,input$sp_x1]))}
					if(log1=="xy"){lm1=lm(log10(lmdata[,input$sp_y1])~log10(lmdata[,input$sp_x1]))}
					if(log1=="y"){lm1=lm(log10(lmdata[,input$sp_y1])~lmdata[,input$sp_x1])}
					abline(lm1$coefficients[1],lm1$coefficients[2],lty=2,lwd=2)
					p1=round(summary(lm1)$coefficients[8],digits=2)
					r21=round(summary(lm1)$adj.r.squared,digits=2)
					if(lm1$coefficients[2]>0){lgdloc="bottomright"}else{lgdloc="topright"}
					legend(lgdloc,legend=c(paste0("r2 = ",r21),paste0("p = ",p1)),bty='n',cex=2)}
				legend("topleft",bty='n',pch=c(21,17),col=c("blue","orange"),legend=c("Surface","Deep"),cex=2)
			plot(gsl_data2[,input$sp_y2]~gsl_data2[,input$sp_x2],pch=NA,log=log2,ylab=names(param_names[param_names==input$sp_y2]),xlab=names(param_names[param_names==input$sp_x2]),cex.lab=1.5,cex.axis=1.5)
				if(("Deep"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Deep",input$sp_y2]~gsl_data2[gsl_data2$Depth=="Deep",input$sp_x2],pch=17,col="orange",cex=2)}
				if(("Surface"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Surface",input$sp_y2]~gsl_data2[gsl_data2$Depth=="Surface",input$sp_x2],pch=21,col="blue",cex=2)}
				legend("topleft",bty='n',pch=c(21,17),col=c("blue","orange"),legend=c("Surface","Deep"),cex=2)
				if(input$sp_show_reg==TRUE){
					if(log2==""){lm2=lm(lmdata[,input$sp_y2]~lmdata[,input$sp_x2])}
					if(log2=="x"){lm2=lm(lmdata[,input$sp_y2]~log10(lmdata[,input$sp_x2]))}
					if(log2=="xy"){lm2=lm(log10(lmdata[,input$sp_y2])~log10(lmdata[,input$sp_x2]))}
					if(log2=="y"){lm2=lm(log10(lmdata[,input$sp_y2])~lmdata[,input$sp_x2])}
					abline(lm2$coefficients[1],lm2$coefficients[2],lty=2,lwd=2)
					p2=round(summary(lm2)$coefficients[8],digits=2)
					r22=round(summary(lm2)$adj.r.squared,digits=2)
					if(lm2$coefficients[2]>0){lgdloc="bottomright"}else{lgdloc="topright"}
					legend(lgdloc,legend=c(paste0("r2 = ",r22),paste0("p = ",p2)),bty='n',cex=2)}
			plot(gsl_data2[,input$sp_y3]~gsl_data2[,input$sp_x3],pch=NA,log=log3,ylab=names(param_names[param_names==input$sp_y3]),xlab=names(param_names[param_names==input$sp_x3]),cex.lab=1.5,cex.axis=1.5)
				if(("Deep"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Deep",input$sp_y3]~gsl_data2[gsl_data2$Depth=="Deep",input$sp_x3],pch=17,col="orange",cex=2)}
				if(("Surface"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Surface",input$sp_y3]~gsl_data2[gsl_data2$Depth=="Surface",input$sp_x3],pch=21,col="blue",cex=2)}
				if(input$sp_show_reg==TRUE){	
					if(log3==""){lm3=lm(lmdata[,input$sp_y3]~lmdata[,input$sp_x3])}
					if(log3=="x"){lm3=lm(lmdata[,input$sp_y3]~log10(lmdata[,input$sp_x3]))}
					if(log3=="xy"){lm3=lm(log10(lmdata[,input$sp_y3])~log10(lmdata[,input$sp_x3]))}
					if(log3=="y"){lm3=lm(log10(lmdata[,input$sp_y3])~lmdata[,input$sp_x3])}
					abline(lm3$coefficients[1],lm3$coefficients[2],lty=2,lwd=2)
					p3=round(summary(lm3)$coefficients[8],digits=2)
					r23=round(summary(lm3)$adj.r.squared,digits=2)
					if(lm3$coefficients[2]>0){lgdloc="bottomright"}else{lgdloc="topright"}
					legend(lgdloc,legend=c(paste0("r2 = ",r23),paste0("p = ",p3)),bty='n',cex=2)}
				legend("topleft",bty='n',pch=c(21,17),col=c("blue","orange"),legend=c("Surface","Deep"),cex=2)
			plot(gsl_data2[,input$sp_y4]~gsl_data2[,input$sp_x4],pch=NA,log=log4,ylab=names(param_names[param_names==input$sp_y4]),xlab=names(param_names[param_names==input$sp_x4]),cex.lab=1.5,cex.axis=1.5)
				if(("Deep"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Deep",input$sp_y4]~gsl_data2[gsl_data2$Depth=="Deep",input$sp_x4],pch=17,col="orange",cex=2)}
				if(("Surface"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Surface",input$sp_y4]~gsl_data2[gsl_data2$Depth=="Surface",input$sp_x4],pch=21,col="blue",cex=2)}
				if(input$sp_show_reg==TRUE){
					if(log4==""){lm4=lm(lmdata[,input$sp_y4]~lmdata[,input$sp_x4])}
					if(log4=="x"){lm4=lm(lmdata[,input$sp_y4]~log10(lmdata[,input$sp_x4]))}
					if(log4=="xy"){lm4=lm(log10(lmdata[,input$sp_y4])~log10(lmdata[,input$sp_x4]))}
					if(log4=="y"){lm4=lm(log10(lmdata[,input$sp_y4])~lmdata[,input$sp_x4])}
					abline(lm4$coefficients[1],lm4$coefficients[2],lty=2,lwd=2)
					p4=round(summary(lm4)$coefficients[8],digits=2)
					r24=round(summary(lm4)$adj.r.squared,digits=2)
					if(lm4$coefficients[2]>0){lgdloc="bottomright"}else{lgdloc="topright"}
					legend(lgdloc,legend=c(paste0("r2 = ",r24),paste0("p = ",p4)),bty='n',cex=2)}
				legend("topleft",bty='n',pch=c(21,17),col=c("blue","orange"),legend=c("Surface","Deep"),cex=2)
		}
	}
	if(input$analysisType==2){
		mar=c(4.1,6.1,2.1,4.1)
		if(input$comp_plot_type==1){
			if(input$complogy==TRUE){
				ymin=min(c(NLA_data[,input$compSPy],NCCA_data[NCCA_data$PARAMETER_2==input$compSPy,"RESULT"],gsl_data[,input$compSPy]),na.rm=T)}else{ymin=0}
			ymax=max(c(NLA_data[,input$compSPy],NCCA_data[NCCA_data$PARAMETER_2==input$compSPy,"RESULT"],gsl_data[,input$compSPy]),na.rm=T)
			if(input$complogx==TRUE){
				xmin=min(c(NLA_data[,input$compSPx],NCCA_data[NCCA_data$PARAMETER_2==input$compSPx,"RESULT"],gsl_data[,input$compSPx]),na.rm=T)}else{xmin=0}
			xmax=max(c(NLA_data[,input$compSPx],NCCA_data[NCCA_data$PARAMETER_2==input$compSPx,"RESULT"],gsl_data[,input$compSPx]),na.rm=T)
			if(input$complogx==TRUE&input$complogy==TRUE){complog='xy'}
			if(input$complogx==TRUE&input$complogy==FALSE){complog='x'}
			if(input$complogx==FALSE&input$complogy==TRUE){complog='y'}
			if(input$complogx==FALSE&input$complogy==FALSE){complog=''}
			plot(NLA_data[,input$compSPy]~NLA_data[,input$compSPx],cex.lab=1.5,cex.axis=1.5,ylim=c(ymin,ymax),xlim=c(xmin,xmax),log=complog,pch=NA,ylab=names(compSP_choices[compSP_choices==input$compSPy]),xlab=names(compSP_choices[compSP_choices==input$compSPx]))
				legend_text=vector()
				legend_pch=vector()
				legend_col=vector()
				if(complog=="xy"){NLA_lm=lm(log10(NLA_data[,input$compSPy])~log10(NLA_data[,input$compSPx]))}
				if(complog=="x"){NLA_lm=lm(NLA_data[,input$compSPy]~log10(NLA_data[,input$compSPx]))}
				if(complog=="y"){NLA_lm=lm(log10(NLA_data[,input$compSPy])~NLA_data[,input$compSPx])}
				if(complog==""){NLA_lm=lm(NLA_data[,input$compSPy]~NLA_data[,input$compSPx])}
				if(("NLA"%in%input$comp_include)==TRUE){
					points(NLA_data[,input$compSPy]~NLA_data[,input$compSPx],col="gray40",cex=2)
					legend_text=append(legend_text,"NLA")
					legend_pch=append(legend_pch,21)
					legend_col=append(legend_col,"gray40")
					abline(NLA_lm$coefficients[1],NLA_lm$coefficients[2],lty=2,lwd=2,col="gray40")
				}
				if(("NCCA"%in%input$comp_include)==TRUE){
					points(NCCA_data_mat[,input$compSPy]~NCCA_data_mat[,input$compSPx],pch=22,col="springgreen4",cex=2)
					legend_text=append(legend_text,"NCCA")
					legend_pch=append(legend_pch,22)
					legend_col=append(legend_col,"springgreen4")
					if(complog=="xy"){NCCA_lm=lm(log10(NCCA_data_mat[,input$compSPy])~log10(NCCA_data_mat[,input$compSPx]))}
					if(complog=="x"){NCCA_lm=lm(NCCA_data_mat[,input$compSPy]~log10(NCCA_data_mat[,input$compSPx]))}
					if(complog=="y"){NCCA_lm=lm(log10(NCCA_data_mat[,input$compSPy])~NCCA_data_mat[,input$compSPx])}
					if(complog==""){NCCA_lm=lm(NCCA_data_mat[,input$compSPy]~NCCA_data_mat[,input$compSPx])}
					abline(NCCA_lm$coefficients[1],NCCA_lm$coefficients[2],lty=2,lwd=2,col="springgreen4")
				}
				if(("GSL deep"%in%input$comp_include)==TRUE){
					points(gsl_data[gsl_data$Depth=="Deep",input$compSPy]~gsl_data[gsl_data$Depth=="Deep",input$compSPx],pch=17,col="orange",cex=2)
					legend_text=append(legend_text,"GSL deep")
					legend_pch=append(legend_pch,17)
					legend_col=append(legend_col,"orange")
				}
				if(("GSL surface"%in%input$comp_include)==TRUE){
					points(gsl_data[gsl_data$Depth=="Surface",input$compSPy]~gsl_data[gsl_data$Depth=="Surface",input$compSPx],pch=19,col="blue",cex=2)
					legend_text=append(legend_text,"GSL surface")
					legend_pch=append(legend_pch,19)
					legend_col=append(legend_col,"blue")
				}
			if(input$compSPy=="TNTP_mol"){abline(h=16,lty=2,lwd=2,col="purple")}
			if(NLA_lm$coefficients[2]>0){
				legend("topleft",pch=legend_pch,col=legend_col,legend=legend_text,bty='n',cex=2)
				}else{
				legend("topright",pch=legend_pch,col=legend_col,legend=legend_text,bty='n',cex=2)
				}
			
		}
		if(input$comp_plot_type==2){
			if(input$complogy==TRUE){compbplogy='y'}else{compbplogy=""}
			bp_result=vector()
			bp_group=vector()
			if(("NLA"%in%input$comp_include)==TRUE){
				result=NLA_data[,input$compSPy]
				group=rep("NLA",length(result))
				bp_result=append(bp_result,result)
				bp_group=append(bp_group,group)
			}
			if(("NCCA"%in%input$comp_include)==TRUE){
				result=NCCA_data_mat[,input$compSPy]
				group=rep("NCCA",length(result))
				bp_result=append(bp_result,result)
				bp_group=append(bp_group,group)
			}
			if(("GSL deep"%in%input$comp_include)==TRUE){
				result=gsl_data[gsl_data$Depth=="Deep",input$compSPy]
				group=rep("GSL deep",length(result))
				bp_result=append(bp_result,result)
				bp_group=append(bp_group,group)
			}
			if(("GSL surface"%in%input$comp_include)==TRUE){
				result=gsl_data[gsl_data$Depth=="Surface",input$compSPy]
				group=rep("GSL surface",length(result))
				bp_result=append(bp_result,result)
				bp_group=append(bp_group,group)
			}
			boxplot(bp_result~bp_group,log=compbplogy,ylab=names(compSP_choices[compSP_choices==input$compSPy]),cex.axis=1.5,cex.lab=1.5)
			if(input$compSPy=="TNTP_mol"){abline(h=16,lty=2,lwd=2,col="purple")}
		}
	}
	if(input$analysisType==3){
		options(scipen=1)
		mass_data2=mass_data[mass_data$Year>=input$mass_years[1]&mass_data$Year<=input$mass_years[2],]
		par(mfrow=c(3,1),mar=c(4.1,6.1,2.1,4.1))
		label1=paste0(names(mass_names[mass_names==input$mass1])," (kg)")
		massPlot(param=input$mass1,mass_data=mass_data2,label=label1)
		label2=paste0(names(mass_names[mass_names==input$mass2])," (kg)")
		massPlot(param=input$mass2,mass_data=mass_data2,label=label2)
		label3=paste0(names(mass_names[mass_names==input$mass3])," (kg)")
		massPlot(param=input$mass3,mass_data=mass_data2,label=label3)
	}
})
}
shinyApp(ui = ui, server = server)
