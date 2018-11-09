#setwd("F:\\Shiny\\GSL_data_explorer")

gsl_data=read.csv(file="data/GSL_WQ.csv")
NLA_data=read.csv(file="data/NLA2012_TrophicData.csv")
NCCA_data=read.csv(file="data/assessed_ncca2010_waterchem.csv")
metalsData=read.csv(file="data/metalsData.csv")
combelv=read.csv(file="data/combined_elev.csv")
combelv$Date=as.Date(combelv$Date,format='%Y-%m-%d')

library(shiny)
library(lubridate)
library(reshape2)
library(plyr)
library(waterData)

colnames(gsl_data)[colnames(gsl_data)=="site_name"]="Site"
options(scipen=999)
gsl_data$Date=as.Date(gsl_data$sample_dt,format='%Y-%m-%d')
gsl_data$Year=year(gsl_data$Date)
gsl_data$Month=month(gsl_data$Date)
gsl_data$Depth=rep(NA,length=dim(gsl_data)[1])
gsl_data$Depth[gsl_data$Sampling_depth_m<=1]="Surface"
gsl_data$Depth[gsl_data$Sampling_depth_m>=4]="Deep"
gsl_data$TNTP_mol=gsl_data$Total_nitrogen_mgL/gsl_data$Total_phosphorus_mgL*2.211323
coords=read.csv(file="data/GSL_BaselineSiteLocations.csv")
gsl_data=merge(gsl_data,coords,by="Site")


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
"Salinity_gL",
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
"Salinity (g/L)",
"TN:TP (mol)")

mass_names=c("TP","FP","TN","FN","salt")
names(mass_names)=c("Total phosphorus","Filtered phosphorus","Total nitrogen","Filtered nitrogen","Salt")

NCCA_data_mat=dcast(NCCA_data,UID~PARAMETER_2,value.var="RESULT")
NCCA_data_mat$Secchi_m=rep(NA,dim(NCCA_data_mat)[1])
NCCA_data_mat$TNTP_mol=NCCA_data_mat$Total_nitrogen_mgL/NCCA_data_mat$Total_phosphorus_mgL*2.211323

source("functions/nut_mass_preprocess.R")

#Manipulate metals data
	#screen out:
		#Q/A data
		#non-target parameters
		#non-sites
	#Convert dates
		#add year
	#Convert ww brine shrimp to dw
	#Convert all units to ug/L or mg/kg (divide ng/L and ng/g dry by 1000)
	#?Calc % meHg (flip to matrix, calculate, flip back to flat, append to data frame)
	
metalsData=metalsData[!is.na(metalsData$SiteID)&metalsData$Matrix!="Water-T"&!is.na(metalsData$Parameter),]
metalsData$Date=as.Date(metalsData$Sample_Date,format='%m/%d/%Y')
metalsData$Year=year(metalsData$Date)
metalsData$Month=month(metalsData$Date)
dw_result=ifelse(metalsData$Units=="mg/kg wet"|metalsData$Units=="ng/g wet",metalsData$Result/0.15,metalsData$Result)
metalsData$Result=dw_result
metalsData$Units[metalsData$Units=="mg/kg wet"]="mg/kg dry"
metalsData$Units[metalsData$Units=="ng/g wet"]="ng/g dry"
result=ifelse(metalsData$Units=="ng/L"|metalsData$Units=="ng/g dry",metalsData$Result/1000,metalsData$Result)
metalsData$Result=result
metalsData$Units[metalsData$Units=="ng/L"]="ug/L"
metalsData$Units[metalsData$Units=="ng/g dry"]="mg/kg dry"
tmeHg_mat=dcast(metalsData[metalsData$Parameter=="Hg"|metalsData$Parameter=="MeHg",],Lab_Name+Sample_Tag+Project_ID+SiteID+Date+WaterColLoc+Matrix+Lab_ID+Year~Parameter,value.var="Result")
pct_meHg=tmeHg_mat$MeHg/tmeHg_mat$Hg*100
tmeHg_mat=cbind(tmeHg_mat,pct_meHg)
pct_meHg_flat=melt(tmeHg_mat,id.vars=c("Lab_Name","Sample_Tag","Project_ID","SiteID","Date","WaterColLoc","Matrix","Lab_ID","Year"))
names(pct_meHg_flat)[names(pct_meHg_flat)=='variable']='Parameter'
names(pct_meHg_flat)[names(pct_meHg_flat)=='value']='Result'
pct_meHg_flat=pct_meHg_flat[pct_meHg_flat$Parameter=="pct_meHg",]
pct_meHg_flat$Parameter="% MeHg"
pct_meHg_flat$Units="pct"
#head(pct_meHg_flat)
metalsData=rbind.fill(metalsData,pct_meHg_flat)
metalsData=metalsData[metalsData$Parameter!="% MeHg",]

source("functions/makeMetalsSummaryTable.R")

load("polyrast/gsl_poly.rdata")

# Define UI for app ----
ui <- fluidPage(
   tags$style("
              body {
    -moz-transform: scale(0.85, 0.85); /* Moz-browsers */
    zoom: 0.85; /* Other non-webkit browsers */
    zoom: 85%; /* Webkit browsers */
	}
              "),
	
#headerPanel(HTML('<img src="deq_dwq_logo.png" height="70" width="199.5"/>'),
headerPanel( title=div(img(src="deq_dwq_logo1.png", height = 125, width = 125*2.89/1.47)),
	windowTitle="GSL data explorer"),

#navbarPage("My application",
	tabsetPanel(id="tabs",
		tabPanel("Lake Levels",value=4),
		tabPanel("Water Quality Data",value=1),
		tabPanel("Metals Data",value=2),
		tabPanel("Water Quality Map",value=3)
		),

#App title ----
titlePanel("",
		tags$head(tags$link(rel = "icon", type = "image/png", href = "dwq_logo_small.png"),
                        tags$title("GSL data explorer"))
	),

conditionalPanel(
	condition="input.tabs==1",
	
  # Sidebar layout with input and output definitions ----
	sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: ----
    radioButtons(inputId = "analysisType",
                label = "Analysis:",
				choiceNames=c("Concentration","Mass","GSL-NARS comparison"),
				choiceValues=c(1,3,2),inline=TRUE),
		conditionalPanel(
			condition="input.analysisType==1",
			radioButtons(inputId="plot_type","Plot type:",inline=TRUE,
				choiceNames=c("Time-series","Boxplots","Scatter plots"),
				choiceValues=c(2,1,3)),
			sliderInput(inputId="years","Year range:",min=min(gsl_data$Year),max=max(gsl_data$Year),value=c(2010,max(gsl_data$Year)),step=1,sep=""),
			checkboxGroupInput("surf_deep","Include:",choiceNames=c("Surface (<=1 m)","Deep (>=4 m)"),choiceValues=c("Surface","Deep"),selected=c("Surface","Deep"),inline=TRUE),
			conditionalPanel(
				condition="input.plot_type==1",
				selectInput("aggby","Aggregate by:", choices=c("Site","Year","Month","Depth"),selected="Year"),
				selectInput("param1","Parameter 1:",choices=param_names,selected="Total_phosphorus_mgL"),
				checkboxInput(inputId="logy1","log(y)",value=TRUE),
				selectInput("param2","Parameter 2:",choices=param_names,selected="Total_nitrogen_mgL"),
				checkboxInput(inputId="logy2","log(y)",value=TRUE),
				selectInput("param3","Parameter 3:",choices=param_names,selected="Salinity_gL"),
				checkboxInput(inputId="logy3","log(y)",value=FALSE)
			),
			conditionalPanel(	
				condition="input.plot_type==2",
				selectInput("param12","Parameter 1:",choices=param_names,selected="Total_phosphorus_mgL"),
				checkboxInput(inputId="logy12","log(y)",value=TRUE),
				selectInput("param22","Parameter 2:",choices=param_names,selected="Total_nitrogen_mgL"),
				checkboxInput(inputId="logy22","log(y)",value=TRUE),
				selectInput("param32","Parameter 3:",choices=param_names,selected="Salinity_gL"),
				checkboxInput(inputId="logy32","log(y)",value=FALSE)
			),
			conditionalPanel(	
				condition="input.plot_type==3",
				fluidRow(
					column(4,selectInput("sp_y1","Y variable 1:",choices=param_names,selected="Total_phosphorus_mgL")),
					column(2,checkboxInput("sp_logy1","log(y1)",value=TRUE)),
					column(4,selectInput("sp_x1","X variable 1:",choices=param_names,selected="Salinity_gL")),
					column(2,checkboxInput("sp_logx1","log(x1)",value=FALSE))
				),

				fluidRow(
					column(4,selectInput("sp_y2","Y variable 2:",choices=param_names,selected="Filtered_phosphorus_mgL")),
					column(2,checkboxInput("sp_logy2","log(y2)",value=TRUE)),
					column(4,selectInput("sp_x2","X variable 2:",choices=param_names,selected="Salinity_gL")),
					column(2,checkboxInput("sp_logx2","log(x2)",value=FALSE))
				),
				fluidRow(
					column(4,selectInput("sp_y3","Y variable 3:",choices=param_names,selected="Total_nitrogen_mgL")),
					column(2,checkboxInput("sp_logy3","log(y3)",value=TRUE)),
					column(4,selectInput("sp_x3","X variable 3:",choices=param_names,selected="Salinity_gL")),
					column(2,checkboxInput("sp_logx3","log(x3)",value=FALSE))
				),
				fluidRow(
					column(4,selectInput("sp_y4","Y variable 4:",choices=param_names,selected="Filtered_nitrogen_mgL")),
					column(2,checkboxInput("sp_logy4","log(y4)",value=TRUE)),
					column(4,selectInput("sp_x4","X variable 4:",choices=param_names,selected="Salinity_gL")),
					column(2,checkboxInput("sp_logx4","log(x4)",value=FALSE))
				),
				checkboxInput(inputId="sp_show_reg","Show regressions",value=TRUE)
			)
		),
		conditionalPanel(
			condition="input.analysisType==2",
			helpText("This tool allows comparison of nutrient and chlorophyll a relationships between GSL and two EPA National Aquatic Resource Surveys (NARS): the National Lakes Assessment (NLA, 2012) and the National Coastal Condition Assessment (NCCA, 2010)."),
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
			helpText("This tool estimates DBL, UBL, and total nutrient masses in Gilbert Bay following Naftz (2017) with small modifications. DBL and UBL volumes (and therefore constituent masses) outside 2011-2014 are extrapolated from models presented in Naftz (2017). Bear River and Farmington Bays and inflows currently excluded from these analyses."),
			sliderInput(inputId="mass_years","Year range:",min=min(mass_data$Year),max=max(mass_data$Year),value=c(min(mass_data$Year),max(mass_data$Year)),step=1,sep=""),
			selectInput("mass1","Parameter 1:",choices=mass_names,selected="TP"),
			selectInput("mass2","Parameter 2:",choices=mass_names,selected="TN"),
			selectInput("mass3","Parameter 3:",choices=mass_names,selected="salt")
		),
		

		helpText("For help with this tool, or to report a bug, please contact Jake Vander Laan, UDWQ, jvander@utah.gov, (801) 536-4350.")
		
		
	),

    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Profile Plots ----
		plotOutput(outputId = "plot1",width="800px",height="800px")
	  )
	)
),


conditionalPanel(
	condition="input.tabs==2",
	
  # Sidebar layout with input and output definitions ----
	sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: ----
    radioButtons(inputId = "metalsPlotType",
                label = "Plot type:",
				choiceNames=c("Time-series","Boxplots"),
				choiceValues=c(1,2),inline=TRUE),
	sliderInput(inputId="years_metals","Year range:",min=min(metalsData$Year),max=max(metalsData$Year),value=c(min(metalsData$Year),max(metalsData$Year)),step=1,sep=""),
	checkboxGroupInput("surf_deep_metals","Include:",choiceNames=c("Surface (0.2 m from surface)","Deep (0.5 m from bottom)"),choiceValues=c("Surface","Deep"),selected=c("Surface","Deep"),inline=TRUE),
	
	#Conditional panel for metals time series
	conditionalPanel(condition="input.metalsPlotType==1",
		radioButtons(inputId = "metalsTS_mat1",label="Matrix 1:",choiceNames=c("Water","Brine shrimp","Bird eggs"),choiceValues=c(1,2,3),inline=TRUE),
		uiOutput("reactiveMetalsChoices1"),
		checkboxInput("log_metalsTS1","log(y1)",value=FALSE),
		radioButtons(inputId = "metalsTS_mat2",label="Matrix 2:",choiceNames=c("Water","Brine shrimp","Bird eggs"),choiceValues=c(1,2,3),inline=TRUE),
		uiOutput("reactiveMetalsChoices2"),
		checkboxInput("log_metalsTS2","log(y2)",value=FALSE),
		radioButtons(inputId = "metalsTS_mat3",label="Matrix 3:",choiceNames=c("Water","Brine shrimp","Bird eggs"),choiceValues=c(1,2,3),inline=TRUE),
		uiOutput("reactiveMetalsChoices3"),
		checkboxInput("log_metalsTS3","log(y3)",value=FALSE),
		helpText("All biological oncentrations reported in dry weight.")

	),
	
	
	#Conditional panel (for metals boxplots inputs)
	conditionalPanel(condition="input.metalsPlotType==2",
		radioButtons(inputId = "metalsBP_mat",label="Matrix:",choiceNames=c("Water","Brine shrimp","Bird eggs"),choiceValues=c(1,2,3),inline=TRUE),
		uiOutput("reactiveMetalsChoicesBP"),
		uiOutput("reactiveAggByBP"),
		checkboxInput("log_metalsBP","log(y)",value=FALSE),
		helpText("AMAV = American Avocet. BNST = Black-necked stilt. AIC = Antelope Island Causeway. ANTI = Antelope Island. FARM = Farmington Bay. OGBA = Ogden Bay. SALT = Saltair. All biological oncentrations reported in dry weight.")
		
	),
	helpText("For help with this tool, or to report a bug, please contact Jake Vander Laan, UDWQ, jvander@utah.gov, (801) 536-4350.")

	),
	
	# Main panel for displaying outputs ----
    mainPanel(
      # Output: Profile Plots ----
		conditionalPanel(
			condition="input.metalsPlotType==1",
			plotOutput(outputId = "plot2",width="800px",height="800px")
		),
		conditionalPanel(
			condition="input.metalsPlotType==2",
			plotOutput(outputId = "plot3",width="800px",height="400px"),
			tableOutput("table3")
		)

	  )
	)
	),
	
	
	conditionalPanel(
		condition="input.tabs==3",
	
		#Sidebar layout with input and output definitions ----
		sidebarLayout(
	
		# Sidebar panel for inputs ----
		sidebarPanel(
		helpText("This tool estimates lake-wide values for selected parameters, depths, and time periods via inverse path distance weighted interpolation."),
		helpText("Please be patient - this type of interpolation accounts for land-barriers between sites, but is  data intensive and therefore slow."),
		helpText("Select desired parameters, then click 'Interpolate' to generate the map."),
		helpText("Numbers plotted on the map show available sample sizes at each site used for interpolation."),
	 		
		# Input: ----
		sliderInput(inputId="wq_map_dates","Date range:",min=as.Date("2010-01-01"),max=max(gsl_data$Date),value=c(as.Date("2010-01-01"),max(gsl_data$Date))),
		checkboxGroupInput("wq_map_depths","Include:",choiceNames=c("Surface (0.2 m from surface)","Deep (0.5 m from bottom)"),choiceValues=c("Surface","Deep"),selected=c("Surface"),inline=TRUE),
		selectInput("wq_map_param","Parameter:",choices=param_names,selected="Salinity_gL"),
		checkboxInput("wq_map_uniform_scale","Maintain uniform z-scale (for time-series and depth comparisons)",value=FALSE),
		actionButton("interpolate", "Interpolate", style='height:50px; padding:4px; font-size:150%',width="200px"),
		br(),
		br(),
		helpText("For help with this tool, or to report a bug, please contact Jake Vander Laan, UDWQ, jvander@utah.gov, (801) 536-4350.")
		),
		
		# Main panel for displaying outputs ----
		mainPanel(
			plotOutput(outputId = "wq_map",width="800px",height="800px")
		)
	)
	),
	
	
	
	conditionalPanel(
		condition="input.tabs==4",
		
		#Sidebar layout with input and output definitions ----
		sidebarLayout(
	
			# Sidebar panel for inputs ----
			sidebarPanel(
			helpText("This tool plots Gilbert and Gunnison Bay lake elevations, estimated volumes, and south to north elevation difference (head). Gilbert Bay volume estimated following Naftz (2017). Gunnison Bay volume estimated following Baskin (2006)"),
			br(),

			# Input: ----
			sliderInput(inputId="elev_plot_dates","Date range:",min=min(as.Date(combelv$Date)),max=max(combelv$Date),value=c(min(combelv$Date),max(gsl_data$Date))),

			radioButtons(inputId = "elev_legend",
                label = "Legends:",
				choiceNames=c("On","Off"),
				choiceValues=c(1,0),selected=1,inline=TRUE),

			radioButtons(inputId = "elev_rulers",
                label = "Rulers:",
				choiceNames=c("On","Off"),
				choiceValues=c(1,0),selected=0,inline=TRUE),
			
			conditionalPanel(
				condition="input.elev_rulers==1",
				sliderInput(inputId="elev_ruler1","Elevation ruler:",
					min=floor(min(min(combelv$value[combelv$variable=="SA_LakeElevFt"|combelv$variable=="NA_LakeElevFt"],na.rm=TRUE),min(combelv$value[combelv$variable=="SA_LakeElevFt"|combelv$variable=="NA_LakeElevFt"],na.rm=TRUE))),
					max=ceiling(max(max(combelv$value[combelv$variable=="SA_LakeElevFt"|combelv$variable=="NA_LakeElevFt"],na.rm=TRUE),max(combelv$value[combelv$variable=="SA_LakeElevFt"|combelv$variable=="NA_LakeElevFt"],na.rm=TRUE))),
					value=mean(combelv$value[combelv$variable=="SA_LakeElevFt"|combelv$variable=="NA_LakeElevFt"],na.rm=TRUE)),
				sliderInput(inputId="elev_ruler2","Volume ruler:",
					min=(floor(min(min(combelv$value[combelv$variable=="SA_vol_acreft"|combelv$variable=="NA_vol_acreft"],na.rm=TRUE),min(combelv$value[combelv$variable=="SA_vol_acreft"|combelv$variable=="NA_vol_acreft"],na.rm=TRUE))))/1000,
					max=(ceiling(max(max(combelv$value[combelv$variable=="SA_vol_acreft"|combelv$variable=="NA_vol_acreft"],na.rm=TRUE),max(combelv$value[combelv$variable=="SA_vol_acreft"|combelv$variable=="NA_vol_acreft"],na.rm=TRUE))))/1000,
					value=mean(combelv$value[combelv$variable=="SA_vol_acreft"|combelv$variable=="NA_vol_acreft"],na.rm=TRUE)/1000),
				sliderInput(inputId="elev_ruler3","Head ruler:",
					min=floor(min(min(combelv$value[combelv$variable=="SmNhead"],na.rm=TRUE),min(combelv$value[combelv$variable=="SmNhead"],na.rm=TRUE))),
					max=ceiling(max(max(combelv$value[combelv$variable=="SmNhead"],na.rm=TRUE),max(combelv$value[combelv$variable=="SmNhead"],na.rm=TRUE))),
					value=0)
			),
			
			
			br(),
			helpText("For help with this tool, or to report a bug, please contact Jake Vander Laan, UDWQ, jvander@utah.gov, (801) 536-4350.")
			),
			
			mainPanel(
				plotOutput(outputId = "elev_plot",width="800px",height="800px")
			)
		)
	)
)


# Define server logic
server <- function(input, output){



#user <- unname(Sys.info()["user"])
#if (user == "shiny") {
#
#  # Set library locations
#  .libPaths(c(
#    "C:/Users/jvander/Documents/R/R-3.4.4/library"
#  )
#  )
#
#}

library(sp)
library(RColorBrewer)
library(raster)
library(rgdal)
library(gdalUtils)
library(mapview)
library(Matrix)
library(gdistance)
library(ipdw)




###
###Plot 1 (input$tabs==1)
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
			par(mfrow=c(3,1),mar=c(4.1,6.1,2.1,4.1))
			if(input$logy12==TRUE){logy12="y"}else{logy12=""}
			if(input$logy22==TRUE){logy22="y"}else{logy22=""}
			if(input$logy32==TRUE){logy32="y"}else{logy32=""}
			if(all(is.na(gsl_data2[,input$param12]))){frame()}else{
			plot(gsl_data2[,input$param12]~gsl_data2$Date,ylab=names(param_names[param_names==input$param12]),xlab="",log=logy12,pch=NA,cex.axis=2,cex.lab=2.25)
				if(("Deep"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Deep",input$param12]~gsl_data2[gsl_data2$Depth=="Deep","Date"],pch=17,col="orange",cex=2)}
				if(("Surface"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Surface",input$param12]~gsl_data2[gsl_data2$Depth=="Surface","Date"],pch=21,col="blue",cex=2)}
				if(input$param12=="TNTP_mol"){abline(h=16,lty=2,lwd=2,col="purple")}
				legend("topleft",bty='n',pch=c(21,17),col=c("blue","orange"),legend=c("Surface","Deep"),cex=2)}
			if(all(is.na(gsl_data2[,input$param22]))){frame()}else{
			plot(gsl_data2[,input$param22]~gsl_data2$Date,pch=NA,ylab=names(param_names[param_names==input$param22]),xlab="",log=logy22,cex.axis=2,cex.lab=2.25)
				if(("Deep"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Deep",input$param22]~gsl_data2[gsl_data2$Depth=="Deep","Date"],pch=17,col="orange",cex=2)}
				if(("Surface"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Surface",input$param22]~gsl_data2[gsl_data2$Depth=="Surface","Date"],pch=21,col="blue",cex=2)}
				if(input$param22=="TNTP_mol"){abline(h=16,lty=2,lwd=2,col="purple")}
				legend("topleft",bty='n',pch=c(21,17),col=c("blue","orange"),legend=c("Surface","Deep"),cex=2)}
			if(all(is.na(gsl_data2[,input$param32]))){frame()}else{
			plot(gsl_data2[,input$param32]~gsl_data2$Date,pch=NA,ylab=names(param_names[param_names==input$param32]),xlab="",log=logy32,cex.axis=2,cex.lab=2.25)
				if(("Deep"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Deep",input$param32]~gsl_data2[gsl_data2$Depth=="Deep","Date"],pch=17,col="orange",cex=2)}
				if(("Surface"%in%input$surf_deep)==TRUE){points(gsl_data2[gsl_data2$Depth=="Surface",input$param32]~gsl_data2[gsl_data2$Depth=="Surface","Date"],pch=21,col="blue",cex=2)}
				if(input$param32=="TNTP_mol"){abline(h=16,lty=2,lwd=2,col="purple")}
				legend("topleft",bty='n',pch=c(c(21,17),17),col=c("blue","orange"),legend=c("Surface","Deep"),cex=2)}
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
			if(dim(gsl_data2)[1]>0&dim(lmdata)[1]>0){
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
			if(any(input$comp_include%in%c("NLA","NCCA","GSL deep","GSL surface"))){
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
			if(length(bp_result)>0){
				boxplot(bp_result~bp_group,log=compbplogy,ylab=names(compSP_choices[compSP_choices==input$compSPy]),cex.axis=1.5,cex.lab=1.5)
				if(input$compSPy=="TNTP_mol"){abline(h=16,lty=2,lwd=2,col="purple")}
			}
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

###
###Plot 2 (input$tabs==2)

#Reactive parameter selectInputs
output$reactiveMetalsChoices1 <- renderUI({
	if(input$metalsTS_mat1==1){
	metalsTS_choices1=unique(metalsData[metalsData$Matrix=="Water",]$Parameter)}
	if(input$metalsTS_mat1==2){
	metalsTS_choices1=unique(metalsData[metalsData$Matrix=="Brine shrimp",]$Parameter)}
	if(input$metalsTS_mat1==3){
	metalsTS_choices1=unique(metalsData[metalsData$Matrix=="Bird eggs",]$Parameter)}
	selectInput("metalsTS_param1", "Parameter 1:", metalsTS_choices1)
	})
output$reactiveMetalsChoices2 <- renderUI({	
	if(input$metalsTS_mat2==1){
	metalsTS_choices2=unique(metalsData[metalsData$Matrix=="Water",]$Parameter)}
	if(input$metalsTS_mat2==2){
	metalsTS_choices2=unique(metalsData[metalsData$Matrix=="Brine shrimp",]$Parameter)}
	if(input$metalsTS_mat2==3){
	metalsTS_choices2=unique(metalsData[metalsData$Matrix=="Bird eggs",]$Parameter)}
	selectInput("metalsTS_param2", "Parameter 2:", metalsTS_choices2)
	})
output$reactiveMetalsChoices3 <- renderUI({	
	if(input$metalsTS_mat3==1){
	metalsTS_choices3=unique(metalsData[metalsData$Matrix=="Water",]$Parameter)}
	if(input$metalsTS_mat3==2){
	metalsTS_choices3=unique(metalsData[metalsData$Matrix=="Brine shrimp",]$Parameter)}
	if(input$metalsTS_mat3==3){
	metalsTS_choices3=unique(metalsData[metalsData$Matrix=="Bird eggs",]$Parameter)}
	selectInput("metalsTS_param3", "Parameter 3:", metalsTS_choices3)
	})
output$reactiveMetalsChoicesBP <- renderUI({
	if(input$metalsBP_mat==1){
	metalsBP_choices=unique(metalsData[metalsData$Matrix=="Water",]$Parameter)}
	if(input$metalsBP_mat==2){
	metalsBP_choices=unique(metalsData[metalsData$Matrix=="Brine shrimp",]$Parameter)}
	if(input$metalsBP_mat==3){
	metalsBP_choices=unique(metalsData[metalsData$Matrix=="Bird eggs",]$Parameter)}
	selectInput("metalsBP_param", "Parameter:", metalsBP_choices)
	})
output$reactiveAggByBP <- renderUI({
	if(input$metalsBP_mat==1|input$metalsBP_mat==2){
	metalsBP_aggby_choices=c("SiteID","Year","Month")
	names(metalsBP_aggby_choices)=c("Site","Year","Month")}
	if(input$metalsBP_mat==3){
	metalsBP_aggby_choices=c("SiteID","Year","Month","Spp")
	names(metalsBP_aggby_choices)=c("Site","Year","Month","Species")}
	selectInput("metalsBP_aggby","Aggregate by:",metalsBP_aggby_choices)
	})


output$plot2=renderPlot({

	#Subset by year input
	metalsData_plots=metalsData[metalsData$Year>=input$years_metals[1]&metalsData$Year<=input$years_metals[2],]
	metalsData_plots=metalsData[metalsData$Year>=input$years_metals[1]&metalsData$Year<=input$years_metals[2],]
	
	#Subset by surface/deep input
	surf_deep_metals=input$surf_deep_metals
	surf_deep_metals[surf_deep_metals=="Surface"]="S"
	surf_deep_metals[surf_deep_metals=="Deep"]="D"
	surf_deep_metals=append(surf_deep_metals,rep("BIO",3-length(surf_deep_metals)))
	metalsData_plots=metalsData_plots[metalsData_plots$WaterColLoc==surf_deep_metals[1]|metalsData_plots$WaterColLoc==surf_deep_metals[2]|metalsData_plots$WaterColLoc==surf_deep_metals[3],]

if(input$metalsPlotType==1){
	matrices_num=c(input$metalsTS_mat1,input$metalsTS_mat2,input$metalsTS_mat3)
	matrices=c("","","")
	matrices[matrices_num==1]="Water"
	matrices[matrices_num==2]="Brine shrimp"
	matrices[matrices_num==3]="Bird eggs"

	xlimit=c(min(metalsData_plots$Date),max(metalsData_plots$Date))
	par(mfrow=c(3,1),mar=c(4.1,6.1,2.1,4.1))
	if(dim(metalsData_plots[metalsData_plots$Parameter==input$metalsTS_param1&metalsData_plots$Matrix==matrices[1],])[1]>0){
		if(input$log_metalsTS1=="TRUE"){logy1="y"}else{logy1=""}
		plotData1=metalsData_plots[metalsData_plots$Parameter==input$metalsTS_param1&metalsData_plots$Matrix==matrices[1],]
		ylabel=paste0(input$metalsTS_param1," ",unique(plotData1$Units))
		plot(Result~Date,plotData1,xlim=xlimit,pch=NA,cex.axis=2,cex.lab=2,xlab="",ylab=ylabel,log=logy1)
		points(Result~Date,plotData1[plotData1$WaterColLoc=="D",],col="Orange",pch=17,cex=2.5)
		points(Result~Date,plotData1[plotData1$WaterColLoc=="S",],col="Blue",pch=21,cex=2.5)
		if(matrices[1]=="Brine shrimp"){bi_col="green"}else{bi_col="purple"}
		points(Result~Date,plotData1[plotData1$WaterColLoc=="BIO",],col=bi_col,pch=19,cex=2.5)
		if(matrices[1]=="Water"){legend("topleft",bty='n',pch=c(21,17),cex=2.5,col=c("blue","orange"),legend=c("Surface","Deep"))}
		if(matrices[1]=="Brine shrimp"){legend("topleft",bty='n',pch=19,cex=2.5,col=c("green"),legend="Brine shrimp")}
		if(matrices[1]=="Bird eggs"){legend("topleft",bty='n',pch=19,cex=2.5,col=c("purple"),legend="Bird eggs")}
		}else{frame()}
	if(dim(metalsData_plots[metalsData_plots$Parameter==input$metalsTS_param2&metalsData_plots$Matrix==matrices[2],])[1]>0){
		if(input$log_metalsTS2=="TRUE"){logy2="y"}else{logy2=""}
		plotData2=metalsData_plots[metalsData_plots$Parameter==input$metalsTS_param2&metalsData_plots$Matrix==matrices[2],]
		ylabel=paste0(input$metalsTS_param2," ",unique(plotData2$Units))
		plot(Result~Date,plotData2,xlim=xlimit,pch=NA,cex.axis=2,cex.lab=2,xlab="",ylab=ylabel,log=logy2)
		points(Result~Date,plotData2[plotData2$WaterColLoc=="D",],col="Orange",pch=17,cex=2.5)
		points(Result~Date,plotData2[plotData2$WaterColLoc=="S",],col="Blue",pch=21,cex=2.5)
		if(matrices[2]=="Brine shrimp"){bi_col="green"}else{bi_col="purple"}
		points(Result~Date,plotData2[plotData2$WaterColLoc=="BIO",],col=bi_col,pch=19,cex=2.5)
		if(matrices[2]=="Water"){legend("topleft",bty='n',pch=c(21,17),cex=2.5,col=c("blue","orange"),legend=c("Surface","Deep"))}
		if(matrices[2]=="Brine shrimp"){legend("topleft",bty='n',pch=19,cex=2.5,col=c("green"),legend="Brine shrimp")}
		if(matrices[2]=="Bird eggs"){legend("topleft",bty='n',pch=19,cex=2.5,col=c("purple"),legend="Bird eggs")}
		}else{frame()}
	if(dim(metalsData_plots[metalsData_plots$Parameter==input$metalsTS_param3&metalsData_plots$Matrix==matrices[3],])[1]>0){
		if(input$log_metalsTS3=="TRUE"){logy3="y"}else{logy3=""}
		plotData3=metalsData_plots[metalsData_plots$Parameter==input$metalsTS_param3&metalsData_plots$Matrix==matrices[3],]
		ylabel=paste0(input$metalsTS_param3," ",unique(plotData3$Units))
		plot(Result~Date,plotData3,xlim=xlimit,pch=NA,cex.axis=2,cex.lab=2,xlab="",ylab=ylabel,log=logy3)
		points(Result~Date,plotData3[plotData3$WaterColLoc=="D",],col="Orange",pch=17,cex=2.5)
		points(Result~Date,plotData3[plotData3$WaterColLoc=="S",],col="Blue",pch=21,cex=2.5)
		if(matrices[3]=="Brine shrimp"){bi_col="green"}else{bi_col="purple"}
		points(Result~Date,plotData3[plotData3$WaterColLoc=="BIO",],col=bi_col,pch=19,cex=2.5)
		if(matrices[3]=="Water"){legend("topleft",bty='n',pch=c(21,17),cex=2.5,col=c("blue","orange"),legend=c("Surface","Deep"))}
		if(matrices[3]=="Brine shrimp"){legend("topleft",bty='n',pch=19,cex=2.5,col=c("green"),legend="Brine shrimp")}
		if(matrices[3]=="Bird eggs"){legend("topleft",bty='n',pch=19,cex=2.5,col=c("purple"),legend="Bird eggs")}
		}else{frame()}
	}
})


output$plot3=renderPlot({

	if(input$metalsPlotType==2){

	#Subset by year input
	metalsData_plots=metalsData[metalsData$Year>=input$years_metals[1]&metalsData$Year<=input$years_metals[2],]
	metalsData_plots=metalsData[metalsData$Year>=input$years_metals[1]&metalsData$Year<=input$years_metals[2],]
	
	#Subset by surface/deep input
	surf_deep_metals=input$surf_deep_metals
	surf_deep_metals[surf_deep_metals=="Surface"]="S"
	surf_deep_metals[surf_deep_metals=="Deep"]="D"
	surf_deep_metals=append(surf_deep_metals,rep("BIO",3-length(surf_deep_metals)))
	metalsData_plots=metalsData_plots[metalsData_plots$WaterColLoc==surf_deep_metals[1]|metalsData_plots$WaterColLoc==surf_deep_metals[2]|metalsData_plots$WaterColLoc==surf_deep_metals[3],]

	metalsData_plots=(metalsData_plots)
		par(mar=c(6.1,6.1,2.1,2.1))
		if(input$metalsBP_mat==1){bp_matrix="Water"}
		if(input$metalsBP_mat==2){bp_matrix="Brine shrimp"}
		if(input$metalsBP_mat==3){bp_matrix="Bird eggs"}
		if(input$log_metalsBP=="TRUE"){logybp="y"}else{logybp=""}
		bp_plotdata=metalsData_plots[metalsData_plots$Parameter==input$metalsBP_param&metalsData_plots$Matrix==bp_matrix,]
		bp_plotdata=droplevels(bp_plotdata)
		ylabel=paste0(input$metalsBP_param," ",unique(bp_plotdata$Units))
		if(dim(bp_plotdata)[1]>0){boxplot(bp_plotdata$Result~bp_plotdata[,input$metalsBP_aggby],log=logybp,ylab=ylabel,cex.lab=1.25,cex.axis=1.25)}
	}
})


output$table3 <- renderUI({
                    output$tableTemp <- renderTable(makeMetalsSummaryTable(metalsData,input$years_metals,input$surf_deep_metals,input$metalsBP_mat,input$metalsBP_param,input$metalsBP_aggby),rownames=TRUE)
                    tableOutput("tableTemp")
	})



source("functions/ipdwMap.R",local=TRUE)
cost=raster("polyrast/costraster3")





output$wq_map <- renderPlot({
	if(input$interpolate>0){
		#input$interpolate
		isolate({
			uniform_scale=input$wq_map_uniform_scale
			parameter=input$wq_map_param
			depths=input$wq_map_depths
			dates=input$wq_map_dates
			date_min=dates[1]
			date_max=dates[2]
			maptitle=names(param_names[param_names==input$wq_map_param])
			ipdwMap(data=gsl_data,parameter,depths,date_min=date_min,date_max=date_max,costraster=cost,maptitle=maptitle,mask_poly=gsl_poly,uniform_scale=uniform_scale)
		})
	}
})


output$elev_plot <- renderPlot({
	par(mfrow=c(3,1),mar=c(4.1,6.1,2.1,4.1))
	elev_plot_data=combelv[combelv$Date>=input$elev_plot_dates[1]&combelv$Date<=input$elev_plot_dates[2],]
	xlimit=c(input$elev_plot_dates[1],input$elev_plot_dates[2])
	
	#Elevation plot
	plot(value~as.Date(Date),elev_plot_data[elev_plot_data$variable=="SA_LakeElevFt"|elev_plot_data$variable=="NA_LakeElevFt",],pch=NA,xlab="",ylab="Elevation (ft)",cex.axis=2,cex.lab=2.25,xlim=xlimit)
	points(value~as.Date(Date),na.omit(elev_plot_data[elev_plot_data$variable=="SA_LakeElevFt",]),type='l',col="blue",lwd=3)
	points(value~as.Date(Date),na.omit(elev_plot_data[elev_plot_data$variable=="NA_LakeElevFt",]),type='l',col="green",lwd=3,lty=3)
	if(input$elev_legend==1){
		legend("topright",bty='n',lwd=3,col=c("blue","green"),lty=c(1,3),legend=c("Gilbert Bay (south arm)","Gunnison Bay (north arm)"),cex=1.75)
	}
	if(input$elev_rulers==1){
		abline(h=input$elev_ruler1,lwd=3,lty=2,col="purple")
	}
	
	#Volume plot
	plot(value/1000~as.Date(Date),elev_plot_data[elev_plot_data$variable=="SA_vol_acreft"|elev_plot_data$variable=="NA_vol_acreft",],pch=NA,xlab="",ylab="Volume (1,000 acre-ft)",cex.axis=2,cex.lab=2.25,xlim=xlimit)
	points(value/1000~as.Date(Date),na.omit(elev_plot_data[elev_plot_data$variable=="SA_vol_acreft",]),type='l',col="blue",lwd=3)
	points(value/1000~as.Date(Date),na.omit(elev_plot_data[elev_plot_data$variable=="NA_vol_acreft",]),type='l',col="green",lwd=3,lty=3)
	if(input$elev_legend==1){
		legend("topright",bty='n',lwd=3,col=c("blue","green"),lty=c(1,3),legend=c("Gilbert Bay (south arm)","Gunnison Bay (north arm)"),cex=1.75)
	}
	if(input$elev_rulers==1){
		abline(h=input$elev_ruler2,lwd=3,lty=2,col="purple")
	}


	#Head plot
	if(min(na.omit(elev_plot_data[elev_plot_data$variable=="SmNhead","value"]))>0){
		plot(value~as.Date(Date),na.omit(elev_plot_data[elev_plot_data$variable=="SmNhead",]),type='l',xlab="",ylab="S-N head (ft)",cex.axis=2,cex.lab=2.25,col="orange",xlim=xlimit,lwd=3,ylim=c(0,max(na.omit(elev_plot_data[elev_plot_data$variable=="SmNhead","value"]))*1.1))
	}else{
		plot(value~as.Date(Date),na.omit(elev_plot_data[elev_plot_data$variable=="SmNhead",]),type='l',xlab="",ylab="S-N head (ft)",cex.axis=2,cex.lab=2.25,col="orange",xlim=xlimit,lwd=3)
	}
	abline(0,0,lty=2,lwd=3)
	if(input$elev_rulers==1){
		abline(h=input$elev_ruler3,lwd=3,lty=2,col="purple")
	}


})

}






shinyApp(ui = ui, server = server)







