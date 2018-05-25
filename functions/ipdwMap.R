
#ipdwMap function
ipdwMap=function(data,parameter,date_min,date_max,depths,costraster,maptitle,mask_poly,uniform_scale=FALSE){
	mapdata=na.omit(data[data$Depth%in%depths&data$Date>=date_min&data$Date<=date_max,c("Site","LatDD","LongDD",parameter)])
		
	if(length(depths)==1){
		maptitle_date_depth=paste0(maptitle," (",date_min," to ",date_max,")"," (",depths,")")
	}
	if(length(depths)==2){
		maptitle_date_depth=paste0(maptitle," (",date_min," : ",date_max,")"," (",depths[1]," & ",depths[2],")")
	}
	if(dim(mapdata)[1]>0){
		withProgress(message="Initializing...",{
			aggmean=aggregate(get(parameter)~Site+LatDD+LongDD,data=mapdata,FUN='mean')
			#aggmean_all=aggregate(get(parameter)~Site+LatDD+LongDD,data=data,FUN='mean')
			aggcount=aggregate(get(parameter)~Site+LatDD+LongDD,data=mapdata,FUN='length')
			names(aggmean)[names(aggmean)=="get(parameter)"]=parameter
			#names(aggmean_all)[names(aggmean_all)=="get(parameter)"]=parameter
			names(aggcount)[names(aggcount)=="get(parameter)"]="n_count"
			dsp <- SpatialPoints(aggmean[,3:2], proj4string=CRS("+proj=longlat +datum=WGS84"))
			dsp <- SpatialPointsDataFrame(dsp, aggmean)
			dsp=merge(dsp,aggcount)
			projection(costraster)=projection(dsp)
			
			setProgress(message="Interpolating...",0.3)
			dsp.ipdw=ipdw(dsp,costraster,paramlist=parameter,range=10)
			#x=pathdistGen(dsp,costraster,range=10)
			
			setProgress(message="Plotting...",0.6)

			maskr=rasterize(mask_poly,dsp.ipdw)
			dsp.ipdw.mask=mask(x=dsp.ipdw,maskr)
			
			setProgress(message="Plotting...",0.9)
			
			projection(mask_poly)=projection(dsp)
			plot(mask_poly,col=NA,border=NA,main=maptitle_date_depth,cex.main=1.75)
			
			setProgress(1)
			
			if(uniform_scale==TRUE){
				breaks=unique(signif(seq(min(data[,parameter],na.rm=TRUE),max(data[,parameter],na.rm=TRUE),length.out=9),2))
				mappalette=brewer.pal(n = length(breaks), name = "OrRd")
				plot(dsp.ipdw.mask,axis.args=c(cex.axis=1.5),col=mappalette,add=T,breaks=breaks)
				#plot(costraster)
			}else{
				mappalette=brewer.pal(n = 9, name = "OrRd")
				plot(dsp.ipdw.mask,axis.args=c(cex.axis=1.5),col=mappalette,add=T)
				#plot(costraster)
			}
			text(dsp,dsp$n_count,cex=1.75,halo=TRUE,hw=0.2)
			axis(1,cex.axis=1.5)
			axis(2,cex.axis=1.5)

		})
	}else{
		frame()
		text(0.5,0.5,"No data for selected parameters.",cex=2.25)
	}
}
#ipdwMap(gsl_data,"Total_nitrogen_mgL",c("Surface","Deep"),date_min="2001-08-01",date_max="2018-04-01")


