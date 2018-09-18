############################
################ raster 


# do the calculation for the regression
calc_raster<-function(x,na.rm=TRUE){
	# remove missing value or if we are masking
	if(any(is.na(x))){
		return(NA)
	}
	else{
		# remove the mask layer
		x<-x[-1]
		return(c(x)%*%mod_coefs)
	}
}


# create the climate mask. return 1 if we are within range of all variables, NA otherwise
mask_climate<-function(x,na.rm=T){
	
	if(sum(is.na(x))>0){
		return(NA)
	}
	
	
	if(all(x>minvals_all) & all(x<maxvals_all)){
		return(1)
	}	
	else{
		return(NA)
	}
}

mask_forest<-function(x,na.rm=T){
	
	if(any(x>min_prop_forest)){
		return(1)
	}	
	else{
		return(NA)
	}
}

mask_both<-function(x,na.rm=T){
	if(any(is.na(x))){
		return(NA)
	}	
	else{
		return(1)
	}
}


# return a raster with NA for points with climate outside of the range of the masking variables
get_climate_mask<-function(mask_vars,raster_layers,prop_extrap,edat){

	min1<-apply(edat[,mask_vars],2,min)
	max1<-apply(edat[,mask_vars],2,max)
	minvals_all<<-min1-(max1-min1)*prop_extrap
	maxvals_all<<-max1+(max1-min1)*prop_extrap
	# get the indices
	inds<-as.numeric(substr(bio_names$layer[match(names(minvals_all),bio_names$name)],4,nchar(bio_names$layer[match(names(minvals_all),bio_names$name)])))
	rsub<-raster_layers[[inds]]
	
	# make the globally masked layer
	mlayer<-stackApply(rsub, indices=rep(1,length(minvals_all)), fun=mask_climate)
	
	return(mlayer)
}

plot_raster<-function(show_plot,col_pallette,leg_title){
	g0<-ggplot(show_plot %>% mutate(index_1=index_1-min(index_1,na.rm=T)) %>% mutate(index_1=index_1/max(index_1,na.rm=T)),
									aes(x=x, y=y)) + geom_tile(aes(fill = index_1)) + coord_equal()+
	geom_polypath(data=lake_clip, aes(x=X, y=Y, group=PID,fill=POS), inherit.aes = F,fill="lightskyblue1")+
	geom_polypath(data=ocean_clip, aes(x=X, y=Y, group=PID,fill=POS), inherit.aes = F,fill="lightskyblue1")	+
	scale_fill_gradientn(colors=rev(col_pallette(100)),name=leg_title)+
	theme(panel.background = element_rect(fill = "white"),panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
		  axis.text=element_text(size=8),
          axis.title=element_blank())+

	return(g0)
}


###################################################################

get.pca<-function(ttdat){

	use.dat<-use.dat.o<-ttdat
	use.dat<-use.dat[complete.cases(use.dat),]
	use.dat<-use.dat[,apply(use.dat,2,var)>0]
	wine.pca <- prcomp(use.dat,scale. = TRUE)
	tpca<-wine.pca$x[,1]
	tpca2<-rep(NA,nrow(use.dat.o))
	names(tpca2)<-rownames(use.dat.o)
	tpca2<-tpca[match(names(tpca2),names(tpca))]
	names(tpca2)<-rownames(use.dat.o)
	tpca2
}


convert_to_pca<-function(x,na.rm=T){
	
	res<-as.numeric(t(t(as.matrix(pr$rotation))%*%as.matrix((x-pr$center)/pr$scale)))
	
	res[res>prmax]<-NA
	res[res<prmin]<-NA
	
	return(res)
}



# for chord diagram
convert_to_df<-function(cor.val){

	nms<-rownames(cor.val)

	k<-length(nms)

	outp<-data.frame(matrix(0,nrow=k^2,ncol=3))

	count<-0
	for(i in 1:k){
		for(j in 1:k){
			if(j>i){
				count<-count+1
				outp[count,1]<-rownames(cor.val)[i]
				outp[count,2]<-colnames(cor.val)[j]
				outp[count,3]<-cor.val[i,j]
			}
		}
	}

	outp<-outp[outp[,3]!=0,]
	outp<-outp[!is.na(outp[,3]),]

	return(outp)
}




### calc the normalized, adjusted R2 attributable to the climate PCAs
relimp_group<-function(fit,keep.env){
	
	if(length(coef(fit))==1 | sum(names(coef(fit))%in%keep.env)==0){
		val<-0
		names(val)<-"G1"
	}
	else{
		if(sum(names(coef(fit))%in%keep.env)==1){
			if(length(names(coef(fit)))==2){
				val<-summary(fit)$r.squared
				names(val)<-"G1"
			}
			else{
				ret<-calc.relimp(fit)$lmg
				val<-ret[names(ret)%in%keep.env]/summary(fit)$r.squared*summary(fit)$adj.r.squared
				names(val)<-"G1"
			}
		}
		else if(all(names(coef(fit))[-1]%in%keep.env)){
			ret<-calc.relimp(fit)$lmg
			val<-sum(ret)/summary(fit)$r.squared*summary(fit)$adj.r.squared
			names(val)<-"G1"
		}
		else{
			val<-calc.relimp(fit,groups=c(names(coef(fit))[names(coef(fit))%in%keep.env]))$lmg["G1"]/summary(fit)$r.squared*summary(fit)$adj.r.squared
		}
	}
	return(val)
}

scale01<-function(x,na.rm=T,flip=F){
	x<-x-min(x,na.rm=T)
	x<-(x/max(x,na.rm=T))
	if(flip){
		x<-(1-x)
	}
	return(x)
}

trim_formula<-function(mod,remove_vars,add_gal=F){
	
	form<-as.character(formula(mod))
	preds<-gsub(" ","",strsplit(form[3],"\\+")[[1]])

	preds<-preds[!preds%in%remove_vars]
	if(length(preds)==0){
		preds<-"1"
	}
	
	if(add_gal){
		preds<-c(preds,"a.gal")
	}
	return(formula(paste(form[[2]],"~",paste(preds,collapse="+"))))
}


### map plotting
plot_fig<-function(mod,climate.pca,lat_lim=c(18,50),lon_lim=c(-130,-62),bstack.o,mask_layer,kc0,forest_mask,lake_fortify,ocean_fortify,
				   minv=0.001,maxv=0.999, col_pallette=plasma,scale_func=scale01,flip=F,return_layer=F){
	
	spec_names<- c("a.gal","a.tab","h.set","m.tre","p.fla","p.rob","p.ruf","s.com")
	spec_inds<-matrix(0,nrow=nrow(climate.pca),ncol=length(spec_names))
	colnames(spec_inds)<-spec_names
	climate.pca<-data.frame(climate.pca,spec_inds)
	
	
	# predict values
	mod_pred<-predict(mod,newdata=climate.pca,se.fit=T)
	pred.vals<-data.frame(index_1=mod_pred$fit,x=climate.pca$x,y=climate.pca$y)
	se.vals<-data.frame(index_1=mod_pred$se.fit,x=climate.pca$x,y=climate.pca$y)
	
	# convert the predictions to a raster layer
	spg <- pred.vals
	coordinates(spg) <- ~ x + y
	gridded(spg) <- TRUE
	rpred <- crop(raster(spg),extent(lon_lim,lat_lim))
	
	
	spg <- se.vals
	coordinates(spg) <- ~ x + y
	gridded(spg) <- TRUE
	rse <- crop(raster(spg),extent(lon_lim,lat_lim))
	
	
	mask_layer <- crop(forest_mask,rpred)
	
	kc<-crop(kc0,rpred)
	
	mask_layer<-resample(mask_layer,kc)
	rpred<-resample(rpred,kc)
	rse<-resample(rse,kc)
	
	## mask the variable based on biome
	bstack<-crop(bstack.o[[c(1,2,3,4,5)]],rpred)
	
	stack_list<-stack(rpred,mask_layer,kc)

	rpred_mask<-stackApply(stack_list,indices = rep(1,dim(stack_list)[3]),function(x,na.rm=F){ifelse(any(is.na(x)),NA,ifelse(x[2]==0,NA,ifelse(x[3]%in%keep_clim,x[1],NA)))})

	# convert back to df
	rpred_crop<- as.data.frame(as(crop(rpred_mask,extent(c(lon_lim,lat_lim))), "SpatialPixelsDataFrame")) %>% filter(!is.na(index_1))
	# rse_crop<- as.data.frame(as(crop(rse_mask,extent(c(lon_lim,lat_lim))), "SpatialPixelsDataFrame"))
	
	
	#clip
	lake_clip<-clipPolys(lake_fortify,xlim=lon_lim,ylim=lat_lim)
	ocean_clip<-clipPolys(ocean_fortify,xlim=lon_lim,ylim=lat_lim)
	mask_layer<-crop(mask_layer,extent(c(lon_lim,lat_lim)))
	
	# mutate variables
	high_lim<-maxv
	low_lim<-minv
	rpred_mask2<-data.frame(rpred_crop) %>% 
						 	mutate(var="main",type="Temp. dominance/tolerance\ntrade-off") %>% 
						 	filter(index_1>quantile(index_1,probs=low_lim), index_1<quantile(index_1,probs=high_lim)) %>% 
						 	mutate(index_1=scale_func(index_1,flip=flip))
						 		   
	
	
	### ggplotting
	col_pallette<-col_pallette
	ratio_val<-1
	water_col<-"lightskyblue1"
	water_border<-"black"
	na_color<-"white"
	strip_bg_col<-"gray70"
	
	### MAP 1
	
	fig1_mask<-rbind(rpred_mask2)	
	if(return_layer){
		return(fig1_mask)
	}
	
	# if(!flip){
		gg1<-ggplot(fig1_mask,aes(x=x, y=y)) + geom_tile(aes(fill = index_1)) + 
			# geom_polypath(data=lake_clip, aes(x=X, y=Y, group=PID,fill=POS,size=0.01), inherit.aes = F,fill=water_col,col=water_border,size=0.1)+
			geom_polypath(data=ocean_clip, aes(x=X, y=Y, group=PID,fill=POS), inherit.aes = F,fill=water_col,col=water_border,size=0.1)	+
			scale_fill_gradientn(colors=rev(col_pallette(100)),name="",na.value = na_color)+
			scale_x_continuous(expand=c(0,0),name="Longitude")+
			scale_y_continuous(expand=c(0,0),name="Latitude")+	
			theme_bw()+
			theme(panel.ontop = FALSE,
				  panel.background = element_rect(fill = na_color),
		          panel.grid.major = element_blank(),
		          panel.grid.minor = element_blank(),
				  strip.text.y=element_blank(),
				  strip.background = element_rect(fill=strip_bg_col))
	# }
	# else{
	# 	gg1<-ggplot(fig1_mask,aes(x=x, y=y)) + geom_tile(aes(fill = index_1)) + 
	# 		# geom_polypath(data=lake_clip, aes(x=X, y=Y, group=PID,fill=POS,size=0.01), inherit.aes = F,fill=water_col,col=water_border,size=0.1)+
	# 		geom_polypath(data=ocean_clip, aes(x=X, y=Y, group=PID,fill=POS), inherit.aes = F,fill=water_col,col=water_border,size=0.1)	+
	# 		scale_fill_gradientn(colors=col_pallette(100),name="",na.value = na_color)+
	# 		scale_x_continuous(expand=c(0,0),name="Longitude")+
	# 		scale_y_continuous(expand=c(0,0),name="Latitude")+	
	# 		theme_bw()+
	# 		theme(panel.ontop = FALSE,
	# 			  panel.background = element_rect(fill = na_color),
	# 	          panel.grid.major = element_blank(),
	# 	          panel.grid.minor = element_blank(),
	# 			  strip.text.y=element_blank(),
	# 			  strip.background = element_rect(fill=strip_bg_col))	
	# }
	
	return(gg1)
}



# making effect size plot
av_gg<-function(mod,var="PC3.p",xlab="Climate PCA 3",ylab="Dominance-tolerance",col_pallette=plasma){

	a<-data.frame(avPlots(fit3.s,terms=var))
	names(a)<-c("x","y")
	a <- a %>% mutate(y=((scale01(y)-0.5)*2))
	gg2<-ggplot(a,aes(x=x,y=y,col=y))+geom_point()+
		geom_smooth(method=lm,aes(color=..y..), size=1.5) +
	  	scale_color_gradientn(colors=rev(col_pallette(100)),name="")+
		xlab(xlab)+ylab(ylab)+  	theme_bw()+
		theme(legend.position = "none")
	
	return(gg2)
}

scale11<-function(x,flip=F){
	
	if(all(x>0)){
		x<-scale01(x)
		x<-x-(0.5)*2
	}
	else{
		minx<-abs(min(x,na.rm=T))
		maxx<-abs(max(x,na.rm=T))
		x<-ifelse(x<0,(x/minx),(x/maxx))
	}
	
	if(flip){
		x<- -x
	}
	
	return(x)

}
	
	
iden<-function(x,scale_min=1,scale_max=1){
	
	return(x^(1/scale_max))
}
	

	