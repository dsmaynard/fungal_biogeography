

rm(list=ls())
options(stringsAsFactors=F)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


library(tidyverse)
library(circlize)
library(ggbiplot)


# source the functions
source("misc_functions.R")


# read in trait data
trait.dat<-read.csv("../fungi_data/Fungal_trait_data.csv")


##### PCA of trait trade-off

# get the relevant variables, identified in clustering
performance.vars<-c("ranking","temp.temp.at.max.rate","water.mpa.at.max.r","temp.niche.width","water.niche.width","rate.0.5","density")
enz.vars<-c("cbh.w","phos.w","perox","abts.w")
vol.vars<-c("X11","X35","X7","X9","X55","X25","X16","X48","X87","X73","X81")


# remove NAs 
use.dat<-trait.dat[,c(performance.vars,enz.vars,vol.vars)] 
use.dat<-use.dat[complete.cases(use.dat),]

# add in the print names
show_names<-read_csv("name_list_shortened.csv")
names(use.dat)<-show_names$show.name[match(names(use.dat),show_names$var.name)]

# get PCA
trait.pca <- prcomp(use.dat,scale. = TRUE)
trait.pca$rotation<--trait.pca$rotation
scores<-trait.pca$x
par(mar=c(0,0,0,0))
g <- ggbiplot::ggbiplot(trait.pca, obs.scale = 1, var.scale = 1)+xlim(-3.5, 3)+ylim(-2.5,3)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
legend.position = 'top')
g <- g + geom_point(size=3,pch=23)
g # note that this can be flipped and rotated




#################################
# trait circle plot

tr.mat<-trait.dat[,c("temp.niche.width","temp.temp.at.max.rate","water.niche.width",
				"water.mpa.at.max.r","ranking","rate.0.5","density")]


cor.all<-cor(tr.mat,use="pairwise.complete.obs",method="kendall")

diag(cor.all)<-0
cor.all<-as.matrix(cor.all)
cor.all[lower.tri(cor.all)]<-0

cor.val<-cor.p<-matrix(NA,nrow=ncol(tr.mat),ncol=ncol(tr.mat))
rownames(cor.val)<-rownames(cor.p)<-colnames(cor.val)<-colnames(cor.p)<-colnames(tr.mat)
col.mat<-matrix(NA,nrow=ncol(tr.mat),ncol=ncol(tr.mat))


# cycle through and get the correlations, assigning colors as well
colv<-colorRampPalette(c("red4","white","navyblue"))
colv<-colv(200)
for(i in 1:ncol(tr.mat)){
	for(j in i:ncol(tr.mat)){
		if(j>i){
			tv<-cor.test(as.numeric(tr.mat[,i]),as.numeric(tr.mat[,j]),na.action="na.omit",method="kendall",exact=F)
			cor.val[i,j]<-tv$estimate
			cor.p[i,j]<-tv$p.value
			col.mat[i,j]<-scales::alpha(colv[round(tv$estimate*100)+100],.8)		
		}
	}
}


# add the display names
show_names$show.name<-gsub("_","\\\n",show_names$show.name)
rownames(col.mat)<-colnames(col.mat)<-rownames(cor.val)
rownames(cor.val)<-show_names$show.name[match(rownames(cor.val),show_names$var.name)]
colnames(cor.val)<-show_names$show.name[match(colnames(cor.val),show_names$var.name)]
rownames(cor.p)<-show_names$show.name[match(rownames(cor.p),show_names$var.name)]
colnames(cor.p)<-show_names$show.name[match(colnames(cor.p),show_names$var.name)]


# remove p values > 0.05
cor.p[cor.p>0.05]<-NA
cor.val[cor.p>0.05 | is.na(cor.p)]<-0

# add missing values to correlation and color matrix
cor.val[cor.val==0]<-NA
col.mat[cor.val==0]<-NA
col.mat[is.na(cor.val)]<-NA
cor.p[is.na(cor.val)]<-NA

# convert p value matrix, correlation, and color matrix to df
cor.p2<-convert_to_df(cor.p)
cor.val2<-convert_to_df(cor.val)
col.val2<-convert_to_df(col.mat)

# adjust p values
cor.p2[,3]<-p.adjust(as.numeric(cor.p2[,3]),method="hochberg")

# choose only those significant
cor.p2[cor.p2[,3]>0.05,3]<-NA
cor.p2[,3]<-round(cor.p2[,3],4)
cor.p2


#plot
chordDiagram(cor.val2,grid.col="black",grid.border="black",col=col.val2[,3],self.link=1)


