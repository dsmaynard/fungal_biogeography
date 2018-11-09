

rm(list=ls())
options(stringsAsFactors=F)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


library(tidyverse)
library(circlize)
library(ggbiplot)
library(ggfortify)
library(corrplot)

# source the functions
source("misc_functions.R")


# read in trait data
trait.dat<-read.csv("../fungi_data/Fungal_trait_data.csv")


##### PCA of trait trade-off

# get the relevant variables, identified in clustering
performance.vars<-c("ranking","temp.temp.at.max.rate","water.mpa.at.max.r","temp.niche.width","water.niche.width","rate.0.5","density")
enz.vars<-c("cbh.w","phos.w","perox","abts.w")
vol.vars<-c("X11","X35","X7","X9","X55","X25","X16","X48","X87","X73","X81")

# remove NAs and add species labels
use.dat<-trait.dat[,c(performance.vars,enz.vars,vol.vars)] 
rownames(use.dat) <- trait.dat$name3
use.dat<-use.dat[complete.cases(use.dat),]
use.dat$species <- substr(rownames(use.dat),1,5)


# get PCA
autoplot(prcomp(use.dat %>% select(-species),scale=T),loadings=T,label=F,loadings.colour = 'black', size=0,alpha=0,
		 loadings.label = TRUE,loadings.label.colour="black",loadings.label.repel=T) + xlim(-.4, .35)+ ylim(-.35,.55)+theme_bw()+
	theme(legend.position = "none", axis.title = element_text(size=16),
		  panel.border = element_rect(colour = "black", fill=NA, size=1))+
    geom_point(aes(fill=use.dat$species),size=4,pch=21,alpha=0.6)+
	xlab("PCA Axis 1 (22.3%)")+
	ylab("PCA Axis 2 (11.4%)")
	




#################################
# trait circle plot

use.dat<-trait.dat[,c(performance.vars,enz.vars,vol.vars)] 


# get all correlations
cor.all<-cor(use.dat,use="pairwise.complete.obs",method="kendall")

# print all, note that only those niche traits significant after multiple comparisons are shown in main text
corrplot(cor.all,method="circle",type="upper",tl.col="black",diag=F,tl.cex=1)

