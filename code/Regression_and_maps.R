rm(list=ls())
options(stringsAsFactors=F)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


library(dplyr)
library(relaimpo)
library(tidyverse)
library(raster)
library(car)

# read in various functions for cleaning/plotting data
source("misc_functions.R")

# read in the data
comb<-read_csv("../fungi_data/Fungal_climate_data.csv") %>% data.frame()

# scale the variables to be within [0,1]
comb$temp01<-scale01(comb$temp.niche.width)
comb$water01<-scale01(comb$water.niche.width)
comb$tw01<-scale01(comb$temp01*comb$water01)

# create the trade-off variables
comb <- comb %>% mutate(dom_tol=ranking-tw01)  %>% 
				 mutate(tdom_tol=ranking-temp01)  %>% 
				 mutate(wdom_tol=ranking-water01)  %>% 
				 mutate(water01=sqrt(water01)) %>% 
				 mutate(temp01=sqrt(temp01)) %>% 	
				 mutate(ranking=sqrt(ranking))  

# main regression models
summary(fit1.s<-lm(ranking~     PC2.p + PC3.p +  PC2.t +    a.gal + a.tab + m.tre + p.fla + p.ruf,data=comb))
summary(fit2.s<-lm(dom_tol~     PC2.p + PC3.p +  PC2.t +    a.gal +         m.tre + p.fla + p.ruf,data=comb))
summary(fit3.s<-lm(tdom_tol~    PC2.p + PC3.p +  PC2.t +    a.gal +                 p.fla + p.ruf,data=comb))
summary(fit4.s<-lm(wdom_tol~    PC2.p + PC3.p +  PC2.t +    a.gal +         m.tre + p.fla + p.ruf,data=comb))
summary(fit5.s<-lm(temp01~      PC2.p + PC3.p +  PC2.t +    a.gal                                ,data=comb))
summary(fit6.s<-lm(water01~     PC2.p + PC3.p +  PC2.t +    a.gal                                ,data=comb))


# fit corresponding null models (keeping species indictors) to quantify relative importance of environmental variables
summary(fit10<-lm(ranking~                                  a.gal + a.tab + m.tre + p.fla + p.ruf,data=comb))
summary(fit20<-lm(dom_tol~                                  a.gal +         m.tre + p.fla + p.ruf,data=comb))
summary(fit30<-lm(tdom_tol~                                 a.gal +                 p.fla + p.ruf,data=comb))
summary(fit40<-lm(wdom_tol~                                 a.gal +         m.tre + p.fla + p.ruf,data=comb))
summary(fit50<-lm(temp01~                                   a.gal                                ,data=comb))
summary(fit60<-lm(water01~                                  a.gal                                ,data=comb))


library(lmtest)
## overall test of environmental variables
print(l1<-anova(fit1.s,fit10))
print(l2<-anova(fit2.s,fit20))
print(l3<-anova(fit3.s,fit30))
print(l4<-anova(fit4.s,fit40))
print(l5<-anova(fit5.s,fit50))
print(l6<-anova(fit6.s,fit60))


# confidence intervals
confint(fit2.s)
confint(fit3.s)
confint(fit4.s)
confint(fit1.s)
confint(fit5.s)
confint(fit6.s)


# get adjusted R2 attributable to the climate PCAs
relimp_group(fit2.s,c("PC2.p","PC3.p","PC2.t"))
relimp_group(fit3.s,c("PC2.p","PC3.p","PC2.t"))
relimp_group(fit4.s,c("PC2.p","PC3.p","PC2.t"))
relimp_group(fit1.s,c("PC2.p","PC3.p","PC2.t"))
relimp_group(fit5.s,c("PC2.p","PC3.p","PC2.t"))
relimp_group(fit6.s,c("PC2.p","PC3.p","PC2.t"))


# model AIC
AIC(fit2.s,fit3.s,fit4.s)
AIC(fit1.s,fit5.s,fit6.s)

#################################################3
########################################################
################## phylogenetic regression
# libraries
library(ape)
library(caper)
library(geiger)
library(phytools)
library(phangorn)
library(nlme)
library(dplyr)


comb.scale<-comb.scale.o<-comb


# get tree (shortcut)
tree <- read.nexus("../fungi_data/LSU_phylogenetic_tree_rooted.nexus")
is.ultrametric(tree)
tree.fungi<-drop.tip(tree, tree$tip.label[!tree$tip.label%in%comb$name3])
rownames(comb.scale.o) <- comb.scale.o$name3

plot(tree.fungi, cex=1)  


# match species in data file and tree
row.names(comb.scale.o) <- comb.scale.o$gen.name2  # make species names row names
name.check(tree.fungi, comb.scale.o)


spec_inds<-c("a.gal","a.tab","h.set","m.tre","p.fla","p.rob","p.ruf","s.com")


set.seed(20)

### fit the models allowing lambda to vary
summary(g1 <- gls(trim_formula(fit1.s,spec_inds), correlation=corPagel(value=1, phy=tree.fungi), data=comb.scale.o,method='REML'))
summary(g2 <- gls(trim_formula(fit2.s,spec_inds), correlation=corPagel(value=1, phy=tree.fungi), data=comb.scale.o,method='REML'))
summary(g3 <- gls(trim_formula(fit3.s,spec_inds), correlation=corPagel(value=1, phy=tree.fungi), data=comb.scale.o,method='REML'))
summary(g4 <- gls(trim_formula(fit4.s,spec_inds), correlation=corPagel(value=1, phy=tree.fungi), data=comb.scale.o,method='REML'))
summary(g5 <- gls(trim_formula(fit5.s,spec_inds), correlation=corPagel(value=1, phy=tree.fungi), data=comb.scale.o,method='REML'))
summary(g6 <- gls(trim_formula(fit6.s,spec_inds), correlation=corPagel(value=1, phy=tree.fungi), data=comb.scale.o,method='REML'))

# force lambda to be zero
g1.lam <- gls(trim_formula(fit1.s,spec_inds), correlation=corPagel(value=0, phy=tree.fungi,fixed=T), data=comb.scale.o,method='REML')
g2.lam <- gls(trim_formula(fit2.s,spec_inds), correlation=corPagel(value=0, phy=tree.fungi,fixed=T), data=comb.scale.o,method='REML')
g3.lam <- gls(trim_formula(fit3.s,spec_inds), correlation=corPagel(value=0, phy=tree.fungi,fixed=T), data=comb.scale.o,method='REML')
g4.lam <- gls(trim_formula(fit4.s,spec_inds), correlation=corPagel(value=0, phy=tree.fungi,fixed=T), data=comb.scale.o,method='REML')
g5.lam <- gls(trim_formula(fit5.s,spec_inds), correlation=corPagel(value=0, phy=tree.fungi,fixed=T), data=comb.scale.o,method='REML')
g6.lam <- gls(trim_formula(fit6.s,spec_inds), correlation=corPagel(value=0, phy=tree.fungi,fixed=T), data=comb.scale.o,method='REML')


# test of lambda=0, can use reml since the fixed parts are the same
print(gt2<-anova(g2,g2.lam))
print(gt3<-anova(g3,g3.lam))
print(gt4<-anova(g4,g4.lam))
print(gt1<-anova(g1,g1.lam))
print(gt5<-anova(g5,g5.lam))
print(gt6<-anova(g6,g6.lam))

# view lambda
as.numeric(attr(g2$apVar,"Pars")[1])
as.numeric(attr(g3$apVar,"Pars")[1])
as.numeric(attr(g4$apVar,"Pars")[1])
as.numeric(attr(g1$apVar,"Pars")[1])
as.numeric(attr(g5$apVar,"Pars")[1])
as.numeric(attr(g6$apVar,"Pars")[1])

# get AIC
AIC(update(g2,method="ML"),fit2.s)
AIC(update(g3,method="ML"),fit3.s)
AIC(update(g4,method="ML"),fit4.s)
AIC(update(g1,method="ML"),fit1.s)
AIC(update(g5,method="ML"),fit5.s) # note, ML won't converge due to < lambda
AIC(update(g6,method="ML"),fit6.s)



##############################################
#########           MAKING MAPS
##############################################

library(sp)
library(ggfortify)
library(ggpolypath)
library(PBSmapping)
library(ggplot2)
library(colorRamps)
library(cowplot)
library(colorspace)
library(viridis)

# set the bounds
lat_lim<-c(0,70)
lon_lim<-c(-155,-50)

# read in climate data
bstack.o<-brick("../geo_data/biome_stack.tif")
lake<-shapefile("../geo_data/lakes/ne_10m_lakes.shp")
ocean<-shapefile("../geo_data/oceans/ne_10m_ocean.shp")
kc0 <- raster("../geo_data/koppen_climate.grd")
forest_mask<-raster("../geo_data/Forest_mask_layer.tif")
climate.pca<-read_csv("../geo_data/bioclim_PCA_variables.csv") %>% mutate(x=long,y=lat) %>% data.frame() 

# clean up the layers
lake_fortify<-fortify(lake) %>% dplyr::rename(X=long, Y=lat, PID=group, POS=order) %>% mutate(piece=as.numeric(piece),PID=as.numeric(PID))
ocean_fortify<-fortify(ocean) %>% dplyr::rename(X=long, Y=lat, PID=group, POS=order) %>% mutate(piece=as.numeric(piece),PID=as.numeric(PID))
clim_ord <- c('Af', 'Am', 'As', 'Aw', 'BSh', 'BSk', 'BWh', 'BWk', 'Cfa', 'Cfb','Cfc', 'Csa', 'Csb', 'Csc', 'Cwa','Cwb', 'Cwc', 'Dfa', 'Dfb', 'Dfc','Dfd', 'Dsa', 'Dsb', 'Dsc', 'Dsd','Dwa', 'Dwb', 'Dwc', 'Dwd', 'EF','ET', 'Ocean')
use_clim<-c("Aw","Am","Af","Cfa","Cfb","Dfa","Dfb")
keep_clim<-(1:length(clim_ord))[clim_ord%in%use_clim]



############ making maps

# main temperature plot
gg3<-plot_fig(mod=fit3.s,climate.pca=climate.pca, bstack.o=bstack.o, kc0=kc0, forest_mask=forest_mask, lake_fortify=lake_fortify, ocean_fortify=ocean_fortify,scale_func=scale11)
gg3av<-av_gg(fit3.s,"PC3.p",xlab="Climate PCA 1")
plot_grid(gg3av,gg3,rel_widths=c(.7,1),labels="AUTO") 


# ranking and temp niche width
gg1<-plot_fig(mod=fit1.s,climate.pca=climate.pca, bstack.o=bstack.o, kc0=kc0, lake_fortify=lake_fortify, ocean_fortify=ocean_fortify)+ggtitle("A. Competitive ability")
gg5<-plot_fig(mod=fit5.s, climate.pca=climate.pca, bstack.o=bstack.o, kc0=kc0, lake_fortify=lake_fortify, ocean_fortify=ocean_fortify,flip=T)+ggtitle("B. Thermal niche width")


# temp x moist and moist dom tol
gg2<-plot_fig(mod=fit2.s, climate.pca=climate.pca, bstack.o=bstack.o, kc0=kc0, lake_fortify=lake_fortify, ocean_fortify=ocean_fortify)+ggtitle("A. Combined trade-off")
gg4<-plot_fig(mod=fit4.s, climate.pca=climate.pca, bstack.o=bstack.o, kc0=kc0, lake_fortify=lake_fortify, ocean_fortify=ocean_fortify)+ggtitle("B. Moisture trade-off")

