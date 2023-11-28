#Jenn Chakravorty, Elizabeth Bach
#Nachusa Grasslands bison exclosure plots 2014-2020
#ordination and vector analysis
#updated February 2023

rm(list=ls())
#load packages
library(tidyverse)
library(vegan)
library(ggplot2)

#Read in data, use "NachusaBisonExComm2014-20.csv"
data.community2<-tibble(read.csv(file.choose(),header=T, na.strings="NA"))

# As our data are abundance-based data, we have to calculate 
#Bray-Curtis distances between samples instead of Euclidean distances. So:

dis <- vegdist(data.community2[,-c(1:5,320)],method="bray")

#Use the adonis2() function to perform PERMANOVA (non-parametric, multivariate analysis of variance) on the communities of each sample
#Again, exclude any columns of metadata (here, columns 1-8)
#First run the full model, look for any interactions with factors
adonis2(dis ~ data.community2$community*data.community2$graze*data.community2$bison_year, data.community2[,-c(1:5,320)], permutations=999)

#no interactions, community and bison year main effects, simplify model to main effects only
adonis2(dis ~ data.community2$community+data.community2$graze+data.community2$bison_year, data.community2[,-c(1:5,320)], permutations=999)

#Df SumOfSqs      R2       F Pr(>F)    
#data.community2$community    3    8.984 0.27617 14.1614  0.001 ***
#  data.community2$graze        1    0.259 0.00796  1.2242  0.215    
#data.community2$bison_year   2    0.661 0.02033  1.5634  0.031 *  
#  Residual                   107   22.628 0.69555                   
#Total                      113   32.533 1.00000           

#Tukey's for which community groups are different
mds.dist<-metaMDSdist(decostand(data.community2[,-c(1:5,320)], "pa"),k=3, index="jaccard", autotransform=FALSE,na.rm=TRUE)
com.stat<-betadisper(mds.dist, data.community2$community, type="median")
com.stat
TukeyHSD(com.stat)

#Tukey's for which year groups are different
yr.stat<-betadisper(mds.dist, data.community2$bison_year, type="median")
yr.stat
TukeyHSD(yr.stat)

#points
mds.ab2<-metaMDS(decostand(data.community2[,-c(1:5,320)],"total" ),distance="bray", k=2,autotransform=FALSE, na.rm=TRUE)
mds.ab2
#stress=0.2150
mds.ab3<-metaMDS(decostand(data.community2[,-c(1:5,320)],"total" ),distance="bray", k=3,autotransform=FALSE, na.rm=TRUE)
mds.ab3
#stress=0.1591
mds.ab4<-metaMDS(decostand(data.community2[,-c(1:5,320)],"total" ),distance="bray", k=4,autotransform=FALSE, na.rm=TRUE)
mds.ab4
#stress=0.1261
mds.ab5<-metaMDS(decostand(data.community2[,-c(1:5,320)],"total" ),distance="bray", k=5,autotransform=FALSE, na.rm=TRUE)
mds.ab5
#stress=0.105

#3 dimensions is a strong stress reduction from 2, smaller improvement with 4, so go with that (also 4 extremely difficult to interpret)

#visualize ordination plot
#Used mds.ab3 generated above
MDS<-data.frame(scores(mds.ab3, choices = c(1,2,3), display = c("sites")))
community<-as.factor(data.community2$community)
bison_year<-as.factor(data.community2$bison_year)
graze<-as.factor(data.community2$graze)
exclosure.NMDS<-data.frame(MDS, community, bison_year, graze)

colors.community<-c(rgb(137,96,179,alpha=255, max=255),
                    rgb(86,174,108,alpha=255, max=255),
                    rgb(186,73,92,alpha=255, max=255),
                    rgb(176,145,59,alpha=255, max=255))
#Use grayscale for Natural Areas Association Journal
#over-riding the colors in the code below 
#If need to produce color figures for presentations, remove "+scale_color_grey()" from each of the graph panels

#Separate out into 3 panels, summarize community type into single point with error range
#Add centroid ellipses with points
#Figure 4
ggplot.NMDS.3pointCent<-function(XX,ZZ,COLORS){
  library(ggplot2)
  library(plyr)
  library(dplyr)
  library(gtable)
  library(gridExtra)
  MDS1<-data.frame(scores(XX,choices = c(1,2,3), display = c("sites")))$NMDS1
  MDS2<-data.frame(scores(XX, choices = c(1,2,3), display = c("sites")))$NMDS2
  MDS3<-data.frame(scores(XX,choices = c(1,2,3), display = c("sites")))$NMDS3
  Treatment<-ZZ
  
  NMDS<-data.frame(MDS1,MDS2,MDS3, Treatment)
  
  NMDS.mean=aggregate(NMDS[,1:3],list(group=Treatment),mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()
for(g in levels(NMDS$Treatment)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                                                   veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

X1<-ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = Treatment),size=4,alpha=0.75) +
  geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), show.legend=FALSE, size=2, linetype=5)+ylab("NMDS2")+xlab("NMDS1")+
  theme_bw()+theme(aspect.ratio=1)+scale_color_grey()+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),legend.position="none",axis.line=element_line(colour="black", size=2), panel.grid = element_blank(), axis.ticks=element_line(colour="black", size=2), axis.ticks.length=unit(0.5,"cm"))
X1
X1a<-ggplotGrob(X1)

df_ell2 <- data.frame()
for(g in levels(NMDS$Treatment)){
  df_ell2 <- rbind(df_ell2, cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                                                   veganCovEllipse(cov.wt(cbind(MDS1,MDS3),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS3)))))
                                ,group=g))
}
X2<-ggplot(data = NMDS, aes(MDS1, MDS3)) + geom_point(aes(color=Treatment), size=4,alpha=0.75)+
  geom_path(data=df_ell2, aes(x=MDS1, y=MDS3,colour=group), show.legend=FALSE, size=2, linetype=5)+ylab("NMDS3")+xlab("NMDS1")+
  theme_bw()+theme(aspect.ratio=1)+scale_color_grey()+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),legend.position="none",axis.line=element_line(colour="black", size=2), panel.grid = element_blank(), axis.ticks=element_line(colour="black", size=2), axis.ticks.length=unit(0.5,"cm"))
X2a<-ggplotGrob(X2)

df_ell3 <- data.frame()
for(g in levels(NMDS$Treatment)){
  df_ell3 <- rbind(df_ell3, cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                                                   veganCovEllipse(cov.wt(cbind(MDS3,MDS2),wt=rep(1/length(MDS3),length(MDS3)))$cov,center=c(mean(MDS3),mean(MDS2)))))
                                ,group=g))
}

X3<-ggplot(data = NMDS, aes(MDS3, MDS2)) + geom_point(aes(color=Treatment), size=4,alpha=0.75)+
  geom_path(data=df_ell3, aes(x=MDS3, y=MDS2,colour=group), size=2, linetype=5)+ylab("NMDS2")+xlab("NMDS3")+annotate("text", x=0.04, y=0.23, label="stress=0.159", size=6)+scale_x_continuous(breaks= c(-0.1, 0, 0.1))+
  theme_bw()+theme(aspect.ratio=1)+scale_color_grey(labels=c('new planting','old planting','native prairie','savanna'))+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),legend.title=element_text(size=20),legend.text=element_text(size=14),panel.grid = element_blank(), axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.line=element_line(colour="black", size=2), axis.ticks=element_line(colour="black", size=2), axis.ticks.length=unit(0.5,"cm"))
X3a<-ggplotGrob(X3)
print(grid::grid.newpage())
print(grid.arrange(X1a, X2a, X3a, nrow=1, widths=c(1.5,1.6,2.3)))
}

community.graph.cent<-ggplot.NMDS.3pointCent(mds.ab3, exclosure.NMDS$community, colors.community)
community.graph.cent

#add burn vector
#use "burn_byplotyear.csv"
exclosure.NMDS$siteID<-as.factor(data.community2$siteID)
exclosure.NMDS$year<-as.integer(data.community2$year)
burn.data<-tibble(read.csv(file.choose(),header=T, na.strings="NA"))
burn.data2<-tibble(burn.data %>% separate_wider_delim(plot_year, delim = "_", names = c("siteID","year")))
burn.data2$year<-as.integer(burn.data2$year)
burn.data3<-left_join(data.community2[,2:6], burn.data2, by=c("siteID","year","graze","community"))
exclosure.NMDSb<-left_join(exclosure.NMDS, burn.data3, by=c("community","graze","siteID","year"))

envectors_burn<-envfit(exclosure.NMDSb[,1:3] ~ exclosure.NMDSb$years_since_burn, choices=1:3, na.rm=TRUE)
envectors_burn

#This is right on the edge of significant
#we do share this result in the manuscript
#***VECTORS

#NMDS1    NMDS2    NMDS3     r2 Pr(>r)  
#exclosure.NMDSb$years_since_burn  0.24902 -0.47983  0.84128 0.0696  0.052 .
---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Permutation: free
#Number of permutations: 999

#Add vector to graph
scores(envectors_burn, display = "vectors")
env.scores<-tibble(x=c(0,0,0), y=c(-0.06570735,0.1266096,-0.2219831))
                
ggplot.NMDS.3pointCent2<-function(XX,ZZ,COLORS){
  library(ggplot2)
  library(plyr)
  library(dplyr)
  library(gtable)
  library(gridExtra)
  MDS1<-data.frame(scores(XX,choices = c(1,2,3), display = c("sites")))$NMDS1
  MDS2<-data.frame(scores(XX, choices = c(1,2,3), display = c("sites")))$NMDS2
  MDS3<-data.frame(scores(XX,choices = c(1,2,3), display = c("sites")))$NMDS3
  Treatment<-ZZ
  
  NMDS<-data.frame(MDS1,MDS2,MDS3, Treatment)
  
  NMDS.mean=aggregate(NMDS[,1:3],list(group=Treatment),mean)
  
  veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }
  
  df_ell <- data.frame()
  for(g in levels(NMDS$Treatment)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                                                     veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                  ,group=g))
  }
  
  X1<-ggplot() + geom_point(data = NMDS, aes(MDS1, MDS2, color = Treatment),size=4,alpha=0.75) + geom_path(aes(x=c(0,-0.0657), y=c(0,0.127)), size=2)+
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), show.legend=FALSE, size=2, linetype=5)+ylab("NMDS2")+xlab("NMDS1")+
    theme_bw()+theme(aspect.ratio=1)+scale_color_grey()+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),legend.position="none",axis.line=element_line(colour="black", size=2), panel.grid = element_blank(), axis.ticks=element_line(colour="black", size=2), axis.ticks.length=unit(0.5,"cm"))
  X1
  X1a<-ggplotGrob(X1)
  
  df_ell2 <- data.frame()
  for(g in levels(NMDS$Treatment)){
    df_ell2 <- rbind(df_ell2, cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                                                       veganCovEllipse(cov.wt(cbind(MDS1,MDS3),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS3)))))
                                    ,group=g))
  }
  X2<-ggplot() + geom_point(data = NMDS, aes(MDS1, MDS3, color=Treatment), size=4,alpha=0.75)+geom_path(aes(x=c(0,-0.0657), y=c(0,-0.222)), size=2)+
    geom_path(data=df_ell2, aes(x=MDS1, y=MDS3,colour=group), show.legend=FALSE, size=2, linetype=5)+ylab("NMDS3")+xlab("NMDS1")+
    theme_bw()+theme(aspect.ratio=1)+scale_color_grey()+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),legend.position="none",axis.line=element_line(colour="black", size=2), panel.grid = element_blank(), axis.ticks=element_line(colour="black", size=2), axis.ticks.length=unit(0.5,"cm"))
  X2a<-ggplotGrob(X2)
  
  df_ell3 <- data.frame()
  for(g in levels(NMDS$Treatment)){
    df_ell3 <- rbind(df_ell3, cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                                                       veganCovEllipse(cov.wt(cbind(MDS3,MDS2),wt=rep(1/length(MDS3),length(MDS3)))$cov,center=c(mean(MDS3),mean(MDS2)))))
                                    ,group=g))
  }
  
  X3<-ggplot() + geom_point(data = NMDS, aes(MDS3, MDS2, color=Treatment), size=4,alpha=0.75)+geom_path(aes(x=c(0,0.127), y=c(0,-0.222)), size=2)+
    geom_path(data=df_ell3, aes(x=MDS3, y=MDS2,colour=group), size=2, linetype=5)+ylab("NMDS2")+xlab("NMDS3")+annotate("text", x=0.04, y=0.23, label="stress=0.159", size=6)+scale_x_continuous(breaks= c(-0.1, 0, 0.1))+
    theme_bw()+theme(aspect.ratio=1)+scale_color_grey(labels=c('new planting','old planting','native prairie','savanna'))+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),legend.title=element_text(size=20),legend.text=element_text(size=14),panel.grid = element_blank(), axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.line=element_line(colour="black", size=2), axis.ticks=element_line(colour="black", size=2), axis.ticks.length=unit(0.5,"cm"))
  X3a<-ggplotGrob(X3)
  print(grid::grid.newpage())
  print(grid.arrange(X1a, X2a, X3a, nrow=1, widths=c(1.5,1.6,2.3)))
}



community.graph.cent2<-ggplot.NMDS.3pointCent2(mds.ab3, exclosure.NMDS$community, colors.community)
community.graph.cent2

####Supplimental Fig. 3 #######
#visualize differences among years
ggplot.NMDS.3pointCent2<-function(XX,ZZ,COLORS){
  library(ggplot2)
  library(plyr)
  library(dplyr)
  library(gtable)
  library(gridExtra)
  MDS1<-data.frame(scores(XX,choices = c(1,2,3), display = c("sites")))$NMDS1
  MDS2<-data.frame(scores(XX, choices = c(1,2,3), display = c("sites")))$NMDS2
  MDS3<-data.frame(scores(XX,choices = c(1,2,3), display = c("sites")))$NMDS3
  Treatment<-ZZ
  
  NMDS<-data.frame(MDS1,MDS2,MDS3, Treatment)
  
  NMDS.mean=aggregate(NMDS[,1:3],list(group=Treatment),mean)
  
  veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }
  
  df_ell <- data.frame()
  for(g in levels(NMDS$Treatment)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                                                     veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                  ,group=g))
  }
  
  X1<-ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = Treatment),size=4,alpha=0.75) +
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), show.legend=FALSE, size=2, linetype=5)+ylab("NMDS2")+xlab("NMDS1")+
    theme_bw()+theme(aspect.ratio=1)+scale_color_manual(values=COLORS)+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),panel.grid = element_blank(), axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),legend.position="none",axis.line=element_line(colour="black", size=2), axis.ticks=element_line(colour="black", size=2), axis.ticks.length=unit(0.5,"cm"))
  X1
  X1a<-ggplotGrob(X1)
  
  df_ell2 <- data.frame()
  for(g in levels(NMDS$Treatment)){
    df_ell2 <- rbind(df_ell2, cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                                                       veganCovEllipse(cov.wt(cbind(MDS1,MDS3),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS3)))))
                                    ,group=g))
  }
  X2<-ggplot(data = NMDS, aes(MDS1, MDS3)) + geom_point(aes(color=Treatment), size=4,alpha=0.75)+
    geom_path(data=df_ell2, aes(x=MDS1, y=MDS3,colour=group), show.legend=FALSE, size=2, linetype=5)+ylab("NMDS3")+xlab("NMDS1")+
    theme_bw()+theme(aspect.ratio=1)+scale_color_manual(values=COLORS)+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),panel.grid = element_blank(), axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),legend.position="none",axis.line=element_line(colour="black", size=2), axis.ticks=element_line(colour="black", size=2), axis.ticks.length=unit(0.5,"cm"))
  X2a<-ggplotGrob(X2)
  
  df_ell3 <- data.frame()
  for(g in levels(NMDS$Treatment)){
    df_ell3 <- rbind(df_ell3, cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                                                       veganCovEllipse(cov.wt(cbind(MDS3,MDS2),wt=rep(1/length(MDS3),length(MDS3)))$cov,center=c(mean(MDS3),mean(MDS2)))))
                                    ,group=g))
  }
  
  X3<-ggplot(data = NMDS, aes(MDS3, MDS2)) + geom_point(aes(color=Treatment), size=4,alpha=0.75)+
    geom_path(data=df_ell3, aes(x=MDS3, y=MDS2,colour=group), size=2, linetype=5)+ylab("NMDS2")+xlab("NMDS3")+annotate("text", x=0.04, y=0.23, label="stress=0.159", size=6)+scale_x_continuous(breaks= c(-0.1, 0, 0.1))+
    theme_bw()+theme(aspect.ratio=1)+scale_color_manual(values=COLORS, labels=c('0 years','3 years','5 years'))+theme(axis.text.x=element_text(size=20),panel.grid = element_blank(), axis.text.y=element_text(size=20),legend.title=element_text(size=20),legend.text=element_text(size=14),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.line=element_line(colour="black", size=2), axis.ticks=element_line(colour="black", size=2), axis.ticks.length=unit(0.5,"cm"))
  X3a<-ggplotGrob(X3)
  print(grid::grid.newpage())
  print(grid.arrange(X1a, X2a, X3a, nrow=1, widths=c(1.5,1.6,2.15)))
}

community.graph.cent2<-ggplot.NMDS.3pointCent2(mds.ab3, exclosure.NMDS$bison_year, colors.year)
community.graph.cent2


#####################################################
#indicator species analysis
library(indicspecies)
#summarize by species, use wide data
graze<-multipatt(data.community2[,-c(1:6)], cluster=data.community2$graze, func = "IndVal", duleg=TRUE)
summary(graze)
#Only significant indicator species for grazing is Symphyotrichum.sericeum is more probable in ungrazed plots
#what a cute indicator!

#break-down by community? Particularly for savanna, are there graze/ungraze indicators
data.community2$community<-as.factor(data.community2$community)
savanna<-droplevels(subset(data.community2, data.community2$community=="savanna"))
savanna.graze<-multipatt(savanna[,-c(1:6)], cluster=savanna$graze, func = "IndVal", duleg=TRUE)
summary(savanna.graze)
#SORNUT and MONFIST indicators for grazing, HELOCC indicator for ungrazed
#Group G  #sps.  2 
#stat p.value   
#Sorghastrum.nutans 0.934   0.005 **
#  Monarda.fistulosa  0.813   0.040 * 
#  
#  Group UG  #sps.  1 
#stat p.value  
#Helianthus.occidentalis 0.738    0.05 *

new.planting<-droplevels(subset(data.community2, data.community2$community=="new"))
new.graze<-multipatt(new.planting[,-c(1:6)], cluster=new.planting$graze, func = "IndVal", duleg=TRUE)
summary(new.graze)
#no significant indicators

old.planting<-droplevels(subset(data.community2, data.community2$community=="old"))
old.graze<-multipatt(old.planting[,-c(1:6)], cluster=old.planting$graze, func = "IndVal", duleg=TRUE)
summary(old.graze)
#LESCAP and EUPALT are indicators for ungrazed
#Group UG  #sps.  2 
#stat p.value   
#Lespedeza.capitata    0.745    0.01 **
#  Eupatorium.altissimum 0.527    0.05 * 

remnant<-droplevels(subset(data.community2, data.community2$community=="remnant"))
remnant.graze<-multipatt(remnant[,-c(1:6)], cluster=remnant$graze, func = "IndVal", duleg=TRUE)
summary(remnant.graze)
#SORNUT and BAPBRA indicators for ungrazed
#Group UG  #sps.  2 
#stat p.value   
#Sorghastrum.nutans  0.839   0.005 **
#  Baptisia.leucophaea 0.660   0.020 * 

#overall community indicators
community<-multipatt(data.community2[,-c(1:6)], cluster=data.community2$community, func = "IndVal", duleg=TRUE)
summary(community)
#LOTS!
#Supplemental table? See Excel table, which is easier to find if we wanted to include.

#####################################################
#pull out NMDS scores and run regressions with envfit
#data_veg is datasciname
#have to run ordination code for full data set before doing this
#will need to bake in these vectors with teh ordination code, since we're working with 3 dimensions

#graze
envectors_graze<-envfit(exclosure.NMDS[,1:3] ~ exclosure.NMDS$graze, na.rm=TRUE)

community.graph.cent2<-ggplot.NMDS.3pointCent(mds.ab3, exclosure.NMDS$community, colors.community)


data.scores <- as.data.frame(scores(mds.ab3, choices = c(1,2,3), display = c("sites")))
data.scores$site_ID <- rownames(data.community2$siteID)

data.scores$grazed <-(data_veg$graze)
data.scores$site_ID <-(data_veg$siteID)
data.scores$plot_year <-(data_veg$plot_year)
data.scores$community <-(data_veg$community)

head(data_veg)
head(data.scores)
#########

MDS1<-data.frame((data.scores))$NMDS1

MDS2<-data.frame(data.scores$NMDS2)

Transect<-data_veg$plot_year

Sample_Year<-as.numeric(paste(data_veg$year))

Community<-data_veg$community

Graze<-data_veg$graze

burn <- data_veg$years_since_burn

full.NMDS<-data.frame(MDS1, MDS2, Transect, Sample_Year,Community,Graze, burn)

head(savanna.scores)

savanna.scores<-droplevels(subset(full.NMDS, full.NMDS$community=="savanna"))

envectors08<-envfit(full.NMDS[,1:2] ~ full.NMDS$burn, na.rm=TRUE)

head(full.NMDS)

head(envectors08)

#Sample_Year P=0.001



#extract vectors, center them at origin

vectors.08<-data.frame(envectors08$vectors[1:2])

library(dplyr)

#add vector to graph

savanna.graph3<-ggplot.NMDS(mds.pa2, data_veg$plot_year, colors.13)+
  annotate("text", x=1.1, y=-1.0, label="stress=0.133", size=6)+labs(title = "Savanna")+
  theme(plot.title = element_text(hjust=0.5, size=20, face="bold"), 
        panel.grid= element_blank())+
  
  geom_segment(data=envectors08, aes(x=0,xend=arrows.MDS1, y=0,yend=arrows.MDS2),
               arrow=arrow(length = unit(0.35, "cm")),colour="darkgrey",
               size=1,inherit_aes=FALSE)

savanna.graph3

####################################################################################



