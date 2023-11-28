#Jenn Chakravorty
#Elizabeth Bach
#Nachusa Grassland bison exclosure - data analysis
#15 Dec 2022

##clear old/other data
rm(list=ls())

### First, let us install and load the "lme4" package

library(lme4)

#Load tidyverse
library(tidyverse)

#read in dataset, use "NachusaBisonExclosureDivAbund.csv"
Nachusa.data2<-read.csv(file.choose(),header = TRUE,na.strings = "NA")

#to run with mixed model, will want bison_year to be categorical
Nachusa.data2$bison_year<-as.factor(Nachusa.data2$bison_year)

#Shannon's diversity
#Run mixed models with bison year rather than calendar year
#Full model
modelF <- lme(shannon ~ graze*community*bison_year, random = ~1|years_since_burn, data=Nachusa.data2)
summary(modelF)
anova(modelF)
#no interactions

model1 <- lme(shannon ~ graze+community+bison_year, random = ~1|years_since_burn, data=Nachusa.data2)
summary(model1)
anova(model1)

#Using lme, AIC value with 3-way interaction is 97, with just main factors, AIC=53
#community is significant P<0.0001
#bison_year is significant P=0.0007

#pairwise comparisons for differences
library(emmeans)
pairs(emmeans(model1, ~graze|bison_year))
#no diff in grazing
pairs(emmeans(model1, ~graze|community))
#graze/ungraze no diff in any community
#Shannon's diversity different within community and across years since bison, no effect of grazing
pairs(emmeans(model1, pairwise~community))
#old and remnant only ones NOT different
pairs(emmeans(model1, pairwise~bison_year))
#Year 5 different from 0 and 3, years 0 and 3 NOT different

#graphs
#summarize by transect
#transform some characters to factors
Nachusa.data2$plot_year<-as.factor(Nachusa.data2$plot_year)
Nachusa.data2$community<-as.factor(Nachusa.data2$community)
Nachusa.data2$graze<-as.factor(Nachusa.data2$graze)

shannon.sum<- Nachusa.data2 %>%
  group_by(community, bison_year) %>%
  summarise(mean=mean(shannon), n=length(shannon), SE=sd(shannon)/sqrt(n))

shannons<-ggplot(shannon.sum, aes(x=bison_year, y=mean, colour=community))+geom_pointrange(aes(ymin=(mean-SE), ymax=(mean+SE)),size=1, position = position_dodge(width=0.25))+ylim(0,4)+ylab("Shannon's Diversity")
shannons

#community main effect
shannon.com<- Nachusa.data2 %>%
  group_by(community) %>%
  summarise(mean=mean(shannon), n=length(shannon), SE=sd(shannon)/sqrt(n))

shannons1.5<-ggplot(shannon.com, aes(x=community, y=mean))+geom_pointrange(aes(ymin=(mean-SE), ymax=(mean+SE)),size=1, position = position_dodge(width=0.25))+ylim(0,3.9)+ylab("Shannon's Diversity")+
  annotate("text", x=1, y=3.65, label="A", size=5)+annotate("text", x=2, y=2.95, label="C", size=5)+annotate("text", x=3.0, y=3.0, label="C", size=5)+annotate("text", x=4, y=3.3, label="B", size=5)+annotate("text", x=0.6, y=3.9, label="A)", size=8)+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, plot.title = element_text(colour="black", size=18, face="bold", hjust=0.5), axis.text.y=element_text(colour="black", size=16), axis.text.x=element_text(colour="black", size=16), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title=element_text(colour="black", size=18, face="bold"), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))
shannons1.5

#Bison_year main effect
shannon.yr<- Nachusa.data2 %>%
  group_by(bison_year) %>%
  summarise(mean=mean(shannon), n=length(shannon), SE=sd(shannon)/sqrt(n))

shannons3<-ggplot(shannon.yr, aes(x=bison_year, y=mean))+geom_pointrange(aes(ymin=(mean-SE), ymax=(mean+SE)),size=1, position = position_dodge(width=0.25))+ylim(0,3.4)+ylab("Shannon's Diversity")+
  annotate("text", x=0, y=3.05, label="B", size=5)+annotate("text", x=3, y=3.1, label="B", size=5)+annotate("text", x=5, y=3.25, label="A", size=5)+xlab("Years since bison")+annotate("text", x=0, y=3.38, label="A)", size=8)+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, plot.title = element_text(colour="black", size=18, face="bold", hjust=0.5), axis.text.y=element_text(colour="black", size=16), axis.text.x=element_text(colour="black", size=16), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title=element_text(colour="black", size=18, face="bold"), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))
shannons3


#nonnative:native ratio
#graphs
nonnative_ratio<-ggplot(Nachusa.data2, aes(x=bison_year, y=ratio_non_native, colour=community))+geom_point()
nonnative_ratio
#outlier in old rest, year 3?
range(Nachusa.data2$ratio_non_native)
outlier<-droplevels(subset(Nachusa.data2, Nachusa.data2$ratio_non_native>3))
outlier
#N8 2017, grazed - this plot had high levels of red clover cover (70-90% in several plots)
#exclude for now
Nachusa.ratio<-droplevels(subset(Nachusa.data2, Nachusa.data2$ratio_non_native<3))
#transform for inverse ratio
Nachusa.ratio$inv_ratio<-(1/(Nachusa.ratio$ratio_non_native+1))

model.ratio <- lme(inv_ratio ~ graze*bison_year*community, random = ~1|years_since_burn, data=Nachusa.ratio)
summary(model.ratio)
anova(model.ratio)

#Full model: significant graze*bison year interaction, AIC -54
#community significant
model.ratio2 <- lme(inv_ratio ~ graze*bison_year+community, random = ~1|years_since_burn, data=Nachusa.ratio)
summary(model.ratio2)
anova(model.ratio2)

emmeans(model.ratio2, pairwise~community)
#new different from old and remnant

emmeans(model.ratio2, pairwise~graze|bison_year)
#graze/ungraze different in year 5 only

#graph
#main effect of community
Non_native_ratio<-tibble(
Nachusa.ratio %>%
  group_by(community) %>%
  summarise(mean=mean(ratio_non_native), n=length(ratio_non_native), se=sd(ratio_non_native)/sqrt(n)))

nonnative_ratio1<-ggplot(Non_native_ratio)+geom_pointrange(aes(x=community, y=mean, ymin=(mean-se), ymax=(mean+se)), size=1)+ylim(0,0.8)+ylab("Non-native:Native ratio")+
  annotate("text",x=1, y=0.7, label="A", size=5)+annotate("text",x=2, y=0.47, label="B", size=5)+annotate("text",x=3, y=0.4, label="B", size=5)+annotate("text",x=4, y=0.5, label="AB", size=5)+annotate("text",x=0.65, y=0.79, label="B)", size=8)+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, plot.title = element_text(colour="black", size=18, face="bold", hjust=0.5), axis.text.y=element_text(colour="black", size=16), axis.text.x=element_text(colour="black", size=16), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title=element_text(colour="black", size=18, face="bold"), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"))
nonnative_ratio1

#graze*bison_year interaction
Non_native_ratio2<-tibble(
  Nachusa.ratio %>%
    group_by(bison_year, graze) %>%
    summarise(mean=mean(ratio_non_native), n=length(ratio_non_native), se=sd(ratio_non_native)/sqrt(n)))

nonnative_ratio3<-ggplot(Non_native_ratio2)+geom_pointrange(aes(x=bison_year, y=mean, ymin=(mean-se), ymax=(mean+se), colour=graze), position=position_dodge(width=0.5), size=1)+ylim(0,0.7)+ylab("Non-native:Native ratio")+xlab("Years since bison")+xlim(-0.2,5.2)+
  annotate("text",x=4.9, y=0.68, label="A", size=5)+annotate("text",x=5.12, y=0.4, label="B", size=5)+annotate("text",x=-0.1, y=0.5, label="ND", size=5)+annotate("text",x=3.1, y=0.5, label="ND", size=5)+scale_color_manual(labels = c("Grazed", "Ungrazed"), values=c("black","gray"))+annotate("text",x=0, y=0.7, label="C)", size=8)+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, plot.title = element_text(colour="black", size=18, face="bold", hjust=0.5), axis.text.y=element_text(colour="black", size=16), axis.text.x=element_text(colour="black", size=16), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title=element_text(colour="black", size=18, face="bold"), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"), legend.title = element_blank(), legend.text = element_text(colour="black", size=14, hjust=0.5, face="bold"), legend.text.align = 0, legend.key = element_blank(), legend.position = c(0.8,0.15))
nonnative_ratio3

#grass:forb ratio
#graphs
GrassForb_ratio<-ggplot(Nachusa.data2, aes(x=bison_year, y=grass_forb, colour=community))+geom_point()
GrassForb_ratio
#fairly uniform distribution, no outliers
#natural log transform

model.GF <- lme(log(grass_forb) ~ graze*community+bison_year, random = ~1|years_since_burn, data=Nachusa.data2)
summary(model.GF)
anova(model.GF)

#graze*community interaction

#Full model AIc=393, same for main effects only?
#community P<0.0001
#bison_year P=0.0101

emmeans(model.GF, pairwise~community)
#New diff from old and remnant, remnant and savanna different

emmeans(model.GF, pairwise~bison_year)
#year 0 diff from 3, year 3 diff from 5

emmeans(model.GF, pairwise~graze|community)
#in savanna only: graze/ungraze difference P=0.0001

savanna<-as.tibble(droplevels(subset(Nachusa.data2, Nachusa.data2$community=="savanna")))
GrassForb_ratio_sav<-ggplot(savanna, aes(x=bison_year, y=grass_forb, colour=graze))+geom_point()
GrassForb_ratio_sav



#graph
GrassForb_ratio_sum<-tibble(
  Nachusa.data2 %>%
    group_by(community, bison_year) %>%
    summarise(mean=mean(grass_forb), n=length(grass_forb), se=sd(grass_forb)/sqrt(n)))

GF_ratio<-ggplot(GrassForb_ratio_sum)+geom_pointrange(aes(x=bison_year, y=mean, ymin=(mean-se), ymax=(mean+se),colour=community),size=1, position=position_dodge(width=0.5))+ylim(0,3.6)+ylab("Grass:Forb ratio")
GF_ratio

#clean up for a supplemental
Sav_grouped<-
  savanna %>%
    dplyr::group_by(bison_year, graze)%>%
    dplyr::summarize(mean=mean(grass_forb), n=length(grass_forb), se=sd(grass_forb)/sqrt(n))

graze.colors<-c("black","gray")                         
Sav_GF_pub<-ggplot(Sav_grouped)+geom_pointrange(aes(x=bison_year, y=mean, ymin=(mean-se), ymax=(mean+se), colour=graze),size=1, lwd=1.2)+ylim(0,3.6)+ylab("Grass:Forb ratio")+scale_color_manual(values=graze.colors, labels=c("grazed","ungrazed"))+
    theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, plot.title = element_text(colour="black", size=18, face="bold", hjust=0.5), axis.text.y=element_text(colour="black", size=16), axis.text.x=element_text(colour="black", size=16), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title=element_text(colour="black", size=18, face="bold"), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"), legend.title = element_blank(), legend.box.background=element_blank(),legend.text = element_text(colour="black", size=14, hjust=0.5, face="bold"), legend.text.align = 0, legend.key = element_blank(), legend.position = c(0.8,0.1))
Sav_GF_pub

#strong decrease in old and remnant sites seem to drive the main effect of year, interesting to see that ratio fluctuate that much (bounces right back in year 5)
#take a look at graze/ungraze, is that driving this?
GrassForb_graze<-Nachusa.data2 %>%
    dplyr::group_by(graze, community) %>%
    dplyr::summarize(mean=mean(grass_forb), n=length(grass_forb), se=sd(grass_forb)/sqrt(n))

GF_ratio_graze2<-ggplot(GrassForb_graze)+geom_pointrange(aes(x=community, y=mean, ymin=(mean-se), ymax=(mean+se),colour=community),size=1, position=position_dodge(width=0.5))+ylim(0,3.6)+ylab("Grass:Forb ratio")
GF_ratio_graze2
#the only graze/ungraze difference is in savanna, the year 3 drop is consistent in old and remnant
#savanna differences consistent across all three sampling points (e.g. the difference pre-dates the bison)

#community level effect for publication
GrassForb_com<-tibble(
  Nachusa.data2 %>%
    group_by(community) %>%
    summarise(mean=mean(grass_forb), n=length(grass_forb), se=sd(grass_forb)/sqrt(n)))

GF_ratio_comm2<-ggplot(GrassForb_com)+geom_pointrange(aes(x=community, y=mean, ymin=(mean-se), ymax=(mean+se)),size=1)+ylim(0,3.6)+ylab("Grass:Forb ratio")+
  annotate("text", x=1, y=1.1, label="B", size=5)+annotate("text", x=2, y=2.7, label="A", size=5)+annotate("text", x=3, y=3.2, label="A", size=5)+annotate("text", x=4, y=2.0, label="B", size=5)+annotate("text", x=0.65, y=3.5, label="C)", size=8)+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, plot.title = element_text(colour="black", size=18, face="bold", hjust=0.5), axis.text.y=element_text(colour="black", size=16), axis.text.x=element_text(colour="black", size=16), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title=element_text(colour="black", size=18, face="bold"), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"), legend.title = element_blank(), legend.text = element_text(colour="black", size=14, hjust=0.5, face="bold"), legend.text.align = 0, legend.key = element_blank(), legend.position = c(0.8,0.15))
GF_ratio_comm2

#year main effect for publication
GrassForb_yr<-tibble(
  Nachusa.data2 %>%
    group_by(bison_year) %>%
    summarise(mean=mean(grass_forb), n=length(grass_forb), se=sd(grass_forb)/sqrt(n)))

GF_ratio_yr2<-ggplot(GrassForb_yr)+geom_pointrange(aes(x=bison_year, y=mean, ymin=(mean-se), ymax=(mean+se)),size=1)+ylim(0,3.5)+ylab("Grass:Forb ratio")+
  annotate("text", x=0, y=2.75, label="A", size=5)+annotate("text", x=3, y=2.0, label="B", size=5)+annotate("text", x=5, y=2.75, label="A", size=5)+xlab("Years since bison")+annotate("text",x=0.1, y=3.4, label="B)", size=8)+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, plot.title = element_text(colour="black", size=18, face="bold", hjust=0.5), axis.text.y=element_text(colour="black", size=16), axis.text.x=element_text(colour="black", size=16), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title=element_text(colour="black", size=18, face="bold"), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"), legend.title = element_blank(), legend.text = element_text(colour="black", size=14, hjust=0.5, face="bold"), legend.text.align = 0, legend.key = element_blank())
GF_ratio_yr2

#robel data only collected in 2019/2020, so run different model on different dataset, no bison_year comparisons!
#use RobelUpdate.csv
robel.data<-read.csv(file.choose(),header = TRUE,na.strings = "NA")

#drop wetlands
robel.data$community_type<-as.factor(robel.data$community_type)
robel.data2<-droplevels(subset(robel.data, robel.data$community_type!="Wetland"))


VOR.data<-ggplot(robel.data2, aes(x=community_type, y=robelavg, colour=graze))+geom_point()
VOR.data

#fairly uniform distribution, none of these seem to pop as an outlier.

model.VOR <- lme(log(robelavg) ~ graze+community_type, random = ~1|years_since_burn, data=robel.data2)
summary(model.VOR)
anova(model.VOR)

#main model AIC 88, no interactions in full, stick with main

#community type the only statistically significant factor: P=0.017
emmeans(model.VOR, pairwise~community_type)
#new and remnant different

#graze is on the cusp, P=0.068
emmeans(model.VOR, pairwise~graze)

robel.comm<-tibble(
  robel.data2 %>%
    group_by(community_type) %>%
    summarise(mean=mean(robelavg), n=length(robelavg), se=sd(robelavg)/sqrt(n)))

robel.graph1<-ggplot(robel.comm)+geom_pointrange(aes(x=community_type, y=mean, ymin=(mean-se), ymax=(mean+se)),size=1)+ylim(0,16)+ylab("Visual Obstruction Ratio")+xlab("community")+
  annotate("text",x=1, y=13.25, label="A", size=5)+annotate("text",x=2, y=8, label="AB", size=5)+annotate("text",x=3, y=5.25, label="B", size=5)+annotate("text",x=4, y=15.5, label="AB", size=5)+annotate("text",x=0.65, y=15.8, label="D)", size=8)+
  theme(axis.line=element_line(colour="black", size=1.5), axis.ticks=element_line(colour="black", size=1), aspect.ratio=1, plot.title = element_text(colour="black", size=18, face="bold", hjust=0.5), axis.text.y=element_text(colour="black", size=16), axis.text.x=element_text(colour="black", size=16), axis.ticks.length=unit(0.2,"cm"), panel.grid.major=element_blank(), panel.grid.major.y=element_line(colour="gray", size=1), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.title=element_text(colour="black", size=18, face="bold"), strip.text = element_text(colour="black", size=14, hjust=0.5, face="bold"), legend.title = element_blank(), legend.text = element_text(colour="black", size=14, hjust=0.5, face="bold"), legend.text.align = 0, legend.key = element_blank())
robel.graph1


#graze
robel.gr<-tibble(
  robel.data2 %>%
    group_by(graze) %>%
    summarise(mean=mean(robelavg), n=length(robelavg), se=sd(robelavg)/sqrt(n)))

robel.graph2<-ggplot(robel.gr)+geom_pointrange(aes(x=graze, y=mean, ymin=(mean-se), ymax=(mean+se)),size=1)+ylim(0,16)+ylab("Visual Obstruction Ratio")
robel.graph2

#just to see the grazing impact, not a statisitcally significant effect
robel.comm.graze<-tibble(
  robel.data2 %>%
    group_by(community_type, graze) %>%
    summarise(mean=mean(robelavg), n=length(robelavg), se=sd(robelavg)/sqrt(n)))

robel.graph2<-ggplot(robel.comm.graze)+geom_pointrange(aes(x=community_type, y=mean, ymin=(mean-se), ymax=(mean+se), colour=graze),size=1)+ylim(0,16)+ylab("Visual Obstruction Ratio")
robel.graph2

#Build merged figures for publication
#Figure A: Community main effects for Shannon's, non-native:native, grass:forb, VOR
library(gridExtra)
grid.arrange(shannons1.5, nonnative_ratio1, GF_ratio_comm2, robel.graph1, newpage = TRUE, ncol=2)

#Fig. B: Year effect, including year*bison grazing interaction
grid.arrange(shannons3, GF_ratio_yr2, nonnative_ratio3, newpage = TRUE, ncol=3)

###########################
