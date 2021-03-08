## Written by Jason M. Matlock and Natalie S. Vande Pol
## January-March 2020

#Set Working Directory to Source Folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
# install.packages('patchwork')


###############################################################################
############################    SOIL EXPERIMENT    ############################
###############################################################################

library(ggpubr)
library(ggsignif)

plot_order <-c("NVP64cu", "NVP64wt", "NVP80cu", "NVP80wt","NoMillet", "Uninoculated")
soil_colors<-c("#33ccff", "#0099cc", "#00cc00", "#009900","#cccccc",  "#999999")

###################################
#####    AERIAL DRY WEIGHT    #####
###################################

soil.raw<-read.table("Endobacteria_Panel.txt", 
                     sep="\t", header=TRUE, 
                     col.names =c("factor","rep","fresh","dry") )


## Fig## - Aerial Dry Biomass (Soil) ##
ggboxplot(data=soil.raw, x="factor", y="dry", ylab = "Dry Weight (g)",xlab="",
          legend="none", order=plot_order, fill="factor", palette=soil_colors,
          add="dotplot")+
  stat_compare_means(comparisons = list(c("NVP64cu", "NVP64wt"),  
                                        c("NVP80cu", "NVP80wt"),
                                        c("NoMillet", "Uninoculated")),
                     label = "p.format",
                     p.adjust.method="holm",
                     method="t.test", label.y=1.45)+
  stat_compare_means(comparisons = list(c("NVP64cu", "Uninoculated"),
                                        c("NVP64wt", "Uninoculated"),
                                        c("NVP80cu", "Uninoculated"),
                                        c("NVP80wt", "Uninoculated")),
                     label = "p.format",
                     p.adjust.method="holm",
                     method="t.test",
                     method.args = list(alternative = "greater"),
                     label.y=c(1.6,1.75,1.9,2.05))



#################################
#####    SEED PRODUCTION    #####
#################################

seed.raw <- read.table("Avg_Seed_Mass.csv", sep=",", header=T)

#### DATA ANALYSIS & VISUALIZATION ####
seed.raw$Avg_Mass_ug <- 1000*seed.raw$Subset_Mass_mg/seed.raw$Count
seed.raw$Total_Seed_Num <- 1000*seed.raw$Total_Mass_g/seed.raw$Avg_Mass_ug


## Fig## - Seed Production (Soil) (panel a) ##
ggboxplot(data=seed.raw, x="Treatment", y="Total_Mass_g",
                         fill="Treatment", ylab = "Total Seed Mass (mg)", 
                         ylim=c(0,500), xlab="", order=plot_order, 
                         palette=soil_colors, legend="none",
          add="dotplot")+
  stat_compare_means(comparisons = list(c("NVP64cu", "NVP64wt"),  
                                        c("NVP80cu", "NVP80wt"),
                                        c("NoMillet", "Uninoculated")),
                     label = "p.format",
                     p.adjust.method="holm",
                     label.y=350)+
  stat_compare_means(comparisons = list(c("NVP64cu", "Uninoculated"),
                                        c("NVP64wt", "Uninoculated"),
                                        c("NVP80cu", "Uninoculated"),
                                        c("NVP80wt", "Uninoculated")),
                     label = "p.format",
                     p.adjust.method="holm",
                     method.args = list(alternative = "greater"),
                     label.y=c(380,410,440,470))


## Fig## - Seed Production (Soil) (panel b) ##
ggboxplot(data=seed.raw, x="Treatment", y="Avg_Mass_ug",
          color="black", fill="Treatment", 
          ylab = expression(paste("Avg Seed Mass (", mu, "g)")),
          order=plot_order, palette=soil_colors, ylim=c(0,30), legend="none",
          add="dotplot")+
  stat_compare_means(comparisons = list(c("NVP64cu", "NVP64wt"),  
                                        c("NVP80cu", "NVP80wt"),
                                        c("NoMillet", "Uninoculated")),
                     label = "p.format",
                     p.adjust.method="holm",
                     label.y=22)+
  stat_compare_means(comparisons = list(c("NVP64cu", "Uninoculated"),
                                        c("NVP64wt", "Uninoculated"),
                                        c("NVP80cu", "Uninoculated"),
                                        c("NVP80wt", "Uninoculated")),
                     label = "p.format",
                     p.adjust.method="holm",
                     method.args = list(alternative = "less"),
                     label.y=c(24,26,28,30))+ 
  xlab("")


## Fig## - Seed Production (Soil) (panel c) ##
ggboxplot(data=seed.raw, x="Treatment", y="Total_Seed_Num",
          color="black", fill="Treatment", 
          ylab = "Total Seed Number", ylim=c(0,45000), 
          xlab="", order=plot_order, palette=soil_colors, legend="none",
          add="dotplot")+
  stat_compare_means(comparisons = list(c("NVP64cu", "NVP64wt"),  
                                        c("NVP80cu", "NVP80wt"),
                                        c("NoMillet", "Uninoculated")),
                     label = "p.format",
                     p.adjust.method="holm",
                     label.y=31000)+
  stat_compare_means(comparisons = list(c("NVP64cu", "Uninoculated"),
                                        c("NVP64wt", "Uninoculated"),
                                        c("NVP80cu", "Uninoculated"),
                                        c("NVP80wt", "Uninoculated")),
                     label = "p.format",
                     p.adjust.method="holm",
                     method.args = list(alternative = "greater"),
                     label.y=c(34000,37000,40000,43000))


#summary(aov(Avg_Mass_ug~Trt, data = seed.raw))
detach("package:ggpubr",unload=TRUE)
detach("package:ggsignif",unload=TRUE)


## SuppFig## - Seed Area Density Plot ##
all.reports.raw<-read.csv("all_reports_raw.csv", stringsAsFactors = F)

all.reports.batch1<- subset(all.reports.raw, grepl("Batch1", Trt))
all.reports.batch1$Rep<-rep("",length(all.reports.batch1$Trt))
all.reports.batch1$Strain<-rep("",length(all.reports.batch1$Trt))
for (j in 1:length(all.reports.batch1$Trt)){
  temp<-strsplit(all.reports.batch1$Trt[j],"_")
  all.reports.batch1$Strain[j]<- temp[[1]][2]
  all.reports.batch1$Rep[j]<-temp[[1]][3]
  
  if (all.reports.batch1$Area[j]>300){
    all.reports.batch1$Area[j]<-NA
  }
  
}

all.reports.batch1$Strain <- as.factor(all.reports.batch1$Strain)
levels(all.reports.batch1$Strain)<-c("NVP64cu","NVP64wt","NVP80cu","NVP80wt",
                                     "NoMillet","Uninoculated")
my_colors_den <-c("#33ccff", "#0099cc", "#00cc00", "#009900","#999999",  "#666666")


library(plyr)
mu <- ddply(all.reports.batch1, "Strain", summarise, grp.mean=mean(Area, na.rm=T))

ssd<-ggplot(all.reports.batch1, aes(x=Area, color=Strain))+
  ylab("Proportion of Seeds Scanned")+
  xlab("Area of Seed (Pixels)")+
  geom_density(size=1)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Strain), linetype="dashed", size=1)+
  scale_color_manual(values=my_colors_den)
ssd



###############################################################################
##########################    BOLTING EXPERIMENT    ###########################
###############################################################################
agar_colors <-c("#999999", "#33ccff", "#0099cc", "#00cc00", "#009900")
data.bolt <- read.csv("Bolting_Panel_Data.csv",header=T,na="NA",
                      col.names = c("Plate", "Treatment", "Age", "Light_Lev"),
                      colClasses= c("Plate"     = "factor",
                                    "Treatment" = "factor",
                                    "Age"       = "numeric",
                                    "Light_Lev" = "numeric"))

ggplot(data.bolt)+
  geom_boxplot(aes(x=Treatment,y=Age))+
  ylim(0,26)
#data are not normally distributed, use kruskal-wallis test, not anova

kruskal.test(Age~Treatment, data=data.bolt)
mean(data.bolt$Age, na.rm=T)
sd(data.bolt$Age, na.rm=T)/sqrt(length(data.bolt$Age))



###############################################################################
############################    AGAR EXPERIMENT    ############################
###############################################################################

#Load Packages
library(lme4)
library(car)

#Load Data from CSV
data.raw<-read.csv("Combined_Data_Stats.csv", 
                   header=TRUE, 
                   col.names= c("Round", "Treatment", "SampleID", "SubID", 
                                "Root_Len","Dry_Wt", "Medium", "Lineage", 
                                "Light_Lev"),
                   colClasses = c("Round"      ="factor",
                                  "Treatment"  ="factor",
                                  "SampleID"   ="factor",
                                  "SubID"      ="factor", 
                                  "Root_Len"   ="numeric",
                                  "Dry_Wt"     ="numeric",
                                  "Medium"     ="factor", 
                                  "Lineage"    ="factor", 
                                  "Light_Lev"  ="numeric"),
                   na="NA")

################################
#####    CURED LINEAGES    #####
################################

#isolate the correct data subset
data.lineage<- data.raw %>% filter(Treatment %in% c("NVP64cu","NVP80cu") ) %>%
  droplevels()

#visualize data  
ggplot(data.lineage) + geom_boxplot(aes (x=Lineage, y=Dry_Wt)) +
  facet_grid(Treatment~Round)
# Data is normally distributed across all treatments

# Extract CuredPanel where there is full replication 
data.lineage_rC<- data.raw %>% filter(Treatment %in% c("NVP64cu","NVP80cu") ) %>%
  filter(Round=="CuredPanel") %>% droplevels()

## Plot average plate dry weight by cured lineage (Supplementary Figure XX) ##
ggplot(data.lineage_rC) + 
  geom_boxplot(aes (x=Lineage, y=Dry_Wt)) +
  facet_grid(~Treatment)+
  ylab("Average Dry Weight (mg)")+
  xlab("Cured Lines")

#Perform One-Way ANOVA for each strain, averaging subsample response values

anova.lineage64<-aov(avg_wt~Lineage, 
                     data=data.lineage_rC %>% filter(Treatment=="NVP64cu")%>%
                       group_by(Lineage, Round, SampleID)%>%
                       summarise(avg_wt=mean(Dry_Wt, na.rm=T))%>%droplevels() )
plot(anova.lineage64, 1) 
summary(anova.lineage64) 

anova.lineage80<-aov(avg_wt~Lineage,
                     data=data.lineage_rC %>% filter(Treatment=="NVP80cu")%>%
                       group_by(Lineage, Round, SampleID)%>%
                       summarise(avg_wt=mean(Dry_Wt, na.rm=T))%>%droplevels() )
plot(anova.lineage80, 1)
summary(anova.lineage80) 

# No Effect of Lineage. Can Be Collapsed in All Future Analyses



##############################
#####    LIGHT LEVEL    ######
##############################

# No Light Level Data in CuredPanel, MediaPanel is only complete design

#Filtered out CuredPanel
data.light<-data.raw %>% filter(Round !="CuredPanel") %>%
                         group_by(Round,Treatment,Medium,Light_Lev, SampleID) %>%
                         summarise (avg_wt=mean(Dry_Wt, na.rm=T))%>%
                         droplevels()

## Plot average plate dry weight by light level (Supplementary Figure XX) ##
ggplot(data.light, aes (x=Light_Lev, y=avg_wt))   + 
                        geom_point() +
                        facet_grid(Treatment~Medium)+
                        geom_smooth(method='lm')+
  xlab("Light Level (umol)")+
  ylab("Average Dry Weight (mg)")

lm.light <-lm(avg_wt~ Light_Lev + Treatment
              + Light_Lev*Treatment + Light_Lev*Medium + Medium*Treatment 
              + Light_Lev*Medium*Treatment,
              data=data.light)

lm.light_nomed <-lm(avg_wt~ Light_Lev + Treatment + Light_Lev*Treatment,
                    data=data.light)

library(emmeans)
anova(lm.light)
anova(lm.light_nomed)

emtrends(lm.light, ~Treatment|Medium, var="Light_Lev")  

#sliced by Medium
pairs(emtrends(lm.light, ~Treatment|Medium, var="Light_Lev"))
CLD(  emtrends(lm.light, ~Treatment|Medium, var="Light_Lev"),
               alpha=0.05, Letters=letters,adjust="tukey")

#true marginal estimates
pairs(emtrends(lm.light, ~Treatment*Medium, var="Light_Lev"))
CLD  (emtrends(lm.light, ~Treatment*Medium, var="Light_Lev"),
      alpha=0.05, Letters=letters,adjust="tukey")

detach("package:emmeans",unload=TRUE)
#Findings: No significant difference between slopes, 95% CI estimates include 0, 
#suggesting no significant effect of light level



########################
#####    MEDIUM    #####
########################

#Filtered out Media Panel
data.media<-data.raw %>% filter(Round =="MediaPanel") %>%
  group_by(Round, Treatment, Medium, SampleID) %>%
  summarise (avg_wt=mean(Dry_Wt, na.rm=T))%>%
  droplevels()

## Plot Avg Plate Dry_Wt by Medium & Treatment ##
ggplot(data.media)+
  geom_boxplot(aes (x=Medium, y=avg_wt))+
  facet_grid(~Treatment)+
  ylab("Average Dry Weight (mg)")+
  xlab("Medium")

lm.media <-lm(avg_wt~ Medium + Treatment + Medium*Treatment,data=data.media)


library(emmeans)
anova(lm.media)

marginal.media<-emmeans(lm.media, pairwise~Medium|Treatment,adjust="tukey") 
marginal.media
CLD(marginal.media$emmeans,alpha=0.05, Letters=letters,adjust="tukey")
detach("package:emmeans",unload=TRUE)

# No Effect of Medium. Can Be Collapsed in All Future Analyses



####################################
#####    ROOT LENGTH ANCOVA    #####
####################################

data.root<- data.raw %>%
  dplyr::select(SampleID, Round, Treatment, SubID, Root_Len, Dry_Wt) %>%
  droplevels()
data.root_rM<- data.root %>% filter(Round=="MediaPanel") %>% droplevels()
data.root_rC<- data.root %>% filter(Round=="CuredPanel") %>% droplevels()

# Plot Root Lenght by Round & Treatment  
ggplot(data.root) + geom_boxplot(aes (x=Treatment, y=Root_Len)) +
         facet_grid(~Round)

## Plot Dry_Wt by Root-Len ##
ggplot(data.root,aes (x=Root_Len, y=Dry_Wt))   +
  geom_point() +
  geom_smooth(method='lm', formula=y~x)+
  facet_grid(Treatment~Round)+
  xlab("Root Length (mm)")+
  ylab("Average Dry Weight (mg)")

# Preliminary Linear Models
lm.root <-lm(Dry_Wt~ Root_Len + Treatment + Round +
                + Root_Len*Treatment 
                + Root_Len*Round 
                + Round*Treatment 
                + Root_Len*Round*Treatment, 
                data=data.root)

lm.root_rM <-lm(Dry_Wt~ Root_Len + Treatment + Root_Len*Treatment, 
                data=data.root_rM)
lm.root_rC <-lm(Dry_Wt~ Root_Len + Treatment + Root_Len*Treatment,
                data=data.root_rC)


library(emmeans)
anova(lm.root_rM)
emtrends(lm.root_rM, ~Treatment, var="Root_Len")  

pairs(emtrends(lm.root_rM, ~Treatment, var="Root_Len"))
CLD(  emtrends(lm.root_rM, ~Treatment, var="Root_Len"),
      alpha=0.05, Letters=letters,adjust="tukey")

anova(lm.root_rC)
emtrends(lm.root_rC, ~Treatment, var="Root_Len")  

pairs(emtrends(lm.root_rC, ~Treatment, var="Root_Len"))
CLD(  emtrends(lm.root_rC, ~Treatment, var="Root_Len"),
      alpha=0.05, Letters=letters,adjust="tukey")


anova(lm.root)
emtrends(lm.root, ~Treatment|Round, var="Root_Len")  

#sliced by round
pairs(emtrends(lm.root, ~Treatment|Round, var="Root_Len"))
CLD(  emtrends(lm.root, ~Treatment|Round, var="Root_Len"),
      alpha=0.05, Letters=letters,adjust="tukey")

#true marginal estimates
pairs(emtrends(lm.root, ~Treatment*Round, var="Root_Len"))
CLD  (emtrends(lm.root, ~Treatment*Round, var="Root_Len"),
      alpha=0.05, Letters=letters,adjust="tukey")

detach("package:emmeans",unload=TRUE)


lm.root_len_dist <- lm(Root_Len ~ Round + Treatment + Round*Treatment,
                       data=data.root)
anova(lm.root_len_dist)


###################################
#####    FULL LINEAR MODEL    #####
###################################

# Round & SampleID are blocking factors
# Treatment and Root_Len are fixed effects

ggplot(data.root)+
  geom_boxplot(aes(x=Treatment, y=Dry_Wt))
# Data appear normally distributed with equal variance

library(lmerTest)

nested.lm <- lmer(Dry_Wt ~ Treatment + Root_Len + (1|Round) + (1|SampleID),
                  data=data.root)
summary(nested.lm)


library(emmeans)
lm.emm <- emmeans(nested.lm, ~Treatment)
pairs(emmeans(nested.lm, ~Treatment))
CLD(lm.emm)


#library(lattice)
# Plot Linear Model Coefficients
#library(coefplot)
#coefplot(nested.lm)

# Preliminary EMM Plot by treatment
plot(lm.emm, horizontal=F, comparisons = TRUE)
detach("package:emmeans",unload=TRUE)


## Figure ## - Estimated Marginal Mean Dry Aerial Biomass by Treatment ##

trt <- c("Ctrl", "NVP64cu", "NVP64wt", "NVP80cu", "NVP80wt")
emm <- c(1.61, 2.17, 2.21, 2.26, 2.29)
cl.l<- c(1.37, 1.94, 1.90, 1.95, 2.06)
cl.u<- c(1.84, 2.41, 2.51, 2.56, 2.52)
emm.plot<-data.frame(cbind(trt, emm, cl.l, cl.u))
emm.plot$emm<-as.numeric(emm)
emm.plot$cl.l<-as.numeric(cl.l)
emm.plot$cl.u<-as.numeric(cl.u)


ggplot(emm.plot, aes(x=trt, y=emm))+ 
  ylim(0,2.75)+
  geom_errorbar(aes(x=trt, ymin=cl.l, ymax=cl.u), width=0.25)+
  geom_point(color= "black", shape=18, size=5, stroke=1)+
  annotate("text", x=1, y=2, label="a")+
  annotate("text", x=2, y=2.7, label="b")+
  annotate("text", x=3, y=2.7, label="b")+
  annotate("text", x=4, y=2.7, label="b")+
  annotate("text", x=5, y=2.7, label="b")+
  ylab("Estimated Marginal Mean (mg)")+
  xlab("Treatment")+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(margin(1,0,0,0)),
        axis.title.y = element_text(margin(0,1,0,0)),
        text = element_text(size=14))

# in color
ggplot(emm.plot, aes(x=trt, y=emm))+ 
  ylim(0,2.75)+
  geom_errorbar(aes(x=trt, ymin=cl.l, ymax=cl.u), width=0.25)+
  geom_point(color= "black", fill=agar_colors, shape=23, size=5, stroke=1)+
  annotate("text", x=1, y=2, label="a")+
  annotate("text", x=2, y=2.7, label="b")+
  annotate("text", x=3, y=2.7, label="b")+
  annotate("text", x=4, y=2.7, label="b")+
  annotate("text", x=5, y=2.7, label="b")+
  ylab("Estimated Marginal Mean (mg)")+
  xlab("")+
  theme(axis.text.x = element_text(size=11),
        axis.text.y = element_text(size=11),
        axis.title.y = element_text(margin(0,0.5,0,0)),
        text = element_text(size=12))
