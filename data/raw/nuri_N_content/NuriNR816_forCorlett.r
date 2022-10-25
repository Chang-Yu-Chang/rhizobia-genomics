

library(readr)
library(tidyverse)
library(ggplot2) 
library(car) 
library(lattice)
library(lme4)
library(lsmeans)
library(glmmADMB)
library(DHARMa)
library(emmeans)
library(devtools)
library(glmmTMB)
library(performance)
library(janitor)
library(dplyr)
library(ISLR)
library(mgcv)
library(visreg)
library(ggeffects)
library(gratia)

# Set contrasts (sum-to-zero rather than R's default treatment contrasts)
# http://rstudio-pubs-static.s3.amazonaws.com/65059_586f394d8eb84f84b1baaf56ffb6b47f.html
options(contrasts=c("contr.sum", "contr.poly"))

setwd("~/Library/CloudStorage/Box-Box")


# Read in full data
df = read.csv("nurirunfullfinal.csv")
xtabs(~Number, df)
dim(df) #194 total items, 42 columns
dim(filter(df, !is.na(Shoot_biomass))) #147 items
dim(filter(df, !is.na(Root_biomass))) #147 items?

#convert to factors (categories)
df$Number <- as.factor(df$Number)
df$Nitrogen <- as.factor(df$Nitrogen)
df$Rhizobia <- as.factor(df$Rhizobia)
df$Run_date <- as.factor(df$Run_date)
df$Rstate <-as.factor(df$Rstate)
df$Fixstate <-as.factor(df$Fixstate)

#convert to numerics so R will recognize as #s, exclude NAs
df$N_percent <-as.numeric(df$N_percent, na.rm=TRUE)
df$C_percent <-as.numeric(df$C_percent, na.rm=TRUE)
df$Total_nod <-as.numeric(df$Total_nod, na.rm=TRUE)
df$Egg_sacs <-as.numeric(df$Egg_sacs, na.rm=TRUE)
df$Galls <-as.numeric(df$Galls, na.rm=TRUE)
df$Shoot_biomass <-as.numeric(df$Shoot_biomass, na.rm=TRUE)
df$Root_biomass <-as.numeric(df$Root_biomass, na.rm=TRUE)  

onetr = read.csv("1randomtechrepNRfinal.csv")
view(onetr)
dim(onetr) #158x42
#convert to factors (categories)
onetr$Number <- as.factor(onetr$Number)
onetr$Nitrogen <- as.factor(onetr$Nitrogen)
onetr$Rhizobia <- as.factor(onetr$Rhizobia)
onetr$Run_date <- as.factor(onetr$Run_date)
onetr$Rstate <-as.factor(onetr$Rstate)
onetr$Fixstate <-as.factor(onetr$Fixstate)

#convert to numerics so R will recognize as #s, exclude NAs
onetr$N_percent <-as.numeric(onetr$N_percent, na.rm=TRUE)
onetr$C_percent <-as.numeric(onetr$C_percent, na.rm=TRUE)
onetr$Total_nod <-as.numeric(onetr$Total_nod, na.rm=TRUE)
onetr$Egg_sacs <-as.numeric(onetr$Egg_sacs, na.rm=TRUE)
onetr$Galls <-as.numeric(onetr$Galls, na.rm=TRUE)
onetr$Shoot_biomass <-as.numeric(onetr$Shoot_biomass, na.rm=TRUE)
onetr$Root_biomass <-as.numeric(onetr$Root_biomass, na.rm=TRUE)

onetr$Total_biomass<- onetr$Root_biomass+onetr$Shoot_biomass
onetr=filter(onetr, !is.na(Number))

onetr$scaledN <- scale(onetr$N_percent)
onetr$scaledNsq <- onetr$scaledN^2
onetr$scaled_biomass <- scale(onetr$Total_biomass)

dim(filter(onetr, !is.na("Shoot_biomass"))) #158x42
dim(filter(onetr, !is.na("Root_biomass"))) #158x42

alive <- onetr[onetr$Mortality_End=='1', ]
tabyl(onetr, Nitrogen, Rhizobia) #all plants alive, regardless of contamination status

xtabs(~Nitrogen+Rhizobia, filter(onetr, !is.na("Total_nod")))

#Making a dataframe for nodule analysis: excluding R- plants that formed nodules
#and also excluding #209, #210 (A145 originally) because they formed a lot of pink nodules
df0nod = rbind(onetr[onetr$Total_nod<1  & onetr$Fixstate=="R-", ], onetr[onetr$Fixstate!="R-", ])
dim(filter(df0nod, !is.na("Total_nod"))) #150x46 for plants that have a total nodule count
dim(filter(df0nod, !is.na("Number"))) #150x46 for plants that have a total nodule count
xtabs(~Nitrogen+Rhizobia, df0nod)
df0nod1=filter(df0nod, Number!=209)
df0nod1=filter(df0nod1, Number!=210)

dim(df0nod) #150x46 - WHY IS IT DECREASING??? FOR WHATEVER REASON THERE ARE STILL NAs that won't go away.
dim(df0nod1) #139x46 #but at least this seems right. Should actually be 139 vs 141?

dim(onetr) #158x46
dim(df0nod) #152x46: good.
#oh no, actually want 148x46 to get rid of 209, 210

dim(filter(df0nod1, !is.na(N_percent))) #if you exclude plants without N_percent vals, down to 98x46
dim(filter(df0nod1, !is.na(Shoot_biomass))) #101x46 for plants that have shoot biomass
dim(filter(df0nod1, !is.na(Root_biomass))) #101x46 for plants with root biomass.
dim(filter(df0nod1, !is.na(Total_biomass))) #101x46 for plants with total biomass. Hopefully this means all shoots also had roots

#make a new category and a new dataframe for R- plants that did form a few nodules
#this dataframe should be when nodule counts/types matter
testdf<-onetr
testdf$Rhizobia<-as.character(testdf$Rhizobia)
Rplus1<-filter(testdf, testdf$Rhizobia=="R-" & testdf$Total_nod>=1)
Rplus1$Rhizobia<-'Rplus'
dim(Rplus1) #8x46
xtabs(~Nitrogen+Rhizobia, Rplus1)
#Plants 209, 210 (A145 originally) should be moved into this undefined Rplus category as well
#Because they formed a high number of pink nodules: assume contamination with a different R strain.
Rplus2<-rbind(testdf[testdf$Number==209, ], testdf[testdf$Number==210, ])
dim(Rplus2) #2x46, as it should be
Rplus2$Rhizobia<-'Rplus'
#bind Rplus1 and Rplus2 for Rplus combined
Rplus<-rbind(Rplus1, Rplus2)
dim(Rplus) #10x46 as it should be

#Rplusmerged is a full dataframe, just reclassified into Rplus (for contaminated R- and A145 pink nod outliers)
Rplusmerged<-rbind(Rplus, df0nod)
Rplusmerged$Rhizobia<-as.factor(Rplusmerged$Rhizobia)
dim(Rplusmerged) #160=150+10!
xtabs(~Nitrogen+Rhizobia, Rplusmerged)

#exclude samples with extremely high N% — only use this dataframe for N_percent analyses. Not for nodule counts
#so make new versions of Rplusmerged and df0nod1
#df0nod1--> df0nodN
df0nodN<-filter(df0nod, !is.na(N_percent))
df0nodN<-df0nodN[df0nodN$N_percent<10,] 
dim(df0nodN) #dim now 98x46
xtabs(~Nitrogen+Rhizobia, df0nodN)

#Rplusmerged --> RplusmergedN
Rplusmerged<-filter(Rplusmerged, !is.na(N_percent))
dim(Rplusmerged) #dim 109x46
RplusmergedN<-Rplusmerged[Rplusmerged$N_percent<10,] 
dim(RplusmergedN) #dim now 107x46
xtabs(~Nitrogen+Rhizobia, RplusmergedN)

#I don't think I need any of this now, but keeping just in case
df5 = rbind(df4[df4$Fixstate!="Fix+", ]) #df5 has only non-contaminated R- plants and Fix- plants (no Fix+ plants), not that I need that right now
df6 = rbind(df3[df3$Fixstate!="Fix+", ]) #df6 has only R- and Fix- plants (no contamination exclusion)
df7 = rbind(df3[df3$Fixstate!="R-", ]) #df6 has only Fix+ and Fix- plants (no contamination exclusion needed)
df8 = rbind(df3[df3$Fixstate!="R-", ]) #df8 just excludes all R- plants

RplusmergedN$Rhizobia <- factor(RplusmergedN$Rhizobia, levels=c("R-", "Rplus", "A145", "Rm41", "Em1021", "Em1022"))
ggplot(data=RplusmergedN, aes(x=N_percent, y=Egg_sacs)) + 
  geom_smooth(aes(group=Rhizobia, color=Rhizobia), method = "loess", se = F, span=1.75, size=1.5) + 
  geom_point(aes(color=Rhizobia), size=2.1) + ylim(-0.5,36) +
  scale_color_manual(values = c(colors()[308], 'lightblue', 'cyan3', 'cyan4', 'steelblue', 'blue'))  + theme_bw() + theme(panel.grid = element_blank(), axis.text = element_text(color="black", size=10))
ggplot(data=RplusmergedN, aes(x=N_percent, y=Egg_sacs)) + 
  geom_smooth(method = "loess", se = T, span=1, size=2, color="black") + 
  geom_point(aes(color=Rhizobia), size=2.1) + ylim(-0.5,36) +
  scale_color_manual(values = c(colors()[308], 'lightblue', 'cyan3', 'cyan4', 'steelblue', 'blue'))  + theme_bw() + theme(panel.grid = element_blank(), axis.text = element_text(color="black", size=10))

ggplot(data=RplusmergedN, aes(x=Total_nod, y=Egg_sacs)) + 
  geom_smooth(method = "loess", se = T, span=10, size=1, color="black") + 
  geom_point(aes(color=Rhizobia), size=1.5) + ylim(-0.5,36) +
  scale_color_manual(values = c(colors()[308], 'lightblue', 'cyan3', 'cyan4', 'steelblue', 'blue'))  + theme_bw() + 
  theme(panel.grid = element_blank(), axis.text = element_text(color="black", size=10)) + 
  facet_grid(~Rhizobia)

#OK, now will fit GAM.
mod_gam = gam(Galls ~ s(N_percent), data = RplusmergedN, family=nb(theta=NULL, link="log")) #works at simplest
summary(mod_gam) #N_percent p-val 0.0519
#rsq = 0.0694, deviance explains 6.54%, -REML=336.42, scale est=1, n=107
plot(mod_gam)
pairs(emmeans(mod_gam, specs="Rhizobia")) #no differences by rhizobia strain alone
mod_gam$Rhizobia <- factor(mod_gam$Rhizobia, levels=c("R-", "Rplus", "A145", "Rm41", "Em1021", "Em1022"))
AIC(mod_gam) #674.2951
summary(mod_gam)$sp.criterion #GCV = 336.4181 
summary(mod_gam)$r.sq #0.0693558

RplusmergedN$Rhizobia <- factor(RplusmergedN$Rhizobia, levels=c("R-", "Rplus", "A145", "Rm41", "Em1021", "Em1022"))
mod_gam1 = gam(Galls ~ s(N_percent, by=Rhizobia) + s(N_percent) + Rhizobia, data = RplusmergedN, family=nb(theta=NULL, link="log")) #works at simplest
Anova(mod_gam1, type=3)
summary(mod_gam1) #interaction effect btween N_percent and R-? (0.0733)
#rsq = 0.0677, deviance explains 16.4%, -REML=334.9, scale est=1, n=107
plot(simulateResiduals(mod_gam1)) #can't yet fit extended.family from mgcv for dharma (see: github.com/florianhartig/DHARMa/issues/11)
plot(mod_gam1)
mod_gam1$Rhizobia <- factor(mod_gam1$Rhizobia, levels=c("R-", "Rplus", "A145", "Rm41", "Em1021", "Em1022"))
visreg(mod_gam1, xvar = "N_percent",
       by = "Rhizobia", data = RplusmergedN,
       method = "REML")
AIC(mod_gam1) #681.6212
summary(mod_gam1)$sp.criterion #GCV = 334.9029
summary(mod_gam1)$r.sq #0.067
pairs(emmeans(mod_gam1, specs="Rhizobia")) #no differences by rhizobia strain alone

mod_gam12=gam(Egg_sacs ~s(N_percent*Shoot_biomass, by=Rhizobia),  data = RplusmergedN, family=nb(theta=NULL, link="log"))
summary(mod_gam12)
plot(mod_gam12)

mod_gam2 = gam(Egg_sacs ~ s(N_percent, by=Rhizobia) + s(N_percent) + Rhizobia + s(Total_biomass), data = RplusmergedN, family=nb(theta=NULL, link="log")) #works at simplest
summary(mod_gam2) 
plot(mod_gam2)
#N_percent:R- 0.01574, N_percent*A145 0.01922, N_percent*Em1021 0.00896, Root_biomass <2e-16
#Rsq = 0.371, -REML=313.56, scale est. =1, n=107, deviance explained = 46.9%
visreg(mod_gam2, xvar = "N_percent",
       by = "Rhizobia", data = RplusmergedN,
       method = "REML")
AIC(mod_gam2) #633.3074
summary(mod_gam2)$sp.criterion #GCV = 315.9881
summary(mod_gam2)$r.sq #0.411693
pairs(emmeans(mod_gam2, specs="Rhizobia")) #no differences by rhizobia strain alone

mod_gam3 = gam(Galls ~ s(N_percent, by=Rhizobia) + s(Total_biomass) +s(Total_nod, by=Rhizobia), data = RplusmergedN, family=nb(theta=NULL, link="log")) #works at simplest
summary(mod_gam3) 
#N_percent:R- 0.01574, N_percent*A145 0.01922, N_percent*Em1021 0.00896, Root_biomass <2e-16
#Rsq = 0.371, -REML=313.56, scale est. =1, n=107, deviance explained = 46.9%
visreg(mod_gam3, xvar = "N_percent",
       by = "Rhizobia", data = RplusmergedN,
       method = "REML")
AIC(mod_gam3) #630.8822
summary(mod_gam3)$sp.criterion #GCV = 315.7649
summary(mod_gam3)$r.sq #0.4199654
pairs(emmeans(mod_gam3, specs="Rhizobia")) #no differences by rhizobia strain alone

plot(ggeffects::ggpredict(mod_gam3), facets = TRUE)
gratia::draw(mod_gam3)

#visualizing GAMS
plot(ggeffects::ggpredict(mod_gam1), facets = TRUE)
gratia::draw(mod_gam1)

mod_gam1$Rhizobia <- factor(mod_gam1$Rhizobia, levels=c("R-", "Rplus", "A145", "Rm41", "Em1021", "Em1022"))
vis.gam(mod_gam1, view=c('Rhizobia', 'N_percent'), type = 'response', plot.type = 'persp', color="topo", phi = 30, theta = 30, n.grid = 500, border = NA)

mod_gam2$Rhizobia <- factor(mod_gam2$Rhizobia, levels=c("R-", "Rplus", "A145", "Rm41", "Em1021", "Em1022"))
vis.gam(mod_gam2, view=c('Rhizobia', 'N_percent'), type = 'response', plot.type = 'persp', color="topo", phi = 30, theta = 30, n.grid = 500, border = NA)
vis.gam(mod_gam2, type = 'response', plot.type = 'persp', phi = 30, theta = 30, n.grid = 500, border = NA)

mod_gam3$Rhizobia <- factor(mod_gam3$Rhizobia, levels=c("R-", "Rplus", "A145", "Rm41", "Em1021", "Em1022"))
vis.gam(mod_gam3, view=c('Rhizobia', 'N_percent'), type = 'response', plot.type = 'persp', color="topo", phi = 30, theta = 30, n.grid = 500, border = NA)
vis.gam(mod_gam3, type = 'response', plot.type = 'persp', phi = 30, theta = 30, n.grid = 500, border = NA)

glmmTMB = glmmTMB(Galls ~ N_percent, RplusmergedN, family="nbinom1")
plot(simulateResiduals(glmmTMB))
summary(glmmTMB)
glm = glm(Galls ~ N_percent, RplusmergedN, family="nbinom1")
#won't implement nbinom1 for glm... for DHARMa, but the glmmTMB looks ok
summary(glm)
#trying nb notation
glmnb=glm.nb(Galls~N_percent, RplusmergedN)
plot(simulateResiduals(glmnb)) #oh that looks ok too
summary(glmnb) #AIC: 678.27
with(summary(glmnb), 1 - deviance/null.deviance) #McFadden r-sq = 0.0000979086. bad.

glmmTMB1 = glmmTMB(Galls ~ N_percent*Rhizobia, RplusmergedN, family="nbinom1")
plot(simulateResiduals(glmmTMB1))
summary(glmmTMB1)
glm1 = glm(Galls ~ N_percent*Rhizobia, RplusmergedN, family="nbinom1")
#won't implement nbinom1 for glm... for DHARMa, but the glmmTMB looks ok
summary(glm1)
#trying nb notation
glmnb1=glm.nb(Galls~N_percent*Rhizobia, RplusmergedN)
plot(simulateResiduals(glmnb1)) #oh that looks ok too
summary(glmnb1)
with(summary(glmnb1), 1 - deviance/null.deviance) #McFadden r-sq = 0.07 ?

glmmTMB2 = glmmTMB(Galls ~ N_percent*Rhizobia + Root_biomass, RplusmergedN, family="nbinom1") 
plot(simulateResiduals(glmmTMB2)) #won't implement nbinom1 for glm for DHARMa, but equiv glmmTMB looks ok
summary(glmmTMB2)
glm2 = glm(Galls ~ N_percent*Rhizobia + Root_biomass, RplusmergedN, family="nbinom1")
summary(glm2)
#trying nb notation
glmnb2=glm.nb(Galls~N_percent*Rhizobia + Root_biomass, RplusmergedN)
plot(simulateResiduals(glmnb2)) #that looks pretty curved.... for model predictions
summary(glmnb2)
with(summary(glmnb2), 1 - deviance/null.deviance) #McFadden r-sq = 0.30258 ?

visreg(glmmTMB1, xvar = "N_percent",
       by = "Rhizobia", data = RplusmergedN,
       method = "REML")
visreg(glmmTMB2, xvar = "N_percent",
       by = "Rhizobia", data = RplusmergedN,
       method = "REML")

glmmTMB3 = glmmTMB(Galls ~ N_percent*Rhizobia + Root_biomass + Total_nod, RplusmergedN, family="nbinom1") 
plot(simulateResiduals(glmmTMB3)) #won't implement nbinom1 for glm for DHARMa, but equiv glmmTMB looks ok
summary(glmmTMB3)
glm3 = glm(Galls ~ N_percent*Rhizobia + Root_biomass + Total_nod, RplusmergedN, family="nbinom1")
summary(glm3)
#trying nb notation
glmnb3=glm.nb(Galls~N_percent*Rhizobia + Root_biomass + Total_nod, RplusmergedN)
plot(simulateResiduals(glmnb3)) #that looks pretty curved.... for model predictions
summary(glmnb3)
with(summary(glmnb3), 1 - deviance/null.deviance) #McFadden r-sq = 0.3305825

anova(mod_gam1, glmnb1, test="Chisq")
anova(mod_gam2, glmnb2, test="Chisq")
anova(mod_gam3, glmnb3, test="Chisq")

view(RplusmergedN)

#for egg sacs:
glmnb=glm.nb(Egg_sacs~N_percent, RplusmergedN)
plot(simulateResiduals(glmnb)) #oh that looks ok too
summary(glmnb) #AIC: 678.27
with(summary(glmnb), 1 - deviance/null.deviance) #McFadden r-sq = 0.0000979086. bad.

glmnb1=glm.nb(Galls~N_percent*Rhizobia, RplusmergedN)
plot(simulateResiduals(glmnb1)) #oh that looks ok too
summary(glmnb1)
with(summary(glmnb1), 1 - deviance/null.deviance) #McFadden r-sq = 0.07 ?

glmnb2=glm.nb(Galls~N_percent*Rhizobia + Total_biomass, RplusmergedN)
plot(simulateResiduals(glmnb2)) #that looks pretty curved.... for model predictions
summary(glmnb2)
with(summary(glmnb2), 1 - deviance/null.deviance) #McFadden r-sq = 0.30258 ?

glmnb3=glm.nb(Galls~N_percent*Rhizobia + Total_biomass + Total_nod, RplusmergedN)
plot(simulateResiduals(glmnb3)) #that looks pretty curved.... for model predictions
summary(glmnb3)
with(summary(glmnb3), 1 - deviance/null.deviance) #McFadden r-sq = 0.3305825

#now trying out residuals:
gmod = glmmTMB(Galls ~ Nitrogen*Rhizobia + Root_biomass, RplusmergedN, family="nbinom1") 
gmod = glm.nb(Galls ~ Nitrogen*Rhizobia + Root_biomass, RplusmergedN)
plot(simulateResiduals(gmod))
summary(gmod)
pairs(emmeans(gmod, specs="Rhizobia")) #no differences by rhizobia strain alone
length(residuals(gmod, type = c("response")))

RplusmergedN$gresid<-residuals(gmod, type=c("response"))
view(RplusmergedN)

emod = glmmTMB(Egg_sacs ~ Nitrogen*Rhizobia + Root_biomass, RplusmergedN, family="nbinom1") 
emod = glm.nb(Egg_sacs ~ Nitrogen*Rhizobia + Root_biomass, RplusmergedN)
plot(simulateResiduals(emod))
summary(emod)
pairs(emmeans(emod, specs="Rhizobia")) #no differences by rhizobia strain alone
length(residuals(emod, type = c("response")))

RplusmergedN$eresid<-residuals(emod, type=c("response"))
view(RplusmergedN)

nmod = glmmTMB(Total_nod ~ Nitrogen*Rhizobia + Root_biomass, RplusmergedN) 
nmod = glm.nb(Total_nod ~ Nitrogen*Rhizobia + Root_biomass, RplusmergedN) #either work but glmmTMB is slightly better
plot(simulateResiduals(nmod))
summary(nmod)
pairs(emmeans(nmod, specs="Rhizobia")) #no differences by rhizobia strain alone
length(residuals(nmod, type = c("response")))

RplusmergedN$nresid<-residuals(nmod, type=c("response"))
view(RplusmergedN)

RplusmergedN$Rhizobia <- factor(RplusmergedN$Rhizobia, levels=c("R-", "Rplus", "A145", "Rm41", "Em1021", "Em1022"))
ggplot(data=RplusmergedN, aes(x=nresid, y=gresid)) + 
  geom_smooth(method = lm, se = T, size=1, color="black") + 
  geom_point(aes(color=Rhizobia), size=1.5) +
  scale_color_manual(values = c(colors()[308], 'lightblue', 'cyan3', 'cyan4', 'steelblue', 'blue'))  + theme_bw() + 
  theme(panel.grid = element_blank(), axis.text = element_text(color="black", size=10))
ggplot(data=RplusmergedN, aes(x=nresid, y=gresid)) + 
  geom_smooth(aes(group=Rhizobia, color=Rhizobia), method = lm, se = F, size=1) + 
  geom_point(aes(color=Rhizobia), size=1.5) +
  scale_color_manual(values = c(colors()[308], 'lightblue', 'cyan3', 'cyan4', 'steelblue', 'blue'))  + theme_bw() + 
  theme(panel.grid = element_blank(), axis.text = element_text(color="black", size=10))
ggplot(data=RplusmergedN, aes(x=nresid, y=gresid)) + 
  geom_smooth(aes(group=Rhizobia, color=Rhizobia), method = lm, se = T, size=1) + 
  geom_point(aes(color=Rhizobia), size=1.5) +
  scale_color_manual(values = c(colors()[308], 'lightblue', 'cyan3', 'cyan4', 'steelblue', 'blue'))  + theme_bw() + 
  theme(panel.grid = element_blank(), axis.text = element_text(color="black", size=10)) + facet_grid(~Rhizobia)

RplusmergedN$Rhizobia <- factor(RplusmergedN$Rhizobia, levels=c("R-", "Rplus", "A145", "Rm41", "Em1021", "Em1022"))
ggplot(data=RplusmergedN, aes(x=nresid, y=eresid)) + 
  geom_smooth(aes(group=Rhizobia, color=Rhizobia), method = lm, se = F, size=1) + 
  geom_point(aes(color=Rhizobia), size=1.5) +
  scale_color_manual(values = c(colors()[308], 'lightblue', 'cyan3', 'cyan4', 'steelblue', 'blue'))  + theme_bw() + 
  theme(panel.grid = element_blank(), axis.text = element_text(color="black", size=10)) + facet_grid(~Nitrogen)

ggplot(data=RplusmergedN, aes(x=nresid, y=eresid)) + 
  geom_smooth(aes(group=Rhizobia, color=Rhizobia), method = lm, se = F, size=1) + 
  geom_point(aes(color=Rhizobia), size=1.5) +
  scale_color_manual(values = c(colors()[308], 'lightblue', 'cyan3', 'cyan4', 'steelblue', 'blue'))  + theme_bw() + 
  theme(panel.grid = element_blank(), axis.text = element_text(color="black", size=10)) + facet_grid(~Rhizobia)

#try looking just at the 2.5N treatment
#do filter function
gmod2 = glmmTMB(Galls ~ Rhizobia*Nitrogen + (1|Gall_Counter) + (1|Block), RplusmergedN, family="nbinom1") 
Anova(gmod2, type=3)
summary(gmod2)
plot(simulateResiduals(gmod2))
pairs(emmeans(gmod2, specs="Rhizobia","Nitrogen"))
pairs(emmeans(gmod2, specs="Nitrogen","Rhizobia"))
pairs(emmeans(ogmod3, specs="Rhizobia"))
pairs(emmeans(ogmod3, specs="Nitrogen"))
emmeans(ogmod3, specs = "Rhizobia", "Nitrogen", type="response")
emm.ogmod3 <- emmeans(ogmod3, specs = "Rhizobia", "Nitrogen", type="response") #backtransform
emm.ogmod3.df <- as.data.frame(emm.ogmod3)
view(emm.ogmod3.df)

#plot of galls adjusted for biomass (from model)
RplusmergedN$Rhizobia <- factor(RplusmergedN$Rhizobia, levels=c("R-", "Rplus", "A145","Rm41", "Em1021", "Em1022"))
emm.ogmod3.df$Rhizobia <- factor(emm.ogmod3.df$Rhizobia, levels=c("R-", "Rplus", "A145","Rm41", "Em1021", "Em1022"))
#Rminusdf=RplusmergedN[df4removed$Rhizobia=='R-' & df4removed$Nitrogen=="2.5",]
#Rminusdf
ggplot() + geom_point(data=RplusmergedN, 
                      aes(x=Nitrogen, y=Galls, fill=Rhizobia), alpha=0.5, position=position_jitterdodge(jitter.width = 0.3, dodge.width=0.8), pch=21) +
  geom_errorbar(data=emm.ogmod3.df, aes(x=Nitrogen, ymin=lower.CL, ymax=upper.CL, group=Rhizobia), 
                width=0, size=0.75, position=position_dodge(width=0.5), lwd=0.5) +
  geom_point(data=emm.ogmod3.df, aes(x=Nitrogen, y=response, fill=Rhizobia), size=5, pch=21, stroke=1, colour="black", position=position_dodge(width=0.5)) +
  scale_colour_manual(values = c(colors()[308], 'skyblue', 'cyan3', 'cyan4', 'steelblue', 'blue')) +
  scale_fill_manual(values = c(colors()[308],'skyblue', 'cyan3', 'cyan4', 'steelblue', 'blue')) +
  ylab("Galls") +
  theme_bw() + theme(panel.grid = element_blank(), axis.text = element_text(color="black", size=18), axis.title.x=element_blank(), axis.title.y=element_blank()) + ylim(-0.5,36)



dim(filter(RplusmergedN, Nitrogen=="2.5")) #38x49
highN <- filter(RplusmergedN, Nitrogen=="2.5")

gmodN = glmmTMB(Galls ~ Total_nod*Rhizobia + Root_biomass, highN, family="nbinom1") 
Anova(gmodN, type=3)
summary(gmodN)
plot(simulateResiduals(gmodN))
pairs(emmeans(gmodN, specs="Rhizobia"))
emmeans(ogmod3, specs = "Rhizobia", "Nitrogen", type="response")
emm.ogmod3 <- emmeans(ogmod3, specs = "Rhizobia", "Nitrogen", type="response") #backtransform
emm.ogmod3.df <- as.data.frame(emm.ogmod3)
view(emm.ogmod3.df)

ggplot() + geom_point(data=highN, 
                      aes(x=Nitrogen, y=Galls, fill=Rhizobia), alpha=0.5, position=position_jitterdodge(jitter.width = 0.3, dodge.width=0.8), pch=21) +
  geom_errorbar(data=emm.ogmod3.df, aes(x=Nitrogen, ymin=lower.CL, ymax=upper.CL, group=Rhizobia), 
                width=0, size=0.75, position=position_dodge(width=0.5), lwd=0.5) +
  geom_point(data=emm.ogmod3.df, aes(x=Nitrogen, y=response, fill=Rhizobia), size=5, pch=21, stroke=1, colour="black", position=position_dodge(width=0.5)) +
  scale_colour_manual(values = c(colors()[308], 'skyblue', 'cyan3', 'cyan4', 'steelblue', 'blue')) +
  scale_fill_manual(values = c(colors()[308],'skyblue', 'cyan3', 'cyan4', 'steelblue', 'blue')) +
  ylab("Galls") +
  theme_bw() + theme(panel.grid = element_blank(), axis.text = element_text(color="black", size=18), axis.title.x=element_blank(), axis.title.y=element_blank()) + ylim(-0.5,36)

