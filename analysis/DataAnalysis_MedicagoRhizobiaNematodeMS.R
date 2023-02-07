#--------------------------------------------------------------#
# Data analysis for the manuscript "Genetic conflict with a
# parasitic nematode disrupts the legume-rhizobia mutualism"

#   Authors:  Corlett Wood, Bonnie Pilkington, Priya Vaidya,
#             Caroline Biel, and John Stinchcombe

#   Corresponding author: Corlett Wood (corlett.wood@utoronto.ca)
#--------------------------------------------------------------#



#--------------------------------------------------------------#
# Prepare workspace
#--------------------------------------------------------------#
# Load libraries
library(car)
library(lattice)
library(lme4)
library(lsmeans)
library(glmmADMB)
# Load script with custom functions to check assumptions of linear models and calculate VIFs
source("~/Documents/Post-doc/Local adaptation/checkAssumptions.R")

# Set working directory
  setwd("~/Documents/Post-doc/Medicago-Nematode MS/Data and code for Dryad/")

# Set contrasts (sum-to-zero rather than R's default treatment contrasts)
  # http://rstudio-pubs-static.s3.amazonaws.com/65059_586f394d8eb84f84b1baaf56ffb6b47f.html
  options(contrasts=c("contr.sum", "contr.poly"))

# Read in the data: EXPERIMENT 1
  experiment1 = read.csv("Experiment1.csv")

# Read in the data: EXPERIMENT 2
  experiment2 = read.csv("Experiment2.csv")



#--------------------------------------------------------------#
# Data analysis: EXPERIMENT 1
#--------------------------------------------------------------#
# THE EFFECT OF NODULES & GALLS ON FITNESS IN CO-INFECTED PLANTS ####
# Aboveground biomass
  mod = lmer(log(AbovegroundMass.g) ~ Nodules + Galls + RootMass.g + scale(Number.of.Nematode.Eggs)
             + Researcher + (1|Block), experiment1)
  Anova(mod, type=3)
  summary(mod)
  check.assumptions(mod) # Assumptions met
  vif.lme(mod)
  # Excluding plants that received too many nematode eggs does not qualitatively change the answer
  mod = lmer(AbovegroundMass.g ~ Nodules + Galls + RootMass.g + scale(Number.of.Nematode.Eggs)
             + Researcher + (1|Block), subset(experiment1, Number.of.Nematode.Eggs<=400))

# Fruit mass
  mod = lmer(FruitMass.g ~ Nodules + Galls + RootMass.g + scale(Number.of.Nematode.Eggs)
             + Researcher + (1|Block), experiment1)
  Anova(mod, type=3)
  summary(mod)
  check.assumptions(mod) # Assumptions met
  vif.lme(mod)
  # Excluding plants that received too many nematode eggs does not qualitatively change the answer
  mod = lmer(FruitMass.g ~ Nodules + Galls + RootMass.g + scale(NumEggs) + Researcher + (1|Block),
             subset(experiment1, Number.of.Nematode.Eggs <=400))


# GENETIC VARIATION IN INFECTIVITY AMONG NEMATODE LINES
# Test for variation among nematode lines
  mod2=lmer(log(Galls+1) ~ (1|NemaGenotype) + (1|PlantGenotype) + (1|Block) + (1|NemaGenotype:PlantGenotype)
            + Researcher + scale(Number.of.Nematode.Eggs) + RootMass.g, experiment1)
  # Neither a Poisson or a negative binomial GLMM fit the data well
    #mod2=glmer(Galls ~ (1|NemaGenotype) + (1|PlantGenotype) + (1|Block) + (1|NemaGenotype:PlantGenotype)
    #           + scale(Number.of.Nematode.Eggs) + RootMass.g, experiment1, family="poisson")
    #mod2=glmmadmb(Galls ~ (1|NemaGenotype) + (1|PlantGenotype) + (1|Block) + (1|NemaGenotype:PlantGenotype)
    #              + scale(Number.of.Nematode.Eggs) + RootMass.g, data = subset(experiment1, !is.na(RootMass.g)),
    #              family="nbinom", zeroInflation = T)
  check.assumptions(mod2)
  sum(residuals(mod2, type="pearson")^2) / df.residual(mod2) # Not overdispersed

# Log-likelihood ratio tests for genetic variation
  # Genetic variation in infectivity among nematode genotypes
  nonema=update(mod2, .~. - (1|NemaGenotype))
  anova(mod2, nonema)

  # Genetic variation in resistance among plant genotypes
  noplant=update(mod2, .~. - (1|PlantGenotype))
  anova(mod2, noplant)

# Exclude the plant-nematode combinations with fewer than three replicates
  experiment1=within(experiment1, {NemaPlant = as.factor(paste(NemaGenotype, PlantGenotype, "_"))})
  reps=data.frame(xtabs(~NemaPlant, experiment1))
  keeps = droplevels(reps[which(reps$Freq>2), "NemaPlant"])

  # Re-run the model without the plant-nematode combinations with fewer than three replicates
  mod2=lmer(log(Galls+1) ~ (1|NemaGenotype) + (1|PlantGenotype) + (1|Block) + (1|NemaGenotype:PlantGenotype)
            + Researcher + scale(Number.of.Nematode.Eggs) + RootMass.g,
            subset(experiment1, NemaPlant %in% keeps))

  # Log-likelihood ratio test for a genotype-by-genotype interaction
  nogxg=update(mod2, .~. - (1|NemaGenotype:PlantGenotype))
  anova(mod2, nogxg)



#--------------------------------------------------------------#
# Data analysis: EXPERIMENT 2
#--------------------------------------------------------------#

#--------------------------------------------------------------#
# GENETIC VARIATION IN PLANT SUSCEPTIBILITY TO NEMATODES
  # Fit the model
  mod.galls.zinb=glmmadmb(Galls ~ (1|PlantGenotype) + (1|Block) + RootMass.g + Researcher,
                          data = subset(experiment2, Treatment=="RN"), family="nbinom", zeroInflation = T)
  check.assumptions(mod.galls.zinb)
  sum(residuals(mod.galls.zinb, type="pearson")^2) / df.residual(mod.galls.zinb)
  # Test fixed effects
  Anova(mod.galls.zinb, type=3)
  # Test for variation among genotypes
  mod.galls.nogeno = update(mod.galls.zinb, .~. - (1|FinalGenotype))
  anova(mod.galls.nogeno, mod.galls.zinb)


#--------------------------------------------------------------#
# GENETIC CORRELATION BETWEEN NODULE AND GALL NUMBER
# GENETIC CORRELATION BETWEEN GALL NUMBER AND THE CHANGE IN NODULE NUMBER
  # Estimate genotype means (i.e., conditional modes or BLUPS) for galls
  mod.galls.zinb=glmmadmb(Galls ~ (1|PlantGenotype) + (1|Block) + RootMass.g + Researcher,
                          data = subset(experiment2, Treatment=="RN"), family="nbinom", zeroInflation = T)
    # Extract genotype means
  means.galls=as.data.frame(ranef(mod.galls.zinb)$PlantGenotype)
    # Extract confidence intervals
  CI.galls=as.data.frame(1.96*ranef(mod.galls.zinb, sd=T)$`1`)
    # Combine into a dataframe
  blups.galls=as.data.frame(cbind(means.galls, CI.galls))
  names(blups.galls) = c("Galls", "CI.galls")
  blups.galls$PlantGenotype = row.names(blups.galls)

  # Estimate genotype means (i.e., conditional modes or BLUPS) for nodules
  # -- Treatment: R (rhizobia-only)
  mod.nod.nb.R = glmmadmb(Nodules ~ (1|PlantGenotype) + (1|Block) + Researcher + RootMass.g,
                          subset(experiment2, Treatment == "R"),
                          family="nbinom", zeroInflation = T)
  check.assumptions(mod.nod.nb.R)
    # Extract genotype means
  means.nod.R=as.data.frame(ranef(mod.nod.nb.R)$PlantGenotype)
    # Extract confidence intervals
  CI.nod.R=as.data.frame(1.96*ranef(mod.nod.nb.R, sd=T)$`1`)
    # Combine into a data frame
  blups.nod.R=as.data.frame(cbind(means.nod.R, CI.nod.R))
  names(blups.nod.R) = c("nod.R", "CI.nod.R")
  blups.nod.R$PlantGenotype = row.names(blups.nod.R)

  # -- Treatment: RN (rhizobia and nematodes)
  mod.nod.nb.RN = glmmadmb(Nodules ~ (1|PlantGenotype) + (1|Block) + Researcher + RootMass.g,
                           subset(experiment2, Treatment == "RN" & Nodules < 200),
                           family="nbinom", zeroInflation = T)
  check.assumptions(mod.nod.nb.RN)
    # Extract genotype means
  means.nod.RN=as.data.frame(ranef(mod.nod.nb.RN)$PlantGenotype)
    # Extract confidence intervals
  CI.nod.RN=as.data.frame(1.96*ranef(mod.nod.nb.RN, sd=T)$`1`)
    # Combine into a data frame
  blups.nod.RN=as.data.frame(cbind(means.nod.RN, CI.nod.RN))
  names(blups.nod.RN) = c("nod.RN", "CI.nod.RN")
  blups.nod.RN$PlantGenotype = row.names(blups.nod.RN)

  # Merge nodule dataframes
  blups.nodules = merge(x=blups.nod.R, y=blups.nod.RN, by="PlantGenotype")
  # Merge nodule dataframe with gall dataframe
  blups.nodgall = merge(x=blups.nodules, y=blups.galls, by="PlantGenotype")
  # Calculate the difference in nodule number in R and RN treatments
  blups.nodgall = within(blups.nodgall, {
    R.minus.RN = nod.R-nod.RN  # Calculate difference
    CI.R.minus.RN = NA  # Create column to store the confidence interval on the difference
  })

  # Resample to create confidence intervals on the difference in nodule number between the two treatments
  # Sample one random point in the confidence interval for each one
  for(geno in 1:nrow(blups.nodgall)){
    R = with(blups.nodgall[geno,], rnorm(mean=nod.R, sd=(CI.nod.R/1.96), n=1000))
    RN = with(blups.nodgall[geno,], rnorm(mean=nod.RN, sd=(CI.nod.RN/1.96), n=1000))
    CI = sd(R-RN)*1.96
    blups.nodgall[geno, "CI.R.minus.RN"] = CI
  }


  # Calculate the correlation between nodule and gall formation
  # -- including the outlier genotype HM170
  with(blups.nodgall, cor.test(nod.R, Galls, method="spearman"))
  # -- excluding the outlier genotype HM170
  with(subset(blups.nodgall, PlantGenotype != "HM170"), cor.test(nod.R, Galls))

  # Correlation between gall formation and the difference in nodule number
  # --> Significant positive correlation (still significant without HM170)
  blups.nodgall = within(blups.nodgall, {deltaNod = nod.R-nod.RN})
  # -- including the outlier genotype HM170
  with(blups.nodgall, cor.test(Galls, deltaNod))
  # -- excluding the outlier genotype HM170
  with(subset(blups.nodgall, PlantGenotype != "HM170"), cor.test(Galls, deltaNod))




#--------------------------------------------------------------#
# EFFECT OF NEMATODES ON THE RHIZOBIA MUTUALISM
# Nodule number
  # Run the analysis
  mod.nod.nb = glmmadmb(Nodules ~ Treatment + (1|PlantGenotype) + (1|Treatment:PlantGenotype)
                        + (1|Block) + (1|Treatment:Block) + Researcher + RootMass.g,
                        experiment2, family="nbinom", zeroInflation = T)
  check.assumptions(mod.nod.nb)
  sum(residuals(mod.nod.nb, type="pearson")^2) / df.residual(mod.nod.nb)
    # -- This model fit poorly -- the residuals

  # Exclude 4 plants with >200 nodules
  # --> Model fits better and the results are qualitatively similar
  mod.nod.nb = glmmadmb(Nodules ~ Treatment + (1|PlantGenotype) + (1|Treatment:PlantGenotype)
                        + (1|Block) + (1|Treatment:Block) + Researcher + RootMass.g,
                        subset(experiment2, Nodules<200), family="nbinom", zeroInflation = T)
  check.assumptions(mod.nod.nb)
  sum(residuals(mod.nod.nb, type="pearson")^2) / df.residual(mod.nod.nb)

  # Test treatment by plant genotype effect
  View(anova(mod.nod.nb, update(mod.nod.nb, .~. -(1|Treatment:PlantGenotype))))
  # Test plant genotype effect
  View(anova(mod.nod.nb, update(mod.nod.nb, .~. -(1|PlantGenotype))))
  # Test treatment effect
  Anova(mod.nod.nb, type=3)


# Mean nodule mass
  mod.meannodmass = lmer(log(MeanNodMass) ~ Treatment + RootMass.g + (1|PlantGenotype)
                         + (1|Treatment:PlantGenotype) + (1|Block) + (1|Treatment:Block), experiment2)
  check.assumptions(mod.meannodmass)
  overdisp_fun(mod.meannodmass)

  # Test the treatment by plant genotype effect
  View(anova(mod.meannodmass, update(mod.meannodmass, .~. - (1|Treatment:PlantGenotype))))
  # Test the plant genotype effect
  View(anova(mod.meannodmass, update(mod.meannodmass, .~. - (1|PlantGenotype))))
  # Test the treatment effect
  Anova(mod.meannodmass, type=3)


# Total nodule biomass
  mod.nodmass = lmer(log(TotalNodMass) ~ Treatment + RootMass.g + (1|PlantGenotype)
                         + (1|Treatment:PlantGenotype) + (1|Block) + (1|Treatment:Block), experiment2)
  check.assumptions(mod.nodmass)
  overdisp_fun(mod.nodmass.ran)

  # Test the treatment by plant genotype interaction
  View(anova(mod.nodmass, update(mod.nodmass, .~. -(1|Treatment:PlantGenotype))))
  # Test the plant genotype effect
  View(anova(mod.nodmass, update(mod.nodmass, .~. -(1|PlantGenotype))))
  # Test the treatment effect
  Anova(mod.nodmass, type=3)


# Aboveground biomass
  mod.above = lmer(log(AbovegroundMass.g) ~ Treatment + (1|PlantGenotype) + (1|Treatment:PlantGenotype)
                   + (1|Block) + (1|Treatment:Block), experiment2)
  check.assumptions(mod.above)
  overdisp_fun(mod.above)

  # Test for a treatment by plant genotype interaction
  View(anova(mod.above, update(mod.above, .~. -(1|Treatment:PlantGenotype))))
  # Test for variation among plant genotypes
  View(anova(mod.above, update(mod.above, .~. -(1|PlantGenotype))))
  # Test the fixed effect of treatment
  Anova(mod.above, type=3)


# Flowering time
  # Subset the genotypes that have 3 replicates that flowered in both treatments
  fr = subset(experiment2, !is.na(FloweringTime))
  df = data.frame(xtabs(~PlantGenotype+Treatment, fr))
  df = subset(df, Freq>2)
  keeps=subset(data.frame(xtabs(~PlantGenotype, df)), Freq==2)$PlantGenotype

  # Run the analysis
  mod.flower = lmer(log(FloweringTime) ~ Treatment + (1|PlantGenotype) + (1|Treatment:PlantGenotype)
                    + (1|Block) + (1|Treatment:Block),
                    subset(experiment2, Flower.YesNo == "Y" & PlantGenotype %in% keeps))
  check.assumptions(mod.flower)
  overdisp_fun(mod.flower)

  # Test for a treatment by genotype interaction
  View(anova(mod.flower, update(mod.flower, .~. -(1|Treatment:PlantGenotype))))
  # Test for a genotype effect
  View(anova(mod.flower, update(mod.flower, .~. -(1|PlantGenotype))))
  # Test treatment effect
  Anova(mod.flower, type=3)


# Total fruit mass
  # Subset the genotypes that have 3 replicates that fruited in both treatments
  fr = subset(experiment2, !is.na(FruitMass.g))
  df = data.frame(xtabs(~PlantGenotype+Treatment, fr))
  df = subset(df, Freq>2)
  keeps=subset(data.frame(xtabs(~PlantGenotype, df)), Freq==2)$PlantGenotype

  # Random effects model
  mod.fruit = lmer(log(FruitMass.g+1) ~ Treatment + (1|PlantGenotype) + (1|Treatment:PlantGenotype)
                   + (1|Block) + (1|Treatment:Block),
                   subset(experiment2, Fruit.YesNo == "Y" & PlantGenotype %in% keeps))
  overdisp_fun(mod.fruit)
  check.assumptions(mod.fruit)

  # Test for a treatment by plant genotype interaction
  View(anova(mod.fruit, update(mod.fruit, .~. -(1|Treatment:PlantGenotype))))
  # Test for a genotype effect
  View(anova(mod.fruit, update(mod.fruit, .~. -(1|PlantGenotype))))
  # Test for an effect of treatment
  Anova(mod.fruit, type=3)
