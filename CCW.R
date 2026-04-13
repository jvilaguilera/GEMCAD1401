###########################################################################.
# Title: Script to Estimate the Per-Protocol Effect on Overall Survival   #
# of Targeted  Therapies at First Line vs Second Line in Patients with    #
# Metastatic Colorectal Cancer using the GEMCAD 1401 registry data        #
# Analysis: Clone-Censor-Weight Analysis                                  #
# Authors: Julia Vila Guilera, Manuel Zamparini, Xabier Garcia de Albeniz #
###########################################################################.

# --------------------------Load required libraries-------------------------
require(dplyr)
require(survival)
require(sandwich)
require(lmtest)
require(ggplot2)
require(survminer)
require(splitstackshape)
require(geepack)
require(mice)
require(purrr)
require(tidyr)
require(splines)
require(rms)
require(sqldf)
require(boot)

# -------------------------------Load data ----------------------------------

## Load data, call it ds
# This data should be cleaned and in the wide form. 
ds<-read.csv(file=path, header=T) 

## Structural variables

# subject_id - id that identifies individuals
# adherence1 - indicator for whether this individual adhered to treatment straetegy "mAB as first line"
# adherence2 - indicator for whether this individual adhered to treatment straetegy "mAB as second line"
# LastStatus - last known status of individual (dead or alive)
# DeathDate - date of death, if person died
# LastFupDate - date of last follow-up
# IndexDate - index date
# DateOfCensoring - Date when a clone stops adhering to their assigned strategy
# time - indicator for each week of follow-up in expanded (long format) data
# fupCI - number of weeks from index date until the development of a contraindication
# fupMAB1 - number of weeks from index date until initiation of mAB as first line
# fupMAB2 - number of weeks from index date until initiation of mAB as second line
# fupQMT2 - number of weeks from index date until initiation of second line chemotherapy
# MAB1_tv - binary indicator, only 1 the week MAB is initiated as first line
# MAB2_tv - binary indicator, only 1 the week MAB is initiated as second line

# timetoadherence - number of weeks from index date until the week of an event that establishes adherence to the assigned strategy occurring 
# adh_indicator - indicator variable, takes value of 0 for the weeks prior to an event establishing adherence to assigned strategy until end of follow-up  
#                                     takes value of 1 the week of the occurrence of an event that establishes adherence to assigned strategy until end of follow-up
#                                     takes value of 2 the weeks posterior to an event establishing adherence to assigned strategy until end of follow-up


## Baseline covariates
# age - age (in years)
# agesq - age squared (in years)
# PS_basal - ECOG Performance score (0, 1, 2)
# Charlson_basal  - Charlson score (<3, >=3)
# RASMutation_basal - RAS mutational status (mutant, wt, NE)
# BRAFMutation_basal  - BRAF mutational status (mutant, wt, NE)
# Microsatel_basal  - Microsatellites (msi, mss, NE)
# LocationPrimaryTumor_basal - Tumor site (left, right)
# PrimarySurgery_basal - Surgery (no, yes)
# NumberOrgansAffected_basal - Number of organs affected (1,>1)
# Liver_basal - Liver metastasis (no liver mts, liver mts + elsewhere, only liver mts 1-3 lesions <=5cm, only liver mts 4-9 lesions <=5cm, only liver mts >10 lesions or >5cm)
# Lung_basal - Lung metastasis (no lung mts, only lung mts, lung mts + elsewhere)
# PeritonealAffected_basal - Peritoneal affected (no, yes)
# NodeAffected_basal - Nodes affected (no, yes)
# LDHnormal_basal - LDH levels (normal, abnormal, NE)

## Time-varying covariates
# PS_tv  - Performance score (0,1,2,3,4)
# Charlson_tv - Charlson score (<3, >=3)
# ToxicidadGrado_tv  - Degree of toxicity (no, yes)
# TimesinceLDH_tv  - Time since last LDH measurement (normal, abnormal, NE)
# LDHnormal_tv - LDH levels (normal, abnormal, NE)


# -------------------------------Cloning ----------------------------------
#Create 2 sets of clones, assign 1 set to strategy 1 and the other to strategy 2
ds_MAB1<-ds  # We create a first set of clones of the dataset
ds_MAB1$arm<-"MAB1" #all the clones are assigned to the MAB1 arm

ds_MAB2<-ds  # We create a second set of clones of the dataset
ds_MAB2$arm<-"MAB2" #all the clones are assigned to the MAB2 arm

# -------------------------------Censoring ----------------------------------

#For strategy 1
ds_MAB1$outcome<- ifelse(adherence1, #if they adhere
                         ds_MAB1$LastStatus, #their last status will be known
                         "censored") #if they don't, their status will be censored
ds_MAB1$fup<- ifelse(adherence1, #if they adhere
      (coalesce(DeathDate, LastFupDate)) - IndexDate +1, #follow-up will be til death or last follow-up date
      pmin(DeathDate, DateOfCensoring, na.rm=TRUE) - IndexDate+1) #if they don't, follow-up will be til death or censoring date (whichever first)

#For strategy 2
ds_MAB2$outcome<- ifelse(adherence2, #if they adhere
                         ds_MAB2$LastStatus, #their last status will be known
                         "censored") #if they don't, their status will be censored
ds_MAB2$fup<- ifelse(adherence2,  #if they adhere
                     (coalesce(DeathDate, LastFupDate)) - IndexDate+1, #follow-up will be til death or last follow-up date
                     pmin(DeathDate, DateOfCensoring, na.rm=TRUE) - IndexDate+1) #if they don't, follow-up will be til death or censoring date (whichever first)


#Cloned, censored population:
dsc<-rbind(ds_MAB1, ds_MAB2)

dsc$death<-ifelse(dsc$outcome=="Dead", 1, 0) #create binary variable for individuals who suffer outcome event

# -------------Estimate subject-specific time-varying nonstabilized inverse-probabilty weights----------------

# Turn data into long format
ds$maxfup<-((coalesce(ds$DeathDate, ds$LastFupDate)-ds$IndexDate)+1)/7 #in weeks, length of maximum follow-up (uncensored)
dslong <- expandRows(ds, "maxfup", drop=F) #turn data into long format

# The denominator of the weights for each week is the probability of initiating mAB given the subjects' covariate history
# To estimate the denominator, we fit two separate models to allow the probabilities to differ according to line of treatment. 

## Model 1 (estimates p1, probability for first line treatment) ----

# Subset data to include only person-weeks when mAB is initiated as first line
dslong15w<-dslong[dslong$time<16,]
# Subset data to only include person-weeks that can be censored (weeks prior 
# to a contraindication (CI) developing or individual initiating mAB) 
dslong15w <- dslong15w[(is.na(dslong15w$fupCI) | dslong15w$time <= dslong15w$fupCI) & 
                            (is.na(dslong15w$fupMAB1) | dslong15w$time <= dslong15w$fupMAB1),]

# Define weight model 1 formula
w.model1 <- MAB1_tv ~ 
  rcs(time, knots=2) +
  age +
  agesq +
  PS_basal +
  Charlson_basal +  
  RASMutation_basal +
  BRAFMutation_basal +   
  Microsatel_basal +
  LocationPrimaryTumor_basal +
  PrimarySurgery_basal +
  NumberOrgansAffected_basal +
  Liver_basal +
  Lung_basal +
  PeritonealAffected_basal +
  NodeAffected_basal +
  LDHnormal_basal +
  PS_tv +
  Charlson_tv +
  ToxicidadGrado_tv +
  TimesinceLDH_tv +
  LDHnormal_tv

# Fit model
fit1 <- glm(w.model1, data = dslong15w, family = binomial()) 

# Compute probabilities of receiving mab at first line treatment
p.cens1.obs <- predict(fit1, type = "response") 
dslong15w$prcens1<-p.cens1.obs

## Model 2 (estimates p2, probability for second line treatment) ----

# Subset data to include only person-weeks when mAB is initiated for second line treatment
dslongQMT2w <- subset(dslong, time == fupQMT2 & !is.na(fupQMT2))

# Define weight model 2 formula
w.model2 <- MAB2_tv ~  
  rcs(time, knots=2) +
  age +
  agesq +
  PS_basal +
  Charlson_basal +
  RASMutation_basal +
  BRAFMutation_basal +
  Microsatel_basal+
  LocationPrimaryTumor_basal +
  PrimarySurgery_basal +
  NumberOrgansAffected_basal +
  Liver_basal +
  Lung_basal +
  PeritonealAffected_basal +
  NodeAffected_basal +
  LDHnormal_basal +
  PS_tv +
  Charlson_tv +
  ToxicidadGrado_tv +
  TimesinceLDH_tv +
  LDHnormal_tv

# Fit model
fit2 <- glm(w.model2 , data = dslongQMT2w, family = binomial()) 

# Compute probabilities of receiving mab at second line treatment
p.cens2.obs <- predict(fit2, type = "response")
dslongQMT2w$prcens2<-p.cens2.obs 

## Assign probabilities to each clone-week ----

#Turn cloned dataset (dsc) into long format
dsc$fup_w<-ceiling(as.numeric(dsc$fup)/7) # follow-up in weeks
dsclong <- expandRows(dsc, "fup_w", drop=F) #expand data into long format
dsclong$time <- sequence(rle(dsclong$subject_id)$lengths) #create weekly indicator

# Merge probabilities estimated in uncloned population (p1 and p2) to cloned dataset
subset_dslong15w <- dslong15w[c("subject_id", "time", "prcens1")]
subset_dslongQMT2w<- dslongQMT2w[c("subject_id", "time", "prcens2")]
merged_data <- merge(dsclong, subset_dslong15w, by = c("subject_id", "time"), all.x = TRUE)
merged_data <- merge(merged_data, subset_dslongQMT2w, by = c("subject_id", "time"), all.x = TRUE)
merged_data <- merged_data[, c(setdiff(names(merged_data), "time"), "time")]
dsclong<-merged_data

dsclong$pmAB<-ifelse(!is.na(dsclong$prcens2), dsclong$prcens2, 
		ifelse(!is.na(dsclong$prcens1), dsclong$prcens1, 0))

## Calculate weights for mAb as first line ----
#  Assuming MAB initiation at the end of the grace period if not initiated earlier

# Calculate factors for the weight at each week
dsclong$factor.w1 <-NA 
dsclong <- dsclong %>%
  mutate(
    factor.w1= case_when(
		adh_indicator == 2 ~ 1, # if adherence has already been established in previous weeks (adh2), factor=1
		time<8 ~ 1, #if adherence not yet established, but still in grace period, factor=1
		time==8 & adh_indicator == 1 & (is.na(fupMAB1) | fupMAB1>timetoadherence) ~ 1, #  if adherence is established that week but not due to initiating MAB, factor=1
		time==8 & adh_indicator == 1 & timetoadherence==fupMAB1 ~ 1/pmAB, # if adherence is established that week due to initiating MAB, factor=1/p
		time>=8 & adh_indicator == 0 ~ 0,  # if adherence not yet established, factor = 0
		TRUE ~ NA))

# Calculate weight, the cumulative product of the factors for the weight at each person-week
dsclong <- dsclong %>%
  group_by(subject_id,arm) %>%
  mutate(cw1 = cumprod(factor.w1))

## Calculate weights for mAb as second line ----
# Assuming MAB initiation at the end of the grace period if not initiated earlier

# Calculate factors for the weight at each week
dsclong$factor.w2<-NA
dsclong <- dsclong %>%
  mutate(
    factor.w2 = case_when(
		adh_indicator == 2 ~ 1, #if adherence has already been established in previous weeks (adh2), factor=1
        adh_indicator == 1 & (is.na(fupMAB2) | fupMAB2>timetoadherence) ~ 1, #the week adherence is established but not due to initiating mab, factor=1
        adh_indicator == 1 & timetoadherence==fupMAB2 ~ 1/pmAB, #the week adherence is established due to initiating mab2, factor=1/p
        adh_indicator==0 & !is.na(fupMAB1) & time==fupMAB1 ~ 0, #the week adherence is broken due to initiation of mAb as first line, factor=0
        adh_indicator==0 & !is.na(fupQMT2) & time==fupQMT2 & is.na(timetoadherence) ~ 1, #the week of qmt2 if adherence is not established, factor=1
        adh_indicator==0 ~ 1/(1-pmAB), #the weeks before starting qmt2
		TRUE ~ NA)) 

dsclong <- dsclong %>%
  group_by(subject_id,arm) %>%
  mutate(cw2 = cumprod(factor.w2))

## Merge censor weights and truncate at p99 ----

dsclong$factor.w<-coalesce(dsclong$factor.w1, dsclong$factor.w2)
dsclong$cw<-coalesce(dsclong$cw1, dsclong$cw2)

# Truncated (99th percentile) 
p99<-quantile(dsclong$cw, probs=0.99)
dsclong$cw.t99<-ifelse(dsclong$cw>p99, p99, dsclong$cw)

# ---------------- Estimate non-parametric survival curves----------------------

# Estimate (non-parametrically) Kaplan-Meier Survival Curves 
# after applying artificial censoring, without any adjustment

## Data pre-processing for survival model ----
dsc$fup_m<-ceiling(as.numeric((dsc$fup+1)/30.4375)) #follow-up time in months
surv_obj <- Surv(dsc$fup_m, dsc$death)

## Fit model ----
km_fit <- survfit(surv_obj ~ arm, data=dsc)

## Plot survival curves for each arm ----
km_plot<-ggsurvplot(km_fit, data=dsc,
                    conf.int=F, 
                    risk.table = T,
                    censor=F,
                    title = "Kaplan-Meier Survival Curves by arm", 
                    xlab = "Months", 
                    ylab = "Survival",
                    xlim = c(0,48),
                    surv.scale="percent",
                    break.x.by = 6,
                    legend.title="Strategy", 
                    legend= "bottom", 
                    legend.labs=c("mAB as first line", "mAB as second line")) 

# ----------- Estimate parametric, unadjusted Survival Curves -----------------

# Estimate (parametrically) survival Curves 
# after applying artificial censoring, without any adjustment

## Data pre-processing for survival model ----
dsc$fup_w<- ceiling(as.numeric((dsc$fup+1)/7)) #follow-up time in weeks
dsc$arm<- as.numeric(recode(dsc$arm, "MAB1" = 0, "MAB2" = 1)) #recode arm variable to integer to allow for interaction term
dsc.surv<-expandRows(dsc, "fup_w", drop=F) #Turn data into long format 
dsc.surv$time <- sequence(rle(dsc.surv$subject_id)$lengths) #create a variable to identify each time unit (in weeks)
dsc.surv$event <- ifelse(dsc.surv$time==dsc.surv$fup_w & dsc.surv$death==1, 1, 0) #create a variable that indicates the week of the outcome event (death)

## Fit model ----
# Fit parametric pooled logistic regression hazards model (modelling the probability of the event not occurring)
unadj.glm.I <- glm(event==0 ~ arm 
                   + rcs(time, knots=2) + arm:rcs(time, knots=2)
                   , family=binomial(), data=dsc.surv)

## Estimate risk of survival for each person-week ----
arm0 <- data.frame(cbind(seq(1, 209),0))
arm1 <- data.frame(cbind(seq(1, 209),1))

colnames(arm0) <- c("time", "arm")
colnames(arm1) <- c("time", "arm")

arm0$p.noevent0 <- predict(unadj.glm.I, arm0, type="response") #predict probability of no-event at each person-week
arm1$p.noevent1 <- predict(unadj.glm.I, arm1, type="response")

arm0$surv0 <- cumprod(arm0$p.noevent0) #computation of cumulative probability of survival for each person-week
arm1$surv1 <- cumprod(arm1$p.noevent1)

arm0$risk0<-1-arm0$surv0 #computation of risk of death for each week
arm1$risk1<-1-arm1$surv1 

unadj.graph <- merge(arm0, arm1, by=c("time"))
unadj.graph$survdiff <- unadj.graph$surv1-unadj.graph$surv0 #Calculate survival risk difference at each week
unadj.graph$riskratio<- unadj.graph$risk1/unadj.graph$risk0 #Calculate risk ratio at each week
unadj.graph$time_mo <- unadj.graph$time / 4.3452  # Time in months

## Plot survival curves for each arm ----
unadj.plot<-ggplot(unadj.graph, aes(x=time_mo, y=surv)) + 
  geom_line(aes(y = surv0, colour = "mAB 1st line")) + 
  geom_line(aes(y = surv1, colour = "mAB 2nd line")) + 
  xlab("Months") + 
  scale_x_continuous(limits = c(0, 48), breaks=seq(0,48,6)) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  ylab("Survival") + 
  ggtitle("Parametric Survival Curve, no adjustment") + 
  labs(colour="Strategy") +
  theme_bw() + 
  theme(legend.position="bottom")

## Print effect estimates at 48-months of follow-up ----
print(unadj.graph[unadj.graph$time == 209, c("surv0", "surv1")]) #survival for each arm at 48 months
print(unadj.graph[unadj.graph$time == 209, c("survdiff")]) #Survival risk difference at 48 months
print(unadj.graph[unadj.graph$time == 209, c("riskratio")]) #Survival Risk ratio at 48 months


# --------- Estimate parametric, baseline-adjusted Survival Curves ------------

# Estimate (parametrically) Survival Curves 
# after applying artificial censoring, with baseline adjustment

## Data pre-processing for survival model ----
dsc$fup_w<- ceiling(as.numeric((dsc$fup+1)/7)) #follow-up time in weeks
dsc$arm <- as.numeric(recode(dsc$arm, "MAB1" = 0, "MAB2" = 1)) #change arm variable to integer to allow for interaction term
dsc.surv <- expandRows(dsc, "fup_w", drop=F) #Turn data into long format (person-week data)
dsc.surv$time <- sequence(rle(dsc.surv$subject_id)$lengths) #create a variable to identify each time unit (in weeks)
dsc.surv$event <- ifelse(dsc.surv$time==dsc.surv$fup_w & dsc.surv$death==1, 1, 0) #create a variable that indicates the week of the event (death)

## Fit model ----
# Fit of parametric pooled logistic regression hazards model 
# with baseline covariates (modelling the probability of the event not occurring)
bsl.adj.glm.I <- glm(event==0 ~ arm 
                     + rcs(time, knots=2) 
                     + arm:rcs(time, knots=2)
                     + age 
                     + agesq
                     + PS_basal 
                     + Charlson_basal 
                     + RASMutation_basal  
                     + BRAFMutation_basal 
                     + Microsatel_basal 
                     + LocationPrimaryTumor_basal  
                     + PrimarySurgery_basal 
                     + NumberOrgansAffected_basal 
                     + Liver_basal 
                     + Lung_basal 
                     + PeritonealAffected_basal 
                     + NodeAffected_basal 
                     + LDHnormal_basal 
                     , data=dsc.surv, family=binomial())


## Estimate risk of survival for each person-week ----
bsl.arm0 <- expandRows(dsc, count=209, count.is.col=F)
bsl.arm0$time <- rep(seq(1, 209), nrow(dsc))
bsl.arm0<- bsl.arm0 %>% group_by(arm) %>% mutate(clone = arm) 
bsl.arm0$arm <- 0

bsl.arm1 <- bsl.arm0
bsl.arm1$arm <- 1

bsl.arm0$p.noevent0 <- predict(bsl.adj.glm.I, bsl.arm0, type="response") #predict pr of no-event at each person-week
bsl.arm1$p.noevent1 <- predict(bsl.adj.glm.I, bsl.arm1, type="response")

bsl.arm0.surv <- bsl.arm0 %>% group_by(subject_id, clone) %>% mutate(surv0 = cumprod(p.noevent0)) #computation of cumulative probability of survival for each person-week
bsl.arm1.surv <- bsl.arm1 %>% group_by(subject_id, clone) %>% mutate(surv1 = cumprod(p.noevent1))

bsl.surv0 <- aggregate(bsl.arm0.surv, by=list(bsl.arm0.surv$time), FUN=mean)[c("arm", "time", "surv0")] #computation of mean probability of survival at each week across all individuals 
bsl.surv1 <- aggregate(bsl.arm1.surv, by=list(bsl.arm1.surv$time), FUN=mean)[c("arm", "time", "surv1")]

bsl.surv0$risk0<-1-bsl.surv0$surv0 #computation of risk of death for each week
bsl.surv1$risk1<-1-bsl.surv1$surv1

bsl.graph <- merge(bsl.surv0, bsl.surv1, by=c("time"))
bsl.graph$survdiff <- bsl.graph$surv1-bsl.graph$surv0
bsl.graph$riskratio<- bsl.graph$risk1/bsl.graph$risk0
bsl.graph$time_mo <- bsl.graph$time / 4.3452  # Time in months

## Plot survival curves for each arm ----
bsl.plot<-ggplot(bsl.graph, aes(x=time_mo, y=surv)) + 
  geom_line(aes(y = surv0, colour = "mAB 1st line")) + 
  geom_line(aes(y = surv1, colour = "mAB 2nd line")) + 
  xlab("Months") + 
  scale_x_continuous(limits = c(0, 48), breaks=seq(0,48,6)) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  ylab("Survival") + 
  ggtitle("Parametric Survival Curve, baseline adjusted") + 
  labs(colour="Strategy") +
  theme_bw() + 
  theme(legend.position="bottom")

## Print effect estimates at 48-months of follow-up ----
print(bsl.graph[bsl.graph$time == 209, c("surv0", "surv1")]) #survival for MAB1 and MAB2 at 48 months
print(bsl.graph[bsl.graph$time == 209, c("survdiff")]) #Survival risk difference at 48 months
print(bsl.graph[bsl.graph$time == 209, c("riskratio")]) #Survival Risk ratio at 48 months

# ------ Estimate parametric, baseline and weight-adjusted Survival Curves -----

# Estimate (parametrically) survival curves 
# after applying artificial censoring, with baseline and weights adjustment

## Data pre-processing for survival model ----
dsclong$event <- ifelse(dsclong$time==dsclong$fup_w & dsclong$death==1, 1, 0)  #create a variable that indicates the week of the event (death)
dsclong$arm <- as.numeric(recode(dsclong$arm, "MAB1" = 0, "MAB2" = 1)) #change arm to integer to allow for interaction term
dsclong$time_m<- dsclong$time/4.3452 #time in months

## Fit model ----
# Fit of parametric pooled logistic regression hazards model with baseline covariates and censoring weights (modelling the probability of the event not occurring)
full.adj.glm.I<-glm(event==0 ~
                    arm +
                    rcs(time, knots=2) +
                    arm:rcs(time, knots=2) +
                    age +
                    agesq +
                    PS_basal +
                    Charlson_basal +
                    RASMutation_basal +
                    BRAFMutation_basal +
                    Microsatel_basal +
                    LocationPrimaryTumor_basal +
                    PrimarySurgery_basal +
                    NumberOrgansAffected_basal +
                    Liver_basal +
                    Lung_basal +
                    PeritonealAffected_basal +
                    NodeAffected_basal +
                    LDHnormal_basal,
                    data = dsclong, family = binomial(), weights=cw.t99)

## Estimate risk of survival for each person-week ----
# We need to predict survival for each subject under arm0 and arm1
bsl.cw.arm0 <- expandRows(dsc, count=209, count.is.col=F) 
bsl.cw.arm0$time <- rep(seq(1, 209), nrow(dsc))
bsl.cw.arm0<- bsl.cw.arm0 %>% group_by(arm) %>% mutate(clone = arm) 
bsl.cw.arm0$arm <- 0

bsl.cw.arm1 <- bsl.cw.arm0
bsl.cw.arm1$arm <- 1

bsl.cw.arm0$p.noevent0 <- predict(full.adj.glm.I, bsl.cw.arm0, type="response") #predict pr of survival at each person-week
bsl.cw.arm1$p.noevent1 <- predict(full.adj.glm.I, bsl.cw.arm1, type="response")

bsl.cw.arm0.surv <- bsl.cw.arm0 %>% group_by(subject_id, clone) %>% mutate(surv0 = cumprod(p.noevent0)) #computation of cumulative probability of survival for each person through time t
bsl.cw.arm1.surv <- bsl.cw.arm1 %>% group_by(subject_id, clone) %>% mutate(surv1 = cumprod(p.noevent1))

bsl.cw.surv0 <- aggregate(bsl.cw.arm0.surv, by=list(bsl.cw.arm0.surv$time), FUN=mean)[c("arm", "time", "surv0")] #computation of mean probability of survival at each week across all individuals   
bsl.cw.surv1 <- aggregate(bsl.cw.arm1.surv, by=list(bsl.cw.arm1.surv$time), FUN=mean)[c("arm", "time", "surv1")]

bsl.cw.surv0$risk0<-1-bsl.cw.surv0$surv0 #computation of risk of death for each week
bsl.cw.surv1$risk1<-1-bsl.cw.surv1$surv1

bsl.cw.graph <- merge(bsl.cw.surv0, bsl.cw.surv1, by=c("time"))
bsl.cw.graph$survdiff <- bsl.cw.graph$surv1-bsl.cw.graph$surv0
bsl.cw.graph$riskratio<- bsl.cw.graph$risk1/bsl.cw.graph$risk0
bsl.cw.graph$time_mo <- bsl.cw.graph$time / 4.3452  # Time in months

## Plot survival curves for each arm ----
bsl.cw.plot<-ggplot(bsl.cw.graph, aes(x=time_mo, y=surv)) + 
  geom_line(aes(y = surv0, colour = "mAB 1st line")) + 
  geom_line(aes(y = surv1, colour = "mAB 2nd line")) + 
  xlab("Months") + 
  scale_x_continuous(limits = c(0, 48), breaks=seq(0,48,6)) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  ylab("Survival") + 
  ggtitle("Parametric Survival Curve, baseline and weights- adjusted") + 
  labs(colour="Strategy") +
  theme_bw() + 
  theme(legend.position="bottom")

## Print effect estimates at 48-months of follow-up ----
print(bsl.cw.graph[bsl.cw.graph$time == 209, c("surv0", "surv1")]) #survival for MAB1 and MAB2 at 48 months
print(bsl.cw.graph[bsl.cw.graph$time == 209, c("survdiff")]) #Survival risk difference at 48 months
print(bsl.cw.graph[bsl.cw.graph$time == 209, c("riskratio")]) #Survival Risk ratio at 48 months

# ------ Bootstraps to obtain 95% CI ----

## Define bootstrap function for the baseline and post-baseline effect estimates ----

btsp_input_ds<-ds

# Create a function that conducts the analysis for the subsetted bootstrap sample,
# and stores the effect estimates at each iteration
boot.func <- function(data, indices, file_name="ccw_bootstraps_full.adj.xlsx") { 
  ds<-data[indices,]
  
  #Create a bootstrap-specific subject custom code:
  ds$subject_id_bts<-ave(ds$subject_id, ds$subject_id, FUN = function(i) paste0(i, '_', seq_along(i)))
  
# Include code for the following sections into the bootstrap function, substituting "subject_id" by "subject_id_bts"
  # Cloning
  # Censoring
  # Estimating IPW
  # Fit of parametric baseline and post-baseline IPW-adjusted model

# Store effect estimates in a spreadsheet (each iteration saved in a new sheet)
 sheet_name_prefix <- "Bootstrap_"  # Create a new workbook and add the first sheet
  if (!file.exists(file_name)) {
  wb <- createWorkbook()
  sheet_name <- paste0(sheet_name_prefix, "1")
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, bsl.cw.graph)
  saveWorkbook(wb, file_name, overwrite = TRUE)
} else {
  # Append to existing workbook
  wb <- loadWorkbook(file_name)
  existing_sheets <- sheets(wb)
  # Find the next available sheet name
  next_sheet_number <- length(existing_sheets) + 1
  sheet_name <- paste0(sheet_name_prefix, next_sheet_number)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, bsl.cw.graph)
  saveWorkbook(wb, file_name, overwrite = TRUE)
 }
}


## Generate bootstrap samples ----

for (i in 1:500) {  
  set.seed(100 + i) #set seed for each iteration
  indices <- sample(nrow(btsp_input_ds), replace = TRUE) #sample subjects for each iteration
  boot.func(data = btsp_input_ds, indices = indices, file_name = "ccw_bootstraps_full.adj.xlsx") #run boot function
}


## Calculate CIs ----

#Retrieve effect estimates at each iteration stored in the xlsx
file_path <- "/ccw_bootstraps_full.adj.xlsx"

sheet_names <- excel_sheets(file_path)

all_data <- lapply(seq_along(sheet_names), function(i) {
  data <- read_excel(file_path, sheet = sheet_names[i])
  data <- data %>% mutate(NBootstrap = i)
  return(data)
})

final_data <- bind_rows(all_data)

#Calculate CIs

calculate_ci <- function(values) {
  round(quantile(values, probs = c(0.025, 0.975), na.rm = TRUE), 3)
}
ci_results <- final_data %>%
  group_by(time) %>%
  summarise(
    surv095CILow = calculate_ci(surv0)[1],
    surv095CIUp = calculate_ci(surv0)[2],
    surv195CILow = calculate_ci(surv1)[1],
    surv195CIUp = calculate_ci(surv1)[2],
    survdiff95CILow = calculate_ci(survdiff)[1],
    survdiff95CIUp = calculate_ci(survdiff)[2],
    RiskRatio95CILow = calculate_ci(riskratio)[1],
    RiskRatio95CIUp = calculate_ci(riskratio)[2]
  ) %>%
  ungroup()

bsl.cw.graph.plusCI <- bsl.cw.graph %>% #Append the effect estimates with the CIs
  left_join(ci_results, by = "time")

# Save effect CIs at specific time points

time_points <- c(27, 53, 79, 105, 131, 157, 183, 209) #define time points at which we want to extract estimates (in weeks)

bsl.cw.adj_timepoints <- bsl.cw.graph.plusCI %>%
  filter(time %in% time_points) %>% 
  mutate(
    RR_CI = paste0(round(riskratio,3)," (", RiskRatio95CILow, ";", RiskRatio95CIUp, ")"),
    RD_CI = paste0(round(survdiff, 3), " (", survdiff95CILow, ";", survdiff95CIUp, ")"),
    R0_CI = paste0(round(surv0,3), " (", surv095CILow, ";", surv095CIUp, ")"),
    R1_CI = paste0(round(surv1,3), " (", surv195CILow, ";", surv195CIUp, ")")
  ) %>%
  select(time, RR_CI, RD_CI, R0_CI, R1_CI) 

write.csv(bsl.cw.adj_timepoints, file = "/ccw_full.adjusted_cis.csv", row.names = FALSE)
  

