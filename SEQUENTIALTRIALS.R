###########################################################################. 
# Project: 0358082 PDA GEMCAD 1401                                        #
# Title: "Observational Study to Evaluate the Use of Targeted             #
#        Therapies in Metastatic Colorectal Cancer - GEMCAD 1401"         #
# Authors: Manuel Zamparini, Julia Vila Guilera, Xabier Garcia de Albeniz #
# Analysis: Sequential Trial Emulation Analysis                           #
###########################################################################. 

# Emulate 8 sequential trials----

## Create 8 sequential trials ----

# Transform data to long format, to have a better view of the weeks
ds12l <- expandRows(ds12, "fupW", drop=F)  

# Exclude all patients not eligible at time 0 of trial 1 (exclude those with a previous progression, toxicity, ps34, or other reasons)
excluded <- ds12l %>% filter(pmin(progression, ps34, otherreasons, na.rm = TRUE) <= start) 
ds12le <- ds12l %>% anti_join(excluded, by = c("subject_custom_code"))
excluded <- c() #Create object that keeps track of excluded subjects in each sequential trial creation

trial_list <- list() # THIS WILL CONTAINS ALL THE TRIALS' DATASET
trial_data <- ds12le # STARTING TRIAL

for (i in 1:8) {
  if (i > 1) {
    trial_data <- ds12le %>%
      anti_join(tibble(subject_custom_code = excluded), by = "subject_custom_code") %>% # BY THE SECOND TRIAL, WE MUST EXCLUDE ALL UN-ELIGIBLE PATIENTS
      arrange(subject_custom_code, week) %>%
      group_by(subject_custom_code) %>%
      slice(-(1:(i-1))) %>% #IN EACH SEQUENTIAL TRIAL WE SHIFT THE STARTING WEEK
      ungroup()
  }


  trial_current <- trial_data %>%
    group_by(subject_custom_code) %>%
    mutate(trial = i,
           arm = ifelse(first(MAB1wRetained) == 1, 0, 1)) %>% #DEFINE ARM: IF MAB1 IS STARTED ON THE FIRST WEEK: ARM=0 (MAB IN FIRST LINE GROUP), ELSE ARM=1 (MAB IN SECOND LINE GROUP)
    mutate(keep_row = ifelse(arm == 0, TRUE, MAB1wRetained == 0)) %>% # IF ARM==0, WE KEEP ALL LINES, ELSE WE KEEP JUST LINES BEFORE STARTING MAB1
    mutate(keep_row = ifelse(arm == 1 & MAB1wRetained == 0 & !is.na(QMT2) & is.na(MAB2), QMT2wRetained == 0, keep_row)) %>% #IF ARM==1 (AND NOT YET STARTED MAB1), WE CENSOR IN THE MOMENT OF QMT2 WITHOUT MAB2
    mutate(keep_row2 = lag(keep_row, default = TRUE)) %>% #THE LAG FUNCTION KEEPS THE FIRST WEEK OF CENSORING 
    filter(keep_row2) %>%
    select(-keep_row, -keep_row2) %>%
    ungroup() %>%
    group_by(subject_custom_code) %>%
    slice_head(n = 209) %>% # MAXIMUM LENGHT OF FOLLOW-UP=209 WEEKS (SHIFTING BETWEEN TRIALS)
    ungroup()
  
    trial_list[[i]] <- trial_current # save the df
  
  
  excluded <- union( # update excluded list
    excluded,
    trial_current %>%
      group_by(subject_custom_code) %>%
      summarise(
        exclude = any((arm == 0 | (!is.na(progression) & progression <= startW & week==i+1) |  #SUBJECTS THAT STARTED MAB1 OR DEVELOPED CONTRO-INDICATIONS BEFORE WEEK X ARE UNELIGIBLE FOR TRIAL X
                         (!is.na(ps34) & ps34 <= startW & week==i+1) | 
                         (!is.na(otherreasons) & otherreasons >= startW-6 & otherreasons <= startW & week==i+1))),
        .groups = 'drop'
      ) %>%
      filter(exclude) %>%
      pull(subject_custom_code)
  )
}


### Introduce baseline variations into the 8 sequential trials ----

# IN EACH SEQUENTIAL TRIALS, THE WEEK CONSIDERED AS BASELINE IS SHIFTED, SO WE NEED TO UPDATE BASELINE VALUES

# Define a function to update baseline values based on the evaluation date and start date


update_basal <- function(data, basal_var, eval_prefix, date_prefix, date_comparator) {
  # Loop through each of the 23 potential evaluation points
  for (i in 1:23) {
    eval_var <- paste0(eval_prefix, i)  # Create variable name for the evaluation.
    date_var <- paste0(date_prefix, i)  # Create variable name for the date.
    # Update the basal variable if the date of evaluation is on or before the comparator date.
    data <- data %>%
      mutate(!!basal_var := case_when(
        get(date_var) <= get(date_comparator)  & !is.na(get(eval_var)) ~ get(eval_var),  # If condition is met, update basal value.
        TRUE ~ .data[[basal_var]]  # Otherwise, keep the original basal value
      ))
  }
  return(data)
}

# Define a function to update trial data with new baseline values.
update_trial <- function(trial_data, trial_number) {
  # Get the first entry for each subject to establish a unique dataset.
  trial_unique <- trial_data %>%
    group_by(subject_custom_code) %>%
    slice(1) %>%
    ungroup()
  # Update basal values for different variables using the update_basal function.
  trial_unique <- update_basal(trial_unique, "PS_basal_imp", "PS_eval_", "FechaEval_", "startW")
  trial_unique <- update_basal(trial_unique, "Charlson2_basal", "Charlson2_eval_", "FechaEval_", "startW")
  trial_unique <- update_basal(trial_unique, "LDH_basal", "LDH_eval_", "FechaEval_", "startW")
#......
  # Merge the updated unique baseline data back into the original trial data.
  trial_updated <- trial_data %>%
    left_join(trial_unique, by = c("subject_custom_code" = "subject_custom_code")) %>%
    mutate(
      PS_basal_impNEW = PS_basal_imp.y,
      Charlson2_basalNEW = Charlson2_basal.y,
      LDHnormal_basalNEW = LDHnormal_basal
      #....
    ) %>%
    select(-matches("\\.y$"))  # Remove duplicated columns from the join.
  # Rename variables to remove suffixes and clean up the names.
  trial_updated <- trial_updated %>%
    rename_with(~ gsub("\\.x", "", .x), ends_with(".x"))
  
  
  trial_updated <- trial_updated %>%
    select(-LDHnormal_basal)
  
  
  # Assign the updated data back to the global environment with a specific trial name.
  assign(paste0("trial", trial_number), trial_updated, envir = .GlobalEnv)
}

# Run the update process for trials 1 through 8.
for (i in 1:8) {
  trial_data <- get(paste0("trial", i))  # Retrieve each trial's data.
  update_trial(trial_data, i)  # Apply updates using the defined function.
}

## Plot time under follow-up ----

for (i in 1:8) {
  ## TRANSFORM EACH SEQUENTIAL TRIALS LONG DF INTO WIDE FORMAT
  trial_data_wide  <- get(paste("trial", i, sep="")) %>%
    group_by(subject_custom_code, arm) %>%                 
    summarise(
      startW_min = min(startW),                  
      endW_max = max(endW),
      MAB1 = min(MAB1),
      MAB2 = min(MAB2),
      QMT2 = min(QMT2),
      death = max(death),
      death_tot = max(death_tot), #Variables _tot are immune from follow-up truncation, useful to describe the total observed period
      fup_w = ceiling(as.numeric(difftime(endW_max, startW_min, units = "weeks"))),
      fup_w_tot = max(fupW),
      fup_m = ceiling(fup_w / 4.345), # Convert weeks to months
      fup_m_tot = ceiling(fup_w_tot / 4.345)
    ) %>%
    ungroup() 
  
    trial_data_wide<- trial_data_wide %>% 
    mutate(
      fup = endW_max-startW_min,
      censored4MAB1 = ifelse(arm==1 & !is.na(MAB1) & MAB1<=endW_max+7,1,0), # MAB IN SECOND LINE PATIENTS ARE CENSORED IF THEY START MAB1
      censored4noMAB2 = ifelse(arm==1 & is.na(MAB2) & !is.na(QMT2) & QMT2<=endW_max+7,1,0) #MAB IN SECOND LINE PATIENTS ARE CENSORED IF THEY START QMT2 BUT NOT MAB2 
    )
  
  # Figure 1. Time under followup after applying the artificial censoring
  
  #Plot the Time under follow up for clones
  surv_obj <- Surv(trial_data_wide$fup_m) #time in months
  km_fit <- survfit(surv_obj ~ arm, data=trial_data_wide)
  tfup_plot<-ggsurvplot(km_fit, data=trial_data_wide,
                        conf.int=TRUE, 
                        risk.table = TRUE,
                        censor=FALSE,
                        title = paste("Trial ",i,sep=""), 
                        xlab = "Months", 
                        ylab = "Survival",
                        xlim = c(0,48),
                        break.x.by = 6,
                        legend.title="Exposure", 
                        legend= "bottom", 
                        legend.labs=c("mAB 1st line", "mAB 2nd line"))
 
  
  ## Estimate non-parametric survival curves, unadjusted (Kaplan Meier) ----
  
  # Figure 2. Kaplan Meier figure of overall survival after applying the artificial censoring, without any adjustment
  
  surv_obj <- Surv(trial_data_wide$fup_m, trial_data_wide$death) 
  km_fit <- survfit(surv_obj ~ arm, data=trial_data_wide)
  summary(km_fit)
  km_plot<-ggsurvplot(km_fit, data=trial_data_wide,
                      conf.int=F, 
                      risk.table = T,
                      censor=F,
                      title = paste("Trial ",i,sep=""), 
                      xlab = "Months", 
                      ylab = "Survival",
                      xlim = c(0,48),
                      surv.scale="percent",
                      break.x.by = 6,
                      legend.title="Strategy", 
                      legend= "bottom", 
                      legend.labs=c("KM MAB1", "KM MAB2")) 
  
  ## Estimate parametric survival curves, unadjusted ----
  
  # Figure 3. Parametric survival curve (pooled log reg) of overall survival after applying the artificial censoring, without any adjustment
  
  ## fit of parametric plr hazards model (modelling the probability of the event occurring)
  unadj.glm.I <- glm(event==0 ~ arm 
                     + rcs(time, knots=2) + arm:rcs(time, knots=2)
                     , family=binomial(), data=trial_data_long)
  
  # creation of dataset with all time points under each treatment level
  arm0 <- data.frame(cbind(seq(0, 209),0,(seq(0, 209))^2))
  arm1 <- data.frame(cbind(seq(0, 209),1,(seq(0, 209))^2))
  
  colnames(arm0) <- c("time", "arm", "timesq")
  colnames(arm1) <- c("time", "arm", "timesq")
  
  # assignment of estimated (1-hazard) to each person-week */
  arm0$p.noevent0 <- predict(unadj.glm.I, arm0, type="response")
  arm1$p.noevent1 <- predict(unadj.glm.I, arm1, type="response")
  
  # computation of survival for each person-week
  arm0$surv0 <- cumprod(arm0$p.noevent0) #cumulative sequential product of values in p.noevent0
  arm1$surv1 <- cumprod(arm1$p.noevent1)
  
  #computation of risk for each person-week
  arm0$risk0<-1-arm0$surv0
  arm1$risk1<-1-arm1$surv1
  
  # some data management to plot estimated survival curves
  unadj.graph <- merge(arm0, arm1, by=c("time", "timesq"))
  unadj.graph$survdiff <- unadj.graph$surv1-unadj.graph$surv0
  unadj.graph$riskratio<- unadj.graph$risk1/unadj.graph$risk0
  unadj.graph$time_mo <- unadj.graph$time / 4.3452  # Time in months
  
  # plot
  unadj.plot<-ggplot(unadj.graph, aes(x=time_mo, y=surv)) + 
    geom_line(aes(y = surv0, colour = "mAB 1st line")) + 
    geom_line(aes(y = surv1, colour = "mAB 2nd line")) + 
    xlab("Months") + 
    scale_x_continuous(limits = c(0, 48), breaks=seq(0,48,6)) +
    scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
    ylab("Survival") + 
    ggtitle(paste("Trial ",i,sep="")) + 
    labs(colour="Strategy") +
    theme_bw() + 
    theme(legend.position="bottom")

  ### Print effect estimates: Unadjusted absolute and relative risks ----
  
  print(i)
  # 48 month survival for MAB1 and MAB2
  print(unadj.graph[unadj.graph$time == 209, c("surv0", "surv1")])
  # Survival difference at 48 months
  print(unadj.graph[unadj.graph$time == 209, c("survdiff")])
  #Risk ratio at 48 months
  print(unadj.graph[unadj.graph$time == 209, c("riskratio")])
  
  ### Combined parametric and non parametric plot for comparison----
  
  combined_plot<- km_plot$plot +
    geom_line(aes(x=time_mo, y = surv0, colour = "PLR MAB1"), data=unadj.graph) + 
    geom_line(aes(x=time_mo, y = surv1, colour = "PLR MAB2"), data=unadj.graph) + 
    xlab("Months") + 
    scale_x_continuous(limits = c(0, 48), breaks=seq(0,48,6)) +
    ylab("Survival") + 
    ggtitle(paste("Trial ",i,sep="")) + 
    labs(colour="Strategy") +
    theme_minimal() + 
    theme(plot.background = element_rect(fill = "white", color = NA),)+
    theme(legend.position="bottom")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    scale_color_manual(values = c("KM MAB1"= "#ff9999", "PLR MAB1" = "#ff9999", "KM MAB2"="#33CCCC","PLR MAB2" = "#33CCCC"))

  ### Estimate unadjusted Hazard ratios (from GLM model without interaction)----

  #To estimate the Hazard ratios, we need to run the same GLM model but without interactions between arm and time
  unadj.glm.noI <- glm(event==0 ~ arm + rcs(time, knots=2), family=binomial(), data=trial_data_long)
  # Extract the HR and the robust standard error
  hr<-1/exp(coef(unadj.glm.noI)[["arm"]]) # to get Hazard Ratio for arm
  robust_se<-coeftest(unadj.glm.noI, vcov=vcovHC(unadj.glm.noI, type="HC1"))[2,2] # to get robust SE estimates for arm
  robust_ci_upper <- 1/exp(log(1/hr) - 1.96 * robust_se)
  robust_ci_lower <- 1/exp(log(1/hr) + 1.96 * robust_se)
  

  
  ##Estimate parametric survival curves, baseline adjusted ----
  
  # Figure 4. Parametric survival curve of overall survival after applying the artificial censoring, adjusted for baseline variables
  
  # fit of parametric hazards model with covariates
  bsl.adj.glm.I <- glm(event==0 ~ arm 
                       + rcs(time, knots=2) 
                       + arm:rcs(time, knots=2)
                       + age #age(in years)
                       + agesq
                       + PS_basal_impNEW #ECOG PS (0, 1, 2)
                       + Charlson2_basalNEW #Charlson (<3, >=3)
                       + RASMutation_basal  # RAS (mutant, wt, NE)
                       + BRAFMutation_basal #BRAF (mutant, wt, NE)
                       + Microsatel_basal #Microsatel (msi, mss, NE)
                       + LocationPrimaryTumor_basal  #Tumor site (left, right)
                       + PrimarySurgery_basal # Surgery (no, yes)
                       + NumberOrgansAffected_basal #Number of organs (1,>1)
                       + Liver2_basal #(no liver mts, liver mts + elsewhere, only liver mts 1-3 lesions <=5cm, only liver mts 4-9 lesions <=5cm, only liver mts >10 lesions or >5cm)
                       + Lung2_basal #(no lung mts, only lung mts, lung mts + elsewhere)
                       + PeritonealAffected_basal # (no, yes)
                       + NodeAffected_basal # (no, yes)
                       + LDHnormal_basalNEW #(normal, abnormal, NE)
                       , data=trial_data_long, family=binomial())
  
  dsc.c0<- trial_data_wide %>%
    filter(arm==0)
  bsl.arm0 <- expandRows(dsc.c0, count=209, count.is.col=F)
  bsl.arm0$time <- rep(seq(1, 209), nrow(dsc.c0))
  
  dsc.c1<- trial_data_wide %>%
    filter(arm==1)
  bsl.arm1 <- expandRows(dsc.c1, count=209, count.is.col=F)
  bsl.arm1$time <- rep(seq(1, 209), nrow(dsc.c1))
  
  # assignment of estimated (1-hazard) to each person-week */
  bsl.arm0$p.noevent0 <- predict(bsl.adj.glm.I, bsl.arm0, type="response") #predict pnoevent at each person-week while each person is followed
  bsl.arm1$p.noevent1 <- predict(bsl.adj.glm.I, bsl.arm1, type="response")
  
  # computation of survival for each person-week
  bsl.arm0.surv <- bsl.arm0 %>% group_by(subject_custom_code) %>% mutate(surv0 = cumprod(p.noevent0))
  bsl.arm1.surv <- bsl.arm1 %>% group_by(subject_custom_code) %>% mutate(surv1 = cumprod(p.noevent1))
  
  bsl.surv0 <- aggregate(bsl.arm0.surv, by=list(bsl.arm0.surv$time), FUN=mean)[c("arm", "time", "surv0")]
  bsl.surv1 <- aggregate(bsl.arm1.surv, by=list(bsl.arm1.surv$time), FUN=mean)[c("arm", "time", "surv1")]
  
  # computation of risk for each person-week
  bsl.surv0$risk0<-1-bsl.surv0$surv0
  bsl.surv1$risk1<-1-bsl.surv1$surv1
  # some data management to plot estimated survival curves
  bsl.graph <- merge(bsl.surv0, bsl.surv1, by=c("time"))
  bsl.graph$survdiff <- bsl.graph$surv1-bsl.graph$surv0
  bsl.graph$riskratio<- bsl.graph$risk1/bsl.graph$risk0
  bsl.graph$time_mo <- bsl.graph$time / 4.3452  # Time in months
  
  # plot
  bsl.plot<-ggplot(bsl.graph, aes(x=time_mo, y=surv)) + 
    geom_line(aes(y = surv0, colour = "mAB 1st line")) + 
    geom_line(aes(y = surv1, colour = "mAB 2nd line")) + 
    xlab("Months") + 
    scale_x_continuous(limits = c(0, 48), breaks=seq(0,48,6)) +
    scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
    ylab("Survival") + 
    ggtitle(paste("Trial ",i,sep="")) + 
    labs(colour="Strategy") +
    theme_bw() + 
    theme(legend.position="bottom")
  
  ### Print effect estimates: baseline-adjusted absolute and relative risks ----
  
  # 48 month survival for MAB1 and MAB2
  print(bsl.graph[bsl.graph$time == 209, c("surv0", "surv1")])
  # Survival difference at 48 months
  print(bsl.graph[bsl.graph$time == 209, c("survdiff")])
  #Risk ratio at 48 months
  print(bsl.graph[bsl.graph$time == 209, c("riskratio")])

  ### Estimate baseline-adj Hazard ratios from GLM model without interaction ----
  
  # fit of parametric hazards model with covariates
  bsl.adj.glm.noI <- glm(event==0 ~ arm
                         + rcs(time, knots=2) 
                         + age #age(in years)
                         + agesq
                         + PS_basal_impNEW #ECOG PS (0, 1, 2, and NAs)
                         + Charlson2_basalNEW #Charlson (<3, >=3)
                         + RASMutation_basal  # RAS (mutant, wt, NE)
                         + BRAFMutation_basal #BRAF (mutant, wt, NE)
                         + Microsatel_basal #Microsatel (msi, mss, NE)
                         + LocationPrimaryTumor_basal  #Tumor site (left, right)
                         + PrimarySurgery_basal # Surgery (no, yes)
                         + NumberOrgansAffected_basal #Number of organs (1,>1)
                         + Liver2_basal #(no liver mts, liver mts + elsewhere, only liver mts 1-3 lesions <=5cm, only liver mts 4-9 lesions <=5cm, only liver mts >10 lesions or >5cm)
                         + Lung2_basal #(no lung mts, only lung mts, lung mts + elsewhere)
                         + PeritonealAffected_basal # (no, yes)
                         + NodeAffected_basal # (no, yes)
                         + LDHnormal_basalNEW #(normal, abnormal, NE)
                         , data=trial_data_long, family=binomial())
  
  # Extract the HR and the robust standard error
  hr<-1/exp(coef(bsl.adj.glm.noI)[["arm"]]) # to get Hazard Ratio for arm
  robust_se<-coeftest(bsl.adj.glm.noI, vcov=vcovHC(bsl.adj.glm.noI, type="HC1"))[2,2]
  robust_ci_upper <- 1/exp(log(1/hr) - 1.96 * robust_se)
  robust_ci_lower <- 1/exp(log(1/hr) + 1.96 * robust_se)
  
  ## Create time varying covariates ----
  
  #PS time varying
  trial_data_long <- trial_data_long %>%
    mutate(
      PS_tv = PS_basal_impNEW,
      PS_tv = case_when(
        time >= fupvisit_1 & (is.na(fupvisit_2) | time < fupvisit_2) ~ PS_eval_1,
        #......................
        time >= fupvisit_22 & (is.na(fupvisit_23) | time < fupvisit_23) ~ PS_eval_22,
        time >= fupvisit_23 ~ PS_eval_23,
        TRUE ~ PS_tv
      )
    )
  
  #Charlson2 time varying
  trial_data_long <- trial_data_long %>%
    mutate(
      Charlson2_tv = Charlson2_basalNEW,
      Charlson2_tv = case_when(
        time >= fupvisit_1 & (is.na(fupvisit_2) | time < fupvisit_2) ~ Charlson2_eval_1,
        #......................
        time >= fupvisit_22 & (is.na(fupvisit_23) | time < fupvisit_23) ~ Charlson2_eval_22,
        time >= fupvisit_23 ~ Charlson2_eval_23,
        TRUE ~ Charlson2_tv
      )
    )
  
  
  #Toxicidad time varying
  trial_data_long <- trial_data_long %>%
    mutate(
      ToxicidadGrado_tv = ToxicidadGrado_eval_1,
      ToxicidadGrado_tv = case_when(
        time >= fupvisit_1 & (is.na(fupvisit_2) | time < fupvisit_2) ~ ToxicidadGrado_eval_1,
        #......................
        time >= fupvisit_22 & (is.na(fupvisit_23) | time < fupvisit_23) ~ ToxicidadGrado_eval_22,
        time >= fupvisit_23 ~ ToxicidadGrado_eval_23,
        TRUE ~ ToxicidadGrado_tv
      )
    )
  
   #LDHnormal time varying
  trial_data_long <- trial_data_long %>%
    mutate(
      LDHnormal_tv = LDHnormal_basalNEW,
      LDHnormal_tv = case_when(
        time >= fupvisit_1 & (is.na(fupvisit_2) | time <= fupvisit_2) ~ LDHnormal_eval_1,
        #......................
        time >= fupvisit_22 & (is.na(fupvisit_23) | time <= fupvisit_23) ~ LDHnormal_eval_22,
        time >= fupvisit_23 ~ LDHnormal_eval_23,
        TRUE ~ LDHnormal_tv
      )
    )
  
  #LDHnormal.cf time varying
  trial_data_long <- trial_data_long %>%
    mutate(
      LDHnormal.cf_tv = LDHnormal_basalNEW,
      LDHnormal.cf_tv = case_when(
        time >= fupvisit_1 & (is.na(fupvisit_2) | time < fupvisit_2) ~ LDHnormal.cf_eval_1,
        #......................
        time >= fupvisit_22 & (is.na(fupvisit_23) | time < fupvisit_23) ~ LDHnormal.cf_eval_22,
        time >= fupvisit_23 ~ LDHnormal.cf_eval_23,
        TRUE ~ LDHnormal.cf_tv
      )
    )
  
  #TIme since last LDH
  trial_data_long <- trial_data_long %>%
    mutate(
      TimesinceLDH_tv = 0,
      TimesinceLDH_tv = case_when(
        time >= fupvisit_1 & (is.na(fupvisit_2) | time < fupvisit_2) ~ TimesinceLDH_1,
        #......................
        time >= fupvisit_22 & (is.na(fupvisit_23) | time < fupvisit_23) ~ TimesinceLDH_22,
        time >= fupvisit_23 ~ TimesinceLDH_23,
        TRUE ~ TimesinceLDH_tv
      )
    )
  
  trial_data_long$TimesinceLDH_tv <- ifelse(trial_data_long$TimesinceLDH_tv == 0, 1, ceiling(trial_data_long$TimesinceLDH_tv/7)) #in weeks, where 0's are included in week 1
  
  
  ## Develop weights----
  
  # Develop weight model (bsl & post-bsl adjusted) 
  # to estimate probability of initiating MAB in the first 15 weeks following QMT1 [model1] 
  # AND probability of initiating MAB in the week of QMT2 [model2]

  
  ### Fit models----
  w.model_vars <- c(
    "time",
    "rcs(time, knots=2)",
    "age", #age(in years),
    "agesq",
    "PS_basal_impNEW", #ECOG PS (0, 1, 2)
    "Charlson2_basalNEW",  #Charlson (<3, >=3)
    "RASMutation_basal", # RAS (mutant, wt, NE)
    "BRAFMutation_basal",  #BRAF (mutant, wt, NE)
    "Microsatel_basal", #Microsatel (msi, mss, NE)
    "LocationPrimaryTumor_basal", #Tumor site (left, right)
    "PrimarySurgery_basal", # Surgery (no, yes)
    "NumberOrgansAffected_basal", #Number of organs (1,>1)
    "Liver2_basal", #(no liver mts, liver mts + elsewhere, only liver mts 1-3 lesions <=5cm, only liver mts 4-9 lesions <=5cm, only liver mts >10 lesions or >5cm)
    "Lung2_basal", #(no lung mts, only lung mts, lung mts + elsewhere)
    "PeritonealAffected_basal", # (no, yes)
    "NodeAffected_basal", # (no, yes)
    "LDHnormal_basalNEW", #(normal, abnormal, NE)
    "PS_tv",  #(0,1,2,3,4)
    "Charlson2_tv", #Charlson (<3, >=3)
    "ToxicidadGrado_tv",  # (no, yes)
    "TimesinceLDH_tv",  #(normal, abnormal, NE)
    "LDHnormal.cf_tv"  #(normal, abnormal, NE)
  )
  
  w.model_vars2 <- c(
    "time",
    "age", #age(in years),
    "agesq",
    "PS_basal_impNEW", #ECOG PS (0, 1, 2)
    "Charlson2_basalNEW",  #Charlson (<3, >=3)
    "RASMutation_basal", # RAS (mutant, wt, NE)
    "BRAFMutation_basal",  #BRAF (mutant, wt, NE)
    "Microsatel_basal", #Microsatel (msi, mss, NE)
    "LocationPrimaryTumor_basal", #Tumor site (left, right)
    "PrimarySurgery_basal", # Surgery (no, yes)
    "NumberOrgansAffected_basal", #Number of organs (1,>1)
    "Liver2_basal", #(no liver mts, liver mts + elsewhere, only liver mts 1-3 lesions <=5cm, only liver mts 4-9 lesions <=5cm, only liver mts >10 lesions or >5cm)
    "Lung2_basal", #(no lung mts, only lung mts, lung mts + elsewhere)
    "PeritonealAffected_basal", # (no, yes)
    "NodeAffected_basal", # (no, yes)
    "LDHnormal_basalNEW", #(normal, abnormal, NE)
    "PS_tv",  #(0,1,2,3,4)
    "Charlson2_tv", #Charlson (<3, >=3)
    "ToxicidadGrado_tv",  # (no, yes)
    "TimesinceLDH_tv",  #(normal, abnormal, NE)
    "LDHnormal.cf_tv"  #(normal, abnormal, NE)
  )
 
  formula <- as.formula(paste("MABincident ~", paste(w.model_vars, collapse = " + ")))
  formula2 <- as.formula(paste("MABincident ~", paste(w.model_vars2, collapse = " + ")))
  # Fit the glm model to estimate weekly probability of initiating MAB1 conditional on baseline and postbaseline vars
  trial_data_arm1 <- subset(trial_data_long, arm == 1)
  trial_data_arm1$PS_basal_impNEW <- relevel(factor(trial_data_arm1$PS_basal_impNEW), ref = "0")
  trial_data_arm1_model1 <- subset(trial_data_arm1, aux_week < 16)
  trial_data_arm1_model2 <- subset(trial_data_arm1, !is.na(QMT2) & startW<=QMT2 & QMT2<=endW)
  cw.fit1 <- glm(formula, data = trial_data_arm1_model1, family = binomial())
  cw.fit2 <- glm(formula2, data = trial_data_arm1_model2, family = binomial())
  weights1[[i]]<-summary(cw.fit1)
  weights2[[i]]<-summary(cw.fit2)
  
  #Inspect and summarise the weekly conditional probabilities of initiating MAB1
  p.mod1.obs <- predict(cw.fit1, newdata = trial_data_arm1_model1, type = "response")
  p.mod2.obs <- predict(cw.fit2, newdata = trial_data_arm1_model2, type = "response")
  
  ### Apply probabilities to each person-week ----
  #APPLY PROBABILITIES FROM 2 MODELS TO SUBJECTS WITH ARM==1 (SUBJECTS WITH ARM==0 WILL WEIGHT 1)
  trial_data_long$prMAB[trial_data_long$arm == 1 & trial_data_long$aux_week<16] <- p.mod1.obs
  trial_data_long$prMAB[trial_data_long$arm == 1 & !is.na(trial_data_long$QMT2) & trial_data_long$startW<=trial_data_long$QMT2 & trial_data_long$QMT2<=trial_data_long$endW] <- p.mod2.obs
  trial_data_long$prMAB[trial_data_long$arm == 1 & is.na(trial_data_long$prMAB)] <- 0

  ### Calculate weight factors----
  
  # Initialize factor.w to NA,
  trial_data_long$factor.w <-NA  
  
  # Calcola factor.w con le nuove regole
  trial_data_long <- trial_data_long %>%
    mutate(
      factor.w = case_when(
        arm == 0 ~ 1,
        arm == 1 & (is.na(QMT2) | (!is.na(QMT2) & QMT2wRetained==0)) ~ 1 / (1 - prMAB),
        arm == 1 & startW<=QMT2 & QMT2<=endW & !is.na(MAB2)~ 1 / prMAB,
        arm == 1 & startW<=QMT2 & QMT2<=endW & is.na(MAB2)~ 0,
        arm == 1 & startW>QMT2 ~ 0
      )
    )
  
  trial_data_long <- trial_data_long %>%
    group_by(subject_custom_code) %>%
    mutate(
      mab_count = cumsum(MAB == 1)  #WE MUST SET FACTOW.W AT 1  FROM THE WEEK AFTER STARTING MAB1
    ) %>%
    ungroup()
  
  trial_data_long <- trial_data_long %>%
    mutate(
      factor.w = ifelse(MAB == 1 & mab_count >=2, 1, factor.w)
    )
  
  trial_data_long <- trial_data_long %>%
    select(-mab_count)
  
    ### Calculate censor weights----
  trial_data_long <- trial_data_long %>%
    group_by(subject_custom_code) %>%
    mutate(cw = cumprod(factor.w))
  
  ### Truncate weights to the 99th percentile ----
  p99<-quantile(trial_data_long$cw, probs=0.99)
  trial_data_long$cw.t99<-ifelse(trial_data_long$cw>p99, p99, trial_data_long$cw)
 
 
  ## Estimate parametric survival curve, baseline and censor weight adjusted ----
  
  # fit of parametric hazards model adjusted for baseline covariates and censoring weights (which include baseline and postbaseline)
  model_vars <- c(
    "arm",
    "time",
    "rcs(time, knots=2)",
    "arm:time",
    "arm:rcs(time, knots=2)",
    "age", #age(in years),
    "agesq",
    "PS_basal_impNEW", #ECOG PS (0, 1, 2)
    "Charlson2_basalNEW",  #Charlson (<3, >=3)
    "RASMutation_basal", # RAS (mutant, wt, NE)
    "BRAFMutation_basal",  #BRAF (mutant, wt, NE)
    "Microsatel_basal", #Microsatel (msi, mss, NE)
    "LocationPrimaryTumor_basal", #Tumor site (left, right)
    "PrimarySurgery_basal", # Surgery (no, yes)
    "NumberOrgansAffected_basal", #Number of organs (1,>1)
    "Liver2_basal", #(no liver mts, liver mts + elsewhere, only liver mts 1-3 lesions <=5cm, only liver mts 4-9 lesions <=5cm, only liver mts >10 lesions or >5cm)
    "Lung2_basal", #(no lung mts, only lung mts, lung mts + elsewhere)
    "PeritonealAffected_basal", # (no, yes)
    "NodeAffected_basal", # (no, yes)
    "LDHnormal_basalNEW" #(normal, abnormal, NE)
  )
  formula <- as.formula(paste("event==0 ~", paste(model_vars, collapse = " + ")))
  full.adj.glm.I<-glm(formula, data = trial_data_long, family = binomial(), weights=cw.t99)
  
  dsc.c0<- trial_data_wide %>%
    filter(arm==0)
  bsl.arm0 <- expandRows(dsc.c0, count=209, count.is.col=F)
  bsl.arm0$time <- rep(seq(1, 209), nrow(dsc.c0))
  
  
  dsc.c1<- trial_data_wide %>%
    filter(arm==1)
  bsl.arm1 <- expandRows(dsc.c1, count=209, count.is.col=F)
  bsl.arm1$time <- rep(seq(1, 209), nrow(dsc.c1))
  
  # assignment of estimated (1-hazard) to each person-week */
  bsl.arm0$p.noevent0 <- predict(full.adj.glm.I, bsl.arm0, type="response") #predict pnoevent at each person-week while each person is followed
  bsl.arm1$p.noevent1 <- predict(full.adj.glm.I, bsl.arm1, type="response")
  
  # computation of survival for each person-week
  bsl.arm0.surv <- bsl.arm0 %>% group_by(subject_custom_code) %>% mutate(surv0 = cumprod(p.noevent0))
  bsl.arm1.surv <- bsl.arm1 %>% group_by(subject_custom_code) %>% mutate(surv1 = cumprod(p.noevent1))
  
  bsl.surv0 <- aggregate(bsl.arm0.surv, by=list(bsl.arm0.surv$time), FUN=mean)[c("arm", "time", "surv0")]
  bsl.surv1 <- aggregate(bsl.arm1.surv, by=list(bsl.arm1.surv$time), FUN=mean)[c("arm", "time", "surv1")]
  
  # computation of risk for each person-week
  bsl.surv0$risk0<-1-bsl.surv0$surv0
  bsl.surv1$risk1<-1-bsl.surv1$surv1
  # some data management to plot estimated survival curves
  bsl.cw.graph <- merge(bsl.surv0, bsl.surv1, by=c("time"))
  bsl.cw.graph$survdiff <- bsl.cw.graph$surv1-bsl.cw.graph$surv0
  bsl.cw.graph$riskratio<- bsl.cw.graph$risk1/bsl.cw.graph$risk0
  bsl.cw.graph$time_mo <- bsl.cw.graph$time / 4.3452  # Time in months
  
  # plot
  bsl.cw.plot<-ggplot(bsl.cw.graph, aes(x=time_mo, y=surv)) + 
    geom_line(aes(y = surv0, colour = "mAB 1st line")) + 
    geom_line(aes(y = surv1, colour = "mAB 2nd line")) + 
    xlab("Months") + 
    scale_x_continuous(limits = c(0, 48), breaks=seq(0,48,6)) +
    scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
    ylab("Survival") + 
    ggtitle(paste("Parametric Survival Curve, baseline and weights- adjusted: Trial ",i,sep="")) + 
    labs(colour="Strategy") +
    theme_bw() + 
    theme(legend.position="bottom")

  ## Combined parametric and non parametric plot for comparison
  surv_obj <- Surv(trial_data_wide$fup_m, trial_data_wide$death) 
  km_fit <- survfit(surv_obj ~ arm, data=trial_data_wide)
  summary(km_fit)
  km_plot<-ggsurvplot(km_fit, data=trial_data_wide,
                      conf.int=F, 
                      risk.table = T,
                      censor=F,
                      title = paste("Trial ",i,sep=""), 
                      xlab = "Months", 
                      ylab = "Survival",
                      xlim = c(0,48),
                      surv.scale="percent",
                      break.x.by = 6,
                      legend.title="Strategy", 
                      legend= "bottom", 
                      legend.labs=c("KM MAB1", "KM MAB2")) 
  
  combined_plot2<- km_plot$plot +
    geom_line(aes(x=time_mo, y = surv0, colour = "PLR MAB1"), data=bsl.cw.graph) +
    geom_line(aes(x=time_mo, y = surv1, colour = "PLR MAB2"), data=bsl.cw.graph) +
    xlab("Months") +
    scale_x_continuous(limits = c(0, 48), breaks=seq(0,48,6)) +
    ylab("Survival") +
    ggtitle("Parametric fit to Non-parametric Survival Curves") +
    labs(colour="Strategy") +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = NA),)+
    theme(legend.position="bottom")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    scale_color_manual(values = c("KM MAB1"= "#ff9999", "PLR MAB1" = "#ff9999", "KM MAB2"="#33CCCC","PLR MAB2" = "#33CCCC"))

  ### Estimate baseline and weights-adjusted hazard ratios----
  full.adj.glm.noI<-glm(formula, data = trial_data_long, family = binomial(), weights=cw.t99)
 
  # Extract the HR and the robust standard error
  hr<-1/exp(coef(full.adj.glm.noI)[["arm"]]) # to get Hazard Ratio for arm
  robust_se<-coeftest(full.adj.glm.noI, vcov=vcovHC(full.adj.glm.noI, type="HC1"))[2,2] # to get robust SE estimates for arm
  
  # Calculate the robust 95% confidence interval
  robust_ci_upper <- 1/exp(log(1/hr) - 1.96 * robust_se)
  robust_ci_lower <- 1/exp(log(1/hr) + 1.96 * robust_se)
  }


#------------------- Pool sequential trials-----------------------------

trial_names_U <- paste0("trial", 1:8,"_U")  
pooled_U <- mget(trial_names_U) %>% bind_rows()
pooled_U <- pooled_U %>%
  mutate(subject_custom_code = paste0(subject_custom_code, "_", trial))
trial_names_W <- paste0("trial", 1:8,"_W")  
pooled_W <- mget(trial_names_W) %>% bind_rows()
pooled_W <- pooled_W %>%
  mutate(subject_custom_code = paste0(subject_custom_code, "_", trial))

for (j in 1:23) {
  pooled_W[[paste0("fupvisit_", j)]]<-ceiling(as.numeric((pooled_W[[paste0("FechaEval_", j)]]-pooled_W$FechaInicioTtoVisitaBasal.3)/7)) #time in weeks
}

## Plot time under follow-up ----
# Figure 1. Time under followup after applying the artificial censoring

#Plot the Time under follow up for clones
surv_obj <- Surv(pooled_W$fup_m) #time in months
km_fit <- survfit(surv_obj ~ arm, data=pooled_W)
tfup_plot<-ggsurvplot(km_fit, data=pooled_W,
                      conf.int=TRUE, 
                      risk.table = TRUE,
                      censor=FALSE,
                      title = "Pooled Trial", 
                      xlab = "Months", 
                      ylab = "Survival",
                      xlim = c(0,48),
                      break.x.by = 6,
                      legend.title="Exposure", 
                      legend= "bottom", 
                      legend.labs=c("mAB 1st line", "mAB 2nd line"))

## Estimate non-parametric survival curves, unadjusted (Kaplan Meier) ----
# Figure 2. Kaplan Meier figure of overall survival after applying the artificial censoring, without any adjustment

surv_obj <- Surv(pooled_W$fup_m, pooled_W$death) 
km_fit <- survfit(surv_obj ~ arm, data=pooled_W)

km_plot<-ggsurvplot(km_fit, data=pooled_W,
                    conf.int=F, 
                    risk.table = T,
                    censor=F,
                    title = "Pooled Trial", 
                    xlab = "Months", 
                    ylab = "Survival",
                    xlim = c(0,48),
                    surv.scale="percent",
                    break.x.by = 6,
                    legend.title="Strategy", 
                    legend= "bottom", 
                    legend.labs=c("KM MAB1", "KM MAB2")) 

## Estimate parametric survival curves, unadjusted ----

# Figure 3. Parametric survival curve (pooled log reg) of overall survival after applying the artificial censoring, without any adjustment
# fit of parametric plr hazards model (modelling the probability of the event NOT occurring)
unadj.glm.I <- glm(event==0 ~ arm + time + I(arm*time) 
                   + rcs(time, knots=2) + arm:rcs(time, knots=2)
                   , family=binomial(), data=pooled_U)
# creation of dataset with all time points under each treatment level
arm0 <- data.frame(cbind(seq(0, 209),0,(seq(0, 209))^2))
arm1 <- data.frame(cbind(seq(0, 209),1,(seq(0, 209))^2))
colnames(arm0) <- c("time", "arm", "timesq")
colnames(arm1) <- c("time", "arm", "timesq")
# assignment of estimated (1-hazard) to each person-week */
arm0$p.noevent0 <- predict(unadj.glm.I, arm0, type="response")
arm1$p.noevent1 <- predict(unadj.glm.I, arm1, type="response")
# computation of survival for each person-week
arm0$surv0 <- cumprod(arm0$p.noevent0) #cumulative sequential product of values in p.noevent0
arm1$surv1 <- cumprod(arm1$p.noevent1)
#computation of risk for each person-week
arm0$risk0<-1-arm0$surv0
arm1$risk1<-1-arm1$surv1
# some data management to plot estimated survival curves
unadj.graph <- merge(arm0, arm1, by=c("time", "timesq"))
unadj.graph$survdiff <- unadj.graph$surv1-unadj.graph$surv0
unadj.graph$riskratio<- unadj.graph$risk1/unadj.graph$risk0
unadj.graph$time_mo <- unadj.graph$time / 4.3452  # Time in months
# plot
unadj.plot<-ggplot(unadj.graph, aes(x=time_mo, y=surv)) + 
  geom_line(aes(y = surv0, colour = "mAB 1st line")) + 
  geom_line(aes(y = surv1, colour = "mAB 2nd line")) + 
  xlab("Months") + 
  scale_x_continuous(limits = c(0, 48), breaks=seq(0,48,6)) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  ylab("Survival") + 
  ggtitle("Pooled Trial") + 
  labs(colour="Strategy") +
  theme_bw() + 
  theme(legend.position="bottom")
### Print effect estimates:unadjusted ----
print("Pooled Trial")
# 48 month survival for MAB1 and MAB2
print(unadj.graph[unadj.graph$time == 209, c("surv0", "surv1")])
# Survival difference at 48 months
print(unadj.graph[unadj.graph$time == 209, c("survdiff")])
#Risk ratio at 48 months
print(unadj.graph[unadj.graph$time == 209, c("riskratio")])
# Combined parametric and non parametric plot for comparison
combined_plot<- km_plot$plot +
  #ggplot(hazards.graph, aes(x=time, y=surv)) + 
  geom_line(aes(x=time_mo, y = surv0, colour = "PLR MAB1"), data=unadj.graph) + 
  geom_line(aes(x=time_mo, y = surv1, colour = "PLR MAB2"), data=unadj.graph) + 
  xlab("Months") + 
  scale_x_continuous(limits = c(0, 48), breaks=seq(0,48,6)) +
  ylab("Survival") + 
  ggtitle("Pooled Trial") + 
  labs(colour="Strategy") +
  theme_minimal() + 
  theme(plot.background = element_rect(fill = "white", color = NA),)+
  theme(legend.position="bottom")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(values = c("KM MAB1"= "#ff9999", "PLR MAB1" = "#ff9999", "KM MAB2"="#33CCCC","PLR MAB2" = "#33CCCC"))


### Estimate unadjusted hazard ratios ----
#To estimate the Hazard ratios, we need to run the same GLM model but without interactions between arm and time
unadj.glm.noI <- glm(event==0 ~ arm + time  + rcs(time, knots=2), family=binomial(), data=pooled_U)
# Extract the HR and the robust standard error
hr<-1/exp(coef(unadj.glm.noI)[["arm"]]) # to get Hazard Ratio for arm
robust_se<-coeftest(unadj.glm.noI, vcov=vcovHC(unadj.glm.noI, type="HC1"))[2,2] # to get robust SE estimates for arm
robust_ci_upper <- 1/exp(log(1/hr) - 1.96 * robust_se)
robust_ci_lower <- 1/exp(log(1/hr) + 1.96 * robust_se)

#Effects: HR and their CI
print(paste("Pooled Trial: Hazard Ratio (HR) [Unadj] =", hr))
print(paste("Pooled Trial: Lower 95% CI [Unadj]=", robust_ci_lower))
print(paste("Pooled Trial: Upper 95% CI [Unadj]=", robust_ci_upper))

## Estimate parametric survival curves, baseline adjusted ----

# Figure 4. Parametric survival curve of overall survival after applying the artificial censoring, adjusted for baseline variables

# fit of parametric hazards model with covariates
bsl.adj.glm.I <- glm(event==0 ~ arm 
                     + rcs(time, knots=5) 
                     + arm:rcs(time, knots=5)
                     + age #age(in years)
                     + agesq
                     + PS_basal_impNEW #ECOG PS (0, 1, 2)
                     + Charlson2_basalNEW #Charlson (<3, >=3)
                     + RASMutation_basal  # RAS (mutant, wt, NE)
                     + BRAFMutation_basal #BRAF (mutant, wt, NE)
                     + Microsatel_basal #Microsatel (msi, mss, NE)
                     + LocationPrimaryTumor_basal  #Tumor site (left, right)
                     + PrimarySurgery_basal # Surgery (no, yes)
                     + NumberOrgansAffected_basal #Number of organs (1,>1)
                     + Liver2_basal #(no liver mts, liver mts + elsewhere, only liver mts 1-3 lesions <=5cm, only liver mts 4-9 lesions <=5cm, only liver mts >10 lesions or >5cm)
                     + Lung2_basal #(no lung mts, only lung mts, lung mts + elsewhere)
                     + PeritonealAffected_basal # (no, yes)
                     + NodeAffected_basal # (no, yes)
                     + LDHnormal_basalNEW #(normal, abnormal, NE)
                     + trial 
                     , data=pooled_U, family=binomial())

dsc.c0<- pooled_W %>%
  filter(arm==0)
bsl.arm0 <- expandRows(dsc.c0, count=209, count.is.col=F)
bsl.arm0$time <- rep(seq(1, 209), nrow(dsc.c0))

dsc.c1<- pooled_W %>%
  filter(arm==1)
bsl.arm1 <- expandRows(dsc.c1, count=209, count.is.col=F)
bsl.arm1$time <- rep(seq(1, 209), nrow(dsc.c1))

# assignment of estimated (1-hazard) to each person-week */
bsl.arm0$p.noevent0 <- predict(bsl.adj.glm.I, bsl.arm0, type="response") #predict pnoevent at each person-week while each person is followed
bsl.arm1$p.noevent1 <- predict(bsl.adj.glm.I, bsl.arm1, type="response")

# computation of survival for each person-week
bsl.arm0.surv <- bsl.arm0 %>% group_by(subject_custom_code) %>% mutate(surv0 = cumprod(p.noevent0))
bsl.arm1.surv <- bsl.arm1 %>% group_by(subject_custom_code) %>% mutate(surv1 = cumprod(p.noevent1))

bsl.surv0 <- aggregate(bsl.arm0.surv, by=list(bsl.arm0.surv$time), FUN=mean)[c("arm", "time", "surv0")]
bsl.surv1 <- aggregate(bsl.arm1.surv, by=list(bsl.arm1.surv$time), FUN=mean)[c("arm", "time", "surv1")]

# computation of risk for each person-week
bsl.surv0$risk0<-1-bsl.surv0$surv0
bsl.surv1$risk1<-1-bsl.surv1$surv1

# some data management to plot estimated survival curves
bsl.graph <- merge(bsl.surv0, bsl.surv1, by=c("time"))
bsl.graph$survdiff <- bsl.graph$surv1-bsl.graph$surv0
bsl.graph$riskratio<- bsl.graph$risk1/bsl.graph$risk0
bsl.graph$time_mo <- bsl.graph$time / 4.3452  # Time in months

# plot
bsl.plot<-ggplot(bsl.graph, aes(x=time_mo, y=surv)) + 
  geom_line(aes(y = surv0, colour = "mAB 1st line")) + 
  geom_line(aes(y = surv1, colour = "mAB 2nd line")) + 
  xlab("Months") + 
  scale_x_continuous(limits = c(0, 48), breaks=seq(0,48,6)) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  ylab("Survival") + 
  ggtitle("Pooled Trial") + 
  labs(colour="Strategy") +
  theme_bw() + 
  theme(legend.position="bottom")

### Print effect estimates: baseline adjusted ----
# 48 month survival for MAB1 and MAB2
print(bsl.graph[bsl.graph$time == 209, c("surv0", "surv1")])
# Survival difference at 48 months
print(bsl.graph[bsl.graph$time == 209, c("survdiff")])
#Risk ratio at 48 months
print(bsl.graph[bsl.graph$time == 209, c("riskratio")])

### Estimate baseline adjusted hazard ratios----
# fit of parametric hazards model with covariates
bsl.adj.glm.noI <- glm(event==0 ~ arm
                       + rcs(time, knots=2) 
                       + age #age(in years)
                       + agesq
                       + PS_basal_impNEW #ECOG PS (0, 1, 2, and NAs)
                       + Charlson2_basalNEW #Charlson (<3, >=3)
                       + RASMutation_basal  # RAS (mutant, wt, NE)
                       + BRAFMutation_basal #BRAF (mutant, wt, NE)
                       + Microsatel_basal #Microsatel (msi, mss, NE)
                       + LocationPrimaryTumor_basal  #Tumor site (left, right)
                       + PrimarySurgery_basal # Surgery (no, yes)
                       + NumberOrgansAffected_basal #Number of organs (1,>1)
                       + Liver2_basal #(no liver mts, liver mts + elsewhere, only liver mts 1-3 lesions <=5cm, only liver mts 4-9 lesions <=5cm, only liver mts >10 lesions or >5cm)
                       + Lung2_basal #(no lung mts, only lung mts, lung mts + elsewhere)
                       + PeritonealAffected_basal # (no, yes)
                       + NodeAffected_basal # (no, yes)
                       + LDHnormal_basalNEW
                       + trial #(normal, abnormal, NE)
                       , data=pooled_U, family=binomial())

# Extract the HR and the robust standard error
hr<-1/exp(coef(bsl.adj.glm.noI)[["arm"]]) # to get Hazard Ratio for arm
robust_se<-coeftest(bsl.adj.glm.noI, vcov=vcovHC(bsl.adj.glm.noI, type="HC1"))[2,2]
robust_ci_upper <- 1/exp(log(1/hr) - 1.96 * robust_se)
robust_ci_lower <- 1/exp(log(1/hr) + 1.96 * robust_se)

## Estimate parametric survival curves, baseline and weights-adjusted ----
# fit of parametric hazards model adjusted for baseline covariates and censoring weights (which include baseline and postbaseline)
model_vars <- c(
  "arm",
  "time",
  "rcs(time, knots=2)",
  "arm:time",
  "arm:rcs(time, knots=2)",
  "age", #age(in years),
  "agesq",
  "PS_basal_impNEW", #ECOG PS (0, 1, 2)
  "Charlson2_basalNEW",  #Charlson (<3, >=3)
  "RASMutation_basal", # RAS (mutant, wt, NE)
  "BRAFMutation_basal",  #BRAF (mutant, wt, NE)
  "Microsatel_basal", #Microsatel (msi, mss, NE)
  "LocationPrimaryTumor_basal", #Tumor site (left, right)
  "PrimarySurgery_basal", # Surgery (no, yes)
  "NumberOrgansAffected_basal", #Number of organs (1,>1)
  "Liver2_basal", #(no liver mts, liver mts + elsewhere, only liver mts 1-3 lesions <=5cm, only liver mts 4-9 lesions <=5cm, only liver mts >10 lesions or >5cm)
  "Lung2_basal", #(no lung mts, only lung mts, lung mts + elsewhere)
  "PeritonealAffected_basal", # (no, yes)
  "NodeAffected_basal", # (no, yes)
  "LDHnormal_basalNEW", #(normal, abnormal, NE)
  "trial"
)
formula <- as.formula(paste("event==0 ~", paste(model_vars, collapse = " + ")))
full.adj.glm.I<-glm(formula, data = pooled_U, family = binomial(), weights=cw.t99)
print(summary(full.adj.glm.I))

dsc.c0<- pooled_W %>%
  filter(arm==0)
bsl.arm0 <- expandRows(dsc.c0, count=209, count.is.col=F)
bsl.arm0$time <- rep(seq(1, 209), nrow(dsc.c0))

dsc.c1<- pooled_W %>%
  filter(arm==1)
bsl.arm1 <- expandRows(dsc.c1, count=209, count.is.col=F)
bsl.arm1$time <- rep(seq(1, 209), nrow(dsc.c1))

# assignment of estimated (1-hazard) to each person-week */
bsl.arm0$p.noevent0 <- predict(full.adj.glm.I, bsl.arm0, type="response") #predict pnoevent at each person-week while each person is followed
bsl.arm1$p.noevent1 <- predict(full.adj.glm.I, bsl.arm1, type="response")

# computation of survival for each person-week
bsl.arm0.surv <- bsl.arm0 %>% group_by(subject_custom_code) %>% mutate(surv0 = cumprod(p.noevent0))
bsl.arm1.surv <- bsl.arm1 %>% group_by(subject_custom_code) %>% mutate(surv1 = cumprod(p.noevent1))

bsl.surv0 <- aggregate(bsl.arm0.surv, by=list(bsl.arm0.surv$time), FUN=mean)[c("arm", "time", "surv0")]
bsl.surv1 <- aggregate(bsl.arm1.surv, by=list(bsl.arm1.surv$time), FUN=mean)[c("arm", "time", "surv1")]

# computation of risk for each person-week
bsl.surv0$risk0<-1-bsl.surv0$surv0
bsl.surv1$risk1<-1-bsl.surv1$surv1
# some data management to plot estimated survival curves
bsl.cw.graph <- merge(bsl.surv0, bsl.surv1, by=c("time"))
bsl.cw.graph$survdiff <- bsl.cw.graph$surv1-bsl.cw.graph$surv0
bsl.cw.graph$riskratio<- bsl.cw.graph$risk1/bsl.cw.graph$risk0
bsl.cw.graph$time_mo <- bsl.cw.graph$time / 4.3452  # Time in months

# plot
bsl.cw.plot<-ggplot(bsl.cw.graph, aes(x=time_mo, y=surv)) + 
  geom_line(aes(y = surv0, colour = "mAB 1st line")) + 
  geom_line(aes(y = surv1, colour = "mAB 2nd line")) + 
  xlab("Months") + 
  scale_x_continuous(limits = c(0, 48), breaks=seq(0,48,6)) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  ylab("Survival") + 
  ggtitle("Parametric Survival Curve, baseline and weights- adjusted: Pooled Trial") + 
  labs(colour="Strategy") +
  theme_bw() + 
  theme(legend.position="bottom")

# Combined parametric and non parametric plot for comparison
surv_obj <- Surv(pooled_W$fup_m, pooled_W$death) 
km_fit <- survfit(surv_obj ~ arm, data=pooled_W)
summary(km_fit)
km_plot<-ggsurvplot(km_fit, data=pooled_W,
                    conf.int=F, 
                    risk.table = T,
                    censor=F,
                    title = "Pooled Trial", 
                    xlab = "Months", 
                    ylab = "Survival",
                    xlim = c(0,48),
                    surv.scale="percent",
                    break.x.by = 6,
                    legend.title="Strategy", 
                    legend= "bottom", 
                    legend.labs=c("KM MAB1", "KM MAB2")) 

combined_plot2<- km_plot$plot +
  geom_line(aes(x=time_mo, y = surv0, colour = "PLR MAB1"), data=bsl.cw.graph) +
  geom_line(aes(x=time_mo, y = surv1, colour = "PLR MAB2"), data=bsl.cw.graph) +
  xlab("Months") +
  scale_x_continuous(limits = c(0, 48), breaks=seq(0,48,6)) +
  ylab("Survival") +
  ggtitle("Parametric fit to Non-parametric Survival Curves: Pooled Trial") +
  labs(colour="Strategy") +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white", color = NA),)+
  theme(legend.position="bottom")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(values = c("KM MAB1"= "#ff9999", "PLR MAB1" = "#ff9999", "KM MAB2"="#33CCCC","PLR MAB2" = "#33CCCC"))

### Print effect estimates: baseline and weights-adjusted ----

print("Completely adjusted")
# Adjusted for  baseline variables and censoring time-varying weights 
# 48 month survival for MAB1 and MAB2
print(bsl.cw.graph[bsl.cw.graph$time == 209, c("surv0", "surv1")])
# Survival difference at 48 months
print(bsl.cw.graph[bsl.cw.graph$time == 209, c("survdiff")])
#Risk ratio at 48 months
print(bsl.cw.graph[bsl.cw.graph$time == 209, c("riskratio")])

### Estimate baseline and weights-adjusted hazards ratio ----
#Baseline-adj Hazard ratios from GLM model without interaction at the different time points


  # Fit glm model without I
  model_vars <- c(
    "arm",
    "time",
    "rcs(time, knots=2)",
    "age", #age(in years),
    "agesq",
    "PS_basal_impNEW", #ECOG PS (0, 1, 2)
    "Charlson2_basalNEW",  #Charlson (<3, >=3)
    "RASMutation_basal", # RAS (mutant, wt, NE)
    "BRAFMutation_basal",  #BRAF (mutant, wt, NE)
    "Microsatel_basal", #Microsatel (msi, mss, NE)
    "LocationPrimaryTumor_basal", #Tumor site (left, right)
    "PrimarySurgery_basal", # Surgery (no, yes)
    "NumberOrgansAffected_basal", #Number of organs (1,>1)
    "Liver2_basal", #(no liver mts, liver mts + elsewhere, only liver mts 1-3 lesions <=5cm, only liver mts 4-9 lesions <=5cm, only liver mts >10 lesions or >5cm)
    "Lung2_basal", #(no lung mts, only lung mts, lung mts + elsewhere)
    "PeritonealAffected_basal", # (no, yes)
    "NodeAffected_basal", # (no, yes)
    "LDHnormal_basalNEW",#(normal, abnormal, NE)
    "trial"
  )
  formula <- as.formula(paste("event==0 ~", paste(model_vars, collapse = " + ")))

full.adj.glm.noI<-glm(formula, data = pooled_U, family = binomial(), weights=cw.t99)
print(summary(full.adj.glm.noI))
# Extract the HR and the robust standard error
hr<-1/exp(coef(full.adj.glm.noI)[["arm"]]) # to get Hazard Ratio for arm
robust_se<-coeftest(full.adj.glm.noI, vcov=vcovHC(full.adj.glm.noI, type="HC1"))[2,2] # to get robust SE estimates for arm

# Calculate the robust 95% confidence interval
robust_ci_lower <- 1/exp(log(1/hr) - 1.96 * robust_se)
robust_ci_upper <- 1/exp(log(1/hr) + 1.96 * robust_se)

#------------------- Bootstrap sequential trials-----------------------------

calculate_ci <- function(values) {
  round(quantile(values, probs = c(0.025, 0.975), na.rm = TRUE), 3)
}

##Bootstraps for the unadjusted model ----


btsp_input_ds<-pooled_W #Store the pooled dataset in the wide data format

boot.func <- function(data, indices, file_name = "bootstrapUnadjW.xlsx") {
  dsc<-data[indices,]
  #Create a bootstrap-specific subject custom code:
  dsc$subject_custom_code_bts<-ave(dsc$subject_custom_code, dsc$subject_custom_code, FUN = function(i) paste0(i, '_', seq_along(i)))

  #Parametric Unadjusted Surv Curves
  dsc$fupweeks<- ceiling(as.numeric(ifelse(dsc$fup==0, 1, dsc$fup/7)))
  dsc.surv<-expandRows(dsc, "fupweeks", drop=F)
  dsc.surv$time <- sequence(rle(dsc.surv$subject_custom_code_bts)$lengths)
  dsc.surv$time_m <- dsc.surv$time/4.3452
  dsc.surv$event <- ifelse(dsc.surv$time==dsc.surv$fupweeks &
                             dsc.surv$death==1, 1, 0)
  dsc.surv$arm<- as.numeric(recode(dsc.surv$arm, "MAB1" = 0, "MAB2" = 1))
   unadj.glm.I <- glm(event==0 ~ arm + time + I(arm*time)
                     + rcs(time, knots=2) + arm:rcs(time, knots=2)
                     , family=binomial(), data=dsc.surv)
    arm0 <- data.frame(cbind(seq(1, 209),0,(seq(1, 209))^2))
  arm1 <- data.frame(cbind(seq(1, 209),1,(seq(1, 209))^2))
  colnames(arm0) <- c("time", "arm", "timesq")
  colnames(arm1) <- c("time", "arm", "timesq")
  arm0$p.noevent0 <- predict(unadj.glm.I, arm0, type="response")
  arm1$p.noevent1 <- predict(unadj.glm.I, arm1, type="response")
  arm0$surv0 <- cumprod(arm0$p.noevent0)
  arm1$surv1 <- cumprod(arm1$p.noevent1)
  arm0$risk0<-1-arm0$surv0
  arm1$risk1<-1-arm1$surv1
  unadj.graph <- merge(arm0, arm1, by=c("time", "timesq"))
  unadj.graph$survdiff <- unadj.graph$surv1-unadj.graph$surv0
  unadj.graph$riskratio<- unadj.graph$risk1/unadj.graph$risk0
  unadj.graph$time_mo <- unadj.graph$time / 4.3452  # Time in months





###########WE BUILT BOOTSTRAP FUNCTIONS TO BE STOPPABLE AND RE-RUNNABLE, IN CASE OF NEEDING

  # Save unadj.graph to an Excel sheet
  sheet_name_prefix <- "Bootstrap_"

  if (!file.exists(file_name)) {
    # Create a new workbook and add the first sheet
    wb <- createWorkbook()
    sheet_name <- paste0(sheet_name_prefix, "1")
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, unadj.graph)
    saveWorkbook(wb, file_name, overwrite = TRUE)
  } else {
    # Append to existing workbook
    wb <- loadWorkbook(file_name)
    existing_sheets <- sheets(wb)
    # Find the next available sheet name
    next_sheet_number <- length(existing_sheets) + 1
    sheet_name <- paste0(sheet_name_prefix, next_sheet_number)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, unadj.graph)
    saveWorkbook(wb, file_name, overwrite = TRUE)
  }

  return(NULL) # Return NULL since the output is saved to a file
}

######IF, FOR ANY REASONS, BOOTSTRAP STOPPED AT i=150 (YOU CAN TELL IT LOOKING AT THE .xlsx FILE AUTOMATICALLY CREATED), YOU JUST HAVE TO PUT i=151 AND RE-RUN IT
for (i in 1:500) {
  set.seed(100 + i)
  indices <- sample(nrow(btsp_input_ds), replace = TRUE)
  boot.func(data = btsp_input_ds, indices = indices, file_name = "bootstrapUnadjW.xlsx")
}

file_path <- "bootstrapUnadjW.xlsx"

sheet_names <- excel_sheets(file_path)

all_data <- lapply(seq_along(sheet_names), function(i) {
  data <- read_excel(file_path, sheet = sheet_names[i])
  data <- data %>% mutate(NBootstrap = i)
  return(data)
})

final_data <- bind_rows(all_data)


ci_results <- final_data %>%
  group_by(time) %>%
  summarise(
    surv095CILow = calculate_ci(surv0)[1],
    surv095CIUp = calculate_ci(surv0)[2],
    risk095CILow = calculate_ci(risk0)[1],
    risk095CIUp = calculate_ci(risk0)[2],
    surv195CILow = calculate_ci(surv1)[1],
    surv195CIUp = calculate_ci(surv1)[2],
    risk195CILow = calculate_ci(risk1)[1],
    risk195CIUp = calculate_ci(risk1)[2],
    survdiff95CILow = calculate_ci(survdiff)[1],
    survdiff95CIUp = calculate_ci(survdiff)[2],
    RiskRatio95CILow = calculate_ci(riskratio)[1],
    RiskRatio95CIUp = calculate_ci(riskratio)[2]
  ) %>%
  ungroup()

unadj.graph2 <- unadj.graph %>%
  left_join(ci_results, by = "time")



unadj.plot <- ggplot(unadj.graph2, aes(x = time_mo)) +
  geom_line(aes(y = surv0, colour = "mAB 1st line")) +
  geom_line(aes(y = surv1, colour = "mAB 2nd line")) +
  #DASHED CIs
  geom_ribbon(aes(ymin = surv095CILow, ymax = surv095CIUp, fill = "mAB 1st line"), alpha = 0.2, linetype = "dashed") +
  geom_ribbon(aes(ymin = surv195CILow, ymax = surv195CIUp, fill = "mAB 2nd line"), alpha = 0.2, linetype = "dashed") +
  xlab("Months") +
  scale_x_continuous(limits = c(0, 48), breaks = seq(0, 48, 6)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  ylab("Survival") +
  ggtitle("Pooled Trial") +
  labs(colour = "Strategy", fill = "Confidence Interval") +
  theme_bw() +
  theme(legend.position = "bottom")




##Bootstraps for the baseline adjusted model----

btsp_input_ds<-pooled_W

# This function has the codes of cloning.R and Param.Bsl.Adj.R,
boot.func <- function(data, indices, file_name = "bootstrapAdjW.xlsx") {
  dsc<-data[indices,]


  #Create a bootstrap-specific subject custom code:
  dsc$subject_custom_code_bts<-ave(dsc$subject_custom_code, dsc$subject_custom_code, FUN = function(i) paste0(i, '_', seq_along(i)))

  #Parametric baseline adjusted Surv Curves
  dsc.c<-dsc
  dsc.c$fupweeks <- ceiling(as.numeric(dsc.c$fup)/7) #fup in weeks
  dsc.c$arm <- as.numeric(recode(dsc.c$arm, "MAB1" = 0, "MAB2" = 1)) #change arm to integer to allow for interaction term
  dsc.c.surv <- expandRows(dsc.c, "fupweeks", drop=F) #If a patient has fup= 35, expandRows will create 35 for this patient with its info repeated
  dsc.c.surv$time <- sequence(rle(dsc.c.surv$subject_custom_code_bts)$lengths)
  dsc.c.surv$time_m<- dsc.c.surv$time/4.3452 #time in months
  dsc.c.surv$event <- ifelse(dsc.c.surv$time==dsc.c.surv$fupweeks &
                               dsc.c.surv$death==1, 1, 0) #create event variable

  bsl.adj.glm.I <- glm(event==0 ~ arm
                       + time
                       + rcs(time, knots=2)
                       + arm:time
                       + arm:rcs(time, knots=2)
                       + age #age(in years)
                       + agesq
                       + PS_basal_impNEW #ECOG PS (0, 1, 2)
                       + Charlson2_basalNEW #Charlson (<3, >=3)
                       + RASMutation_basal  # RAS (mutant, wt, NE)
                       + BRAFMutation_basal #BRAF (mutant, wt, NE)
                       + Microsatel_basal #Microsatel (msi, mss, NE)
                       + LocationPrimaryTumor_basal  #Tumor site (left, right)
                       + PrimarySurgery_basal # Surgery (no, yes)
                       + NumberOrgansAffected_basal #Number of organs (1,>1)
                       + Liver2_basal #(no liver mts, liver mts + elsewhere, only liver mts 1-3 lesions <=5cm, only liver mts 4-9 lesions <=5cm, only liver mts >10 lesions or >5cm)
                       + Lung2_basal #(no lung mts, only lung mts, lung mts + elsewhere)
                       + PeritonealAffected_basal # (no, yes)
                       + NodeAffected_basal # (no, yes)
                       + LDHnormal_basalNEW
                       + trial #(normal, abnormal, NE)
                       , data=dsc.c.surv, family=binomial())





  dsc.c0<- dsc %>%
    filter(arm==0)
  bsl.arm0 <- expandRows(dsc.c0, count=209, count.is.col=F)
  bsl.arm0$time <- rep(seq(1, 209), nrow(dsc.c0))

  dsc.c1<- dsc %>%
    filter(arm==1)
  bsl.arm1 <- expandRows(dsc.c1, count=209, count.is.col=F)
  bsl.arm1$time <- rep(seq(1, 209), nrow(dsc.c1))

  # assignment of estimated (1-hazard) to each person-week */
  bsl.arm0$p.noevent0 <- predict(bsl.adj.glm.I, bsl.arm0, type="response") #predict pnoevent at each person-week while each person is followed
  bsl.arm1$p.noevent1 <- predict(bsl.adj.glm.I, bsl.arm1, type="response")

  # computation of survival for each person-week
  bsl.arm0.surv <- bsl.arm0 %>% group_by(subject_custom_code_bts) %>% mutate(surv0 = cumprod(p.noevent0))
  bsl.arm1.surv <- bsl.arm1 %>% group_by(subject_custom_code_bts) %>% mutate(surv1 = cumprod(p.noevent1))

  bsl.surv0 <- aggregate(bsl.arm0.surv, by=list(bsl.arm0.surv$time), FUN=mean)[c("arm", "time", "surv0")]
  bsl.surv1 <- aggregate(bsl.arm1.surv, by=list(bsl.arm1.surv$time), FUN=mean)[c("arm", "time", "surv1")]

  # computation of risk for each person-week
  bsl.surv0$risk0<-1-bsl.surv0$surv0
  bsl.surv1$risk1<-1-bsl.surv1$surv1
  bsl.graph <- merge(bsl.surv0, bsl.surv1, by=c("time"))
  bsl.graph$survdiff <- bsl.graph$surv1-bsl.graph$surv0
  bsl.graph$riskratio<- bsl.graph$risk1/bsl.graph$risk0
  bsl.graph$time_mo <- bsl.graph$time / 4.3452  # Time in months




  # Save bsl.graph to an Excel sheet
  sheet_name_prefix <- "Bootstrap_"

  if (!file.exists(file_name)) {
    # Create a new workbook and add the first sheet
    wb <- createWorkbook()
    sheet_name <- paste0(sheet_name_prefix, "1")
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, bsl.graph)
    saveWorkbook(wb, file_name, overwrite = TRUE)
  } else {
    # Append to existing workbook
    wb <- loadWorkbook(file_name)
    existing_sheets <- sheets(wb)
    # Find the next available sheet name
    next_sheet_number <- length(existing_sheets) + 1
    sheet_name <- paste0(sheet_name_prefix, next_sheet_number)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, bsl.graph)
    saveWorkbook(wb, file_name, overwrite = TRUE)
  }

  return(NULL) # Return NULL since the output is saved to a file
}


for (i in 1:500) {
  set.seed(100 + i)
  indices <- sample(nrow(btsp_input_ds), replace = TRUE)
  boot.func(data = btsp_input_ds, indices = indices, file_name = "bootstrapAdjW.xlsx")
}

file_path <- "bootstrapAdjW.xlsx"

sheet_names <- excel_sheets(file_path)

all_data <- lapply(seq_along(sheet_names), function(i) {
  data <- read_excel(file_path, sheet = sheet_names[i])
  data <- data %>% mutate(NBootstrap = i)
  return(data)
})

final_data <- bind_rows(all_data)

ci_results <- final_data %>%
  group_by(time) %>%
  summarise(
    surv095CILow = calculate_ci(surv0)[1],
    surv095CIUp = calculate_ci(surv0)[2],
    risk095CILow = calculate_ci(risk0)[1],
    risk095CIUp = calculate_ci(risk0)[2],
    surv195CILow = calculate_ci(surv1)[1],
    surv195CIUp = calculate_ci(surv1)[2],
    risk195CILow = calculate_ci(risk1)[1],
    risk195CIUp = calculate_ci(risk1)[2],
    survdiff95CILow = calculate_ci(survdiff)[1],
    survdiff95CIUp = calculate_ci(survdiff)[2],
    RiskRatio95CILow = calculate_ci(riskratio)[1],
    RiskRatio95CIUp = calculate_ci(riskratio)[2]
  ) %>%
  ungroup()

bsl.graph2 <- bsl.graph %>%
  left_join(ci_results, by = "time")

# plot
bsl.plot<-ggplot(bsl.graph2, aes(x=time_mo)) +
  geom_line(aes(y = surv0, colour = "mAB 1st line")) +
  geom_line(aes(y = surv1, colour = "mAB 2nd line")) +
  #DASHED CIs
  geom_ribbon(aes(ymin = surv095CILow, ymax = surv095CIUp, fill = "mAB 1st line"), alpha = 0.2, linetype = "dashed") +
  geom_ribbon(aes(ymin = surv195CILow, ymax = surv195CIUp, fill = "mAB 2nd line"), alpha = 0.2, linetype = "dashed") +
  xlab("Months") +
  scale_x_continuous(limits = c(0, 48), breaks=seq(0,48,6)) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  ylab("Survival") +
  ggtitle("Pooled Trial") +
  labs(colour = "Strategy", fill = "Confidence Interval") +
  theme_bw() +
  theme(legend.position="bottom")

ggsave(file.path("C:/Users/mzamparini/OneDrive - Research Triangle Institute/Desktop/GEMCAD/output",
                 "Parametricsurvival_bdladj_Pooled_Trial_CIs.png"),bsl.plot,
       width = 8, height = 6)




##Bootstraps for the fully adjusted model ----


btsp_input_ds<-pooled_W 
boot.func <- function(data, indices, file_name = "bootstrapFullAdjW975.xlsx") { 
  dsc<-data[indices,]
  #Create a bootstrap-specific subject custom code:
  dsc$subject_custom_code_bts<-ave(dsc$subject_custom_code, dsc$subject_custom_code, FUN = function(i) paste0(i, '_', seq_along(i)))
  dsclong <- expandRows(dsc, "fup_w", drop=F) 
  dsclong$time <- sequence(rle(dsclong$subject_custom_code_bts)$lengths)
  dsclong$event <- ifelse(dsclong$time==dsclong$fup_w &
                            dsclong$death==1, 1, 0) #create event variable
 
  dsclong <- dsclong %>%
    left_join(pooled_U %>% select(subject_custom_code, aux_week, cw.t99), 
              by = c("subject_custom_code" = "subject_custom_code", "time" = "aux_week"))
  
  model_vars <- c(
    "arm",
    "time",
    "rcs(time, knots=2)",
    "arm:time",
    "arm:rcs(time, knots=2)",
    "age", 
    "agesq",
    "PS_basal_impNEW", 
    "Charlson2_basalNEW",  
    "RASMutation_basal",
    "BRAFMutation_basal", 
    "Microsatel_basal",
    "LocationPrimaryTumor_basal", 
    "PrimarySurgery_basal", 
    "NumberOrgansAffected_basal", 
    "Liver2_basal", 
    "Lung2_basal", 
    "PeritonealAffected_basal", 
    "NodeAffected_basal", 
    "LDHnormal_basalNEW",
    "trial"
  )
  formula <- as.formula(paste("event==0 ~", paste(model_vars, collapse = " + ")))
  full.adj.glm.I<-glm(formula, data = dsclong, family = binomial(), weights=cw.t99)
  dsc.c0<- dsc %>%
    filter(arm==0)
  bsl.cw.arm0 <- expandRows(dsc.c0, count=209, count.is.col=F)
  bsl.cw.arm0$time <- rep(seq(1, 209), nrow(dsc.c0))
  bsl.cw.arm0$timesq <- bsl.cw.arm0$time^2
  bsl.cw.arm0$arm <- 0
  dsc.c1<- dsc %>%
    filter(arm==1)
  bsl.cw.arm1 <- expandRows(dsc.c1, count=209, count.is.col=F)
  bsl.cw.arm1$time <- rep(seq(1, 209), nrow(dsc.c1))
  bsl.cw.arm1$timesq <- bsl.cw.arm1$time^2
  bsl.cw.arm1$arm <- 1
  bsl.cw.arm0$p.noevent0 <- predict(full.adj.glm.I, bsl.cw.arm0, type="response")
  bsl.cw.arm1$p.noevent1 <- predict(full.adj.glm.I, bsl.cw.arm1, type="response")
  bsl.cw.arm0.surv <- bsl.cw.arm0 %>% group_by(subject_custom_code_bts) %>% mutate(surv0 = cumprod(p.noevent0))
  bsl.cw.arm1.surv <- bsl.cw.arm1 %>% group_by(subject_custom_code_bts) %>% mutate(surv1 = cumprod(p.noevent1))
  bsl.cw.surv0 <- aggregate(bsl.cw.arm0.surv, by=list(bsl.cw.arm0.surv$time), FUN=mean)[c("arm", "time", "surv0")]
  bsl.cw.surv1 <- aggregate(bsl.cw.arm1.surv, by=list(bsl.cw.arm1.surv$time), FUN=mean)[c("arm", "time", "surv1")]
  bsl.cw.surv0$risk0<-1-bsl.cw.surv0$surv0
  bsl.cw.surv1$risk1<-1-bsl.cw.surv1$surv1
  bsl.cw.graph <- merge(bsl.cw.surv0, bsl.cw.surv1, by=c("time"))
  bsl.cw.graph$survdiff <- bsl.cw.graph$surv1-bsl.cw.graph$surv0
  bsl.cw.graph$riskratio<- bsl.cw.graph$risk1/bsl.cw.graph$risk0
  bsl.cw.graph$time_mo <- bsl.cw.graph$time / 4.3452  

    # Save bsl.cw.graph to an Excel sheet
  sheet_name_prefix <- "Bootstrap_"
  
  if (!file.exists(file_name)) {
    # Create a new workbook and add the first sheet
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
  
  return(NULL) # Return NULL since the output is saved to a file
}


for (i in 1:500) {  
  set.seed(100 + i)
  indices <- sample(nrow(btsp_input_ds), replace = TRUE)
  boot.func(data = btsp_input_ds, indices = indices, file_name = "bootstrapFullAdjW975.xlsx")
}  

file_path <- "bootstrapFullAdjW995.xlsx"

sheet_names <- excel_sheets(file_path)

all_data <- lapply(seq_along(sheet_names), function(i) {
  data <- read_excel(file_path, sheet = sheet_names[i])
  data <- data %>% mutate(NBootstrap = i)
  return(data)
})

final_data <- bind_rows(all_data)


ci_results <- final_data %>%
  group_by(time) %>%
  summarise(
    surv095CILow = calculate_ci(surv0)[1],
    surv095CIUp = calculate_ci(surv0)[2],
    risk095CILow = calculate_ci(risk0)[1],
    risk095CIUp = calculate_ci(risk0)[2],
    surv195CILow = calculate_ci(surv1)[1],
    surv195CIUp = calculate_ci(surv1)[2],
    risk195CILow = calculate_ci(risk1)[1],
    risk195CIUp = calculate_ci(risk1)[2],
    survdiff95CILow = calculate_ci(survdiff)[1],
    survdiff95CIUp = calculate_ci(survdiff)[2],
    RiskRatio95CILow = calculate_ci(riskratio)[1],
    RiskRatio95CIUp = calculate_ci(riskratio)[2]
  ) %>%
  ungroup()

bsl.cw.graph2 <- bsl.cw.graph %>%
  left_join(ci_results, by = "time")

bsl.cw.plot<-ggplot(bsl.cw.graph2, aes(x=time_mo)) +
  geom_line(aes(y = surv0, colour = "mAB 1st line")) +
  geom_line(aes(y = surv1, colour = "mAB 2nd line")) +
  # DASHED CIs
  geom_ribbon(aes(ymin = surv095CILow, ymax = surv095CIUp, fill = "mAB 1st line"), alpha = 0.2, linetype = "dashed") +
  geom_ribbon(aes(ymin = surv195CILow, ymax = surv195CIUp, fill = "mAB 2nd line"), alpha = 0.2, linetype = "dashed") +
  xlab("Months") +
  scale_x_continuous(limits = c(0, 48), breaks=seq(0,48,6)) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  ylab("Survival") +
  ggtitle("Parametric Survival Curve, baseline and weights- adjusted: Pooled Trial") +
  labs(colour = "Strategy", fill = "Confidence Interval") +
  theme_bw() +
  theme(legend.position="bottom")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## Extract outputs for all the timepoints of interest ----

time_points <- c(27, 53, 79, 105, 131, 157, 183, 209) 
 
unadj_timepoints <- unadj.graph2 %>%
  filter(time %in% time_points) %>%
  mutate(
    R0_CI = paste0(round(surv0,3), " (", surv095CILow, ";", surv095CIUp, ")"),
    R1_CI = paste0(round(surv1,3), " (", surv195CILow, ";", surv195CIUp, ")"),
    RD_CI = paste0(round(survdiff,3), " (", survdiff95CILow, ";", survdiff95CIUp, ")"),
    RR_CI = paste0(round(riskratio,3), " (", RiskRatio95CILow, ";", RiskRatio95CIUp, ")")
  ) %>%
  select(time, R0_CI, R1_CI,RD_CI, RR_CI)

bsladj_timepoints <- bsl.graph2 %>%
  filter(time %in% time_points) %>%
  mutate(
    R0_CI = paste0(round(surv0,3), " (", surv095CILow, ";", surv095CIUp, ")"),
    R1_CI = paste0(round(surv1,3), " (", surv195CILow, ";", surv195CIUp, ")"),
    RD_CI = paste0(round(survdiff,3), " (", survdiff95CILow, ";", survdiff95CIUp, ")"),
    RR_CI = paste0(round(riskratio,3), " (", RiskRatio95CILow, ";", RiskRatio95CIUp, ")")
  ) %>%
  select(time, R0_CI, R1_CI,RD_CI, RR_CI)


fulladj_timepoints <- bsl.cw.graph2 %>%
  filter(time %in% time_points) %>%
  mutate(
    R0_CI = paste0(round(surv0,3), " (", surv095CILow, ";", surv095CIUp, ")"),
    R1_CI = paste0(round(surv1,3), " (", surv195CILow, ";", surv195CIUp, ")"),
    RD_CI = paste0(round(survdiff,3), " (", survdiff95CILow, ";", survdiff95CIUp, ")"),
    RR_CI = paste0(round(riskratio,3), " (", RiskRatio95CILow, ";", RiskRatio95CIUp, ")")
  ) %>%
  select(time, R0_CI, R1_CI,RD_CI, RR_CI)


# Meta-analysis of sequential trials----

library(meta)

## Unadjusted----
df <- data.frame(
  trial = 1:8,
  hr =       h_u,
  ci_lower = l_u,
  ci_upper = u_u)

df$se <- (log(df$ci_upper) - log(df$ci_lower)) / (2 * 1.96)

result <- metagen(
  TE = log(df$hr),             
  seTE = df$se,
  studlab = paste("Trial", df$trial),
  sm = "HR",                   
  method.tau = "REML",         
  hakn = TRUE                  
)

## Baseline-adjusted ----
df <- data.frame(
  trial = 1:8,
  hr =       h_a,
  ci_lower = l_a,
  ci_upper = u_a)

df$se <- (log(df$ci_upper) - log(df$ci_lower)) / (2 * 1.96)

result <- metagen(
  TE = log(df$hr),             
  seTE = df$se,
  studlab = paste("Trial", df$trial),
  sm = "HR",                   
  method.tau = "REML",         
  hakn = TRUE                  
)

## Fully-adjusted ----

df <- data.frame(
  trial = 1:8,
  hr =       h_ca,
  ci_lower = l_ca,
  ci_upper = u_ca)


df$se <- (log(df$ci_upper) - log(df$ci_lower)) / (2 * 1.96)

result <- metagen(
  TE = log(df$hr),             
  seTE = df$se,
  studlab = paste("Trial", df$trial),
  sm = "HR",                   
  method.tau = "REML",         
  hakn = TRUE                  
)
