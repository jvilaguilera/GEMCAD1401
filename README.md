# About

This repository contains the R scripts developed to conduct the analyses of the GEMCAD 1401 registry data. 
Code was developed by Julia Vila Guilera, Manuel Zamparini, and Xabier Garcia de Albeniz.

**Study Title**:Studying When to Add Biological Therapy to Cytotoxic Therapy in Advanced Cancer Care Using Real World Data (RWD): 
Challenges and Solutions to Avoid Biases that Generate Immortal Time

**Overview**: We conducted a study that aimed to estimate the effect of initiating monoclonal antibody (mAB) therapies along with first 
line chemotherapy vs initiating them with second line chemotherapy in the overall survival of patients diagnosed
with metastatic colorectal cancer (mCRC). To investigate this effect, we used data from 1014 eligible patients from 
the Spanish Multidisciplinary Group on Digestive Cancer (GEMCAD) Registry data. Careful study design choices were needed to 
estimate the effect of interest without introducing immortal time bias due to selection or due to treatment misclassification.  

**Methods**: To avoid design choices that would introduce immortal time, we specified a target trial to estimate the effect 
of of “initiating mAB within 8 weeks (grace period) of starting first line chemotherapy” versus “deferring their addition
 to second line chemotherapy” on the overall survival of mCRC patients. We emulated the trial by aligning eligibility criteria
with time zero (corresponding to the date of initiation of first line chemotherapy) and by classifying patients into treatment 
strategies via two methods: 

- Via cloning, censoring and weighting
- Via the emulation of 8 sequential trials, starting at each of the grace period weeks.

Then, survival curves adjusted for baseline and post-baseline covariates were estimated via weighted pooled logistic regression.

**Analytical scripts**

The script for the clone-censor-weight analysis available in this repository contains the following analytical steps: 
- Loading of required libraries and data (unavailable)
- Specification of study variables
- Cloning of individuals and assignment to treatment strategies
- Censoring of individuals when they deviate from assigned strategy
- Estimation of subject-specific time-varying non-stabilized inverse-probability weights
- Estimation of non-parametric (unadjusted) survival curves
- Estimation of parametric (unadjusted) survival curves
- Estimation of parametric baseline-adjusted survival curves
- Estimation of parametric baseline and time-varying IPW-adjusted survival curves
- Bootstrapping procedures to estimate variance

The script for the sequential trial emulation analysis available in this repository contains the following analytical steps:
- Loading of required libraries and data (unavailable)
- Specification of study variables
- Trials creation and assignment to treatment strategies (censoring individuals when they deviate from assigned strategy)
- Function for updating baseline variables based on the trial number
- Estimation of non-parametric (unadjusted) survival curves in single trials
- Estimation of parametric (unadjusted) survival curves in single trials
- Estimation of parametric baseline-adjusted survival curves in single trials
- Estimation of subject-specific time-varying non-stabilized inverse-probability weights
- Estimation of parametric baseline and time-varying IPW-adjusted survival curves in single trials
- Estimation of all the previous survival curves in the Pooled trial
- Bootstrapping procedures to estimate variance in the Pooled trial
