**Title:** GEMCAD Real-World Data Study. When to Start Targeted Therapies in Metastatic Colorectal Cancer

**Overview:** This repository contains the R scripts developed by Julia Vila Guilera, Manuel Zamparini, and Xabier Garcia de Albeniz for the GEMCAD Real-World Data Study. The GEMCAD Real-World Data Study consists of the analysis of a dataset of metastatic colorectal cancer patients in Spain to estimate the causal effect of initiating targeted therapies at different lines of treatment on their overall survival. To draw causal inference estimates without introducing bias, a clone-censor-weight analysis was implemented. Additionally, a sequential trial emulation analysis was carried out to contrast the results obtained. 

**Data:** The Spanish Multidisciplinary Group on Digestive Cancer (GEMCAD) Registry collects data from colorectal cancer patients across Spain. Data from 1014 eligible patients was used to implement these analyses. 

**Clone-censor-weight analysis:** Using real-world data to quantify the effect of treatment initiation at different times on survival outcomes can introduce bias since only those who live for a longer time can receive treatment at later stages. Cloning-censoring-weighting is an analytical approach that eliminates immortal time bias by: 1-cloning people to assign 1 clone to each different treatment strategy, 2- censoring if and when the clones deviate from their assigned treatment strategy, 3- developing inverse probability weights to adjust for the selection bias introduced by informative censoring. 

**Sequential trial emulation analysis:** An alternative to address immortal time bias due to initiation of treatment at different times is the emulation of a sequence of hypothetical trials, each with a baseline at a different time point. Briefly, sequential trial emulation consists of creating a series of trials starting at each sequential time point. At each time point, patient's eligibility is assessed and they are assigned to the strategy their observed data is consistent with at that time point. Once the patient deviates from their strategy, they are censored.

**Analysis steps:**

- setup: this script loads required packages and imports data
- prep: this script cleans and preprocesses data to prepare the analytical dataset (ds12)

*Clone-censor-weight analysis*
- cloning: this script clones the population, assigns each clone to a treatment strategy, and censors if and when a clone stops adhereing to their assigned strategy
- param.unadj: outputs the unadjusted effect estimates
- param.bsl.adj: outputs the baseline adjusted effect estimates
- weightmodel: estimates the probability of receiving treatment each week conditional on baseline and postbaseline patient characteristics and develops censoring weights for each clone-week. 
- param.full.adj: outputs the baseline and postbaseline censoring weights-adjusted effect estimates

*Sequential trial emulation analysis*
- Trial setup: Considered variations in baseline and maximum follow-up truncation
- Trials Analysis:
  - Assessed time under follow-up
  - Non-parametric survival curves
  - Parametric models with different adjustments:
    - Unadjusted
    - Baseline-adjusted
    - Baseline and censoring weight-adjusted (Weight estimation: Derived from the probability of receiving treatment each week, based on baseline and post-baseline patient characteristics)
  - Risk estimation:
    - Absolute and relative risk
    - Hazard ratios
- Pooled analysis: Repeated the same analysis on pooled trials, combining all sequential trials into one dataset
- Confidence intervals: Estimated using 500 bootstrap replications
- Meta-analysis: Conducted to investigate potential anomalies between trials
