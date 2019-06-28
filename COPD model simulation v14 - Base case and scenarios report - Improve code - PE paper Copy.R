###############################################################
########## BI/iMTA COPD simulation model ######################
###############################################################


###################################
######## LOG version ##############
###################################

### Modificaitions with respect to previous versions

### November 2017 ###

# 1. Allow running the model as a function of several inputs (that will be listed below).
# 2. Implement PSA.
# 3. Remove Rayleigh and Log-logistic from mortality distirbutions.
# 4. Implement treatment effect as an "improvement" controlled by random seed
# 5. Correct mortality factor based on remaining life expectancy
# 6. Add probability of severe eacerbation as treatment effect --> DONT forget to add it to the PSA
# 7. Constraints on fevppa <<-- ask Martine about logical constraints - for the time being add 0

### December 2017 ###

# 1. Reduce the number of if statements
# 2. Remove parts regarding "no_event_time"
# 3. Include treatment effect scenarios as input parameters

### July 2018 ###

# 1. Add new metric for symptoms for the ViH paper. Note this metric is already implmented in the interfaced version 
#    of the model.


#######################################
########## Code starts below ##########
#######################################


### Remove previous objects from R session
rm(list = ls())  


init_time <- proc.time()

### Install and load packages - Update this: ask Frederick!
library(lattice)
library(MASS)
library(survival)
library(plyr)
library(triangle)


### Set working directory: can this be automathized?
wd <- setwd("C:/Users/Isaac/Dropbox/COPD")

##########################
### AUXILIAR FUNCTIONS ###
##########################
 
### Read auxiliar function from a different file
source("COPD model simulation - auxiliar functions.R")


########################
### SIMULATION LOGIC ###
########################

### Simulation parameters
### 1. run_obs_input = 1 uses regression equations based on observed data. Otherwise, usees regression equations based on observed+predicted data
### 2. patient_size_input = number of patients included in the simulation.
### 3. run_PSA_input = runs the model in probabilistic mode. Otherwise, deterministic.
### 4. exac_treatment_effect_tte_input = variable to indicate increase (or decrease) in time to exacerbation. Default should be 1.
### 5. fev1_treatment_effect_input = variable to indicate a reduction or increase in lung function decline
### 6. subgroup_input = variable string to indicate which set of patients should be used in the analysis. This should be linked
###    to a dropdown list or similar. Default should be all patients (base-case) so no subgorup.


COPD_model_simulation <- function(run_obs_input,
                                  patient_size_input,
                                  run_PSA_input,
                                  exac_treatment_effect_tte_input,
                                  exac_treatment_effect_sevexaprob_input,
                                  fev1_treatment_effect_input,
                                  cwe_treatment_effect_input,
                                  sgact_treatment_effect_input,
                                  coughsputum_treatment_effect_input,
                                  breathless_treatment_effect_input,
                                  sgtot_treatment_effect_input,
                                  subgroup_input,
                                  seed_input){
  
  
  
  
  ### Step 1: Read regression coefficients depending on the data used to predict them.
  ###         since the code is reading from csv files in certain folders, it is important to name these
  ###         folders as suggested by the code.
  
  if(run_obs_input==1){
    
    lung_function_regression_coef         <- (read.csv(paste0(wd,c("/Model - regression coefficients/Lung function/lung_function_regression_coef_observed_data2.csv")),sep=";"))$Value
    lung_function_cov_matrix              <- read.csv(paste0(wd,c("/Model - regression coefficients/Lung function/lung_function_cov_matrix_observed_data2.csv")),sep=";")
    
    cwe_tot_regression_coef               <- (read.csv(paste0(wd,c("/Model - regression coefficients/Exercise capacity/cwe_tot_regression_coef_observed_data_v3.csv")),sep=","))$Value
    cwe_tot_cov_matrix                    <- read.csv(paste0(wd,c("/Model - regression coefficients/Exercise capacity/cwe_tot_cov_matrix_observed_data_v3.csv")),sep=";")
    
    breathless_regression_coef            <- (read.csv(paste0(wd,c("/Model - regression coefficients/Symptoms/breathless_regression_coef_observed_data_v2.csv")),sep=";"))$Estimate
    breathless_cov_matrix                 <- read.csv(paste0(wd,c("/Model - regression coefficients/Symptoms/breathless_cov_matrix_observed_data_v2.csv")),sep=";")
    
    coughsputum_regression_coef           <- (read.csv(paste0(wd,c("/Model - regression coefficients/Symptoms/coughsputum_regression_coef_observed_data_v2.csv")),sep=";"))$Estimate
    coughsputum_cov_matrix                <- read.csv(paste0(wd,c("/Model - regression coefficients/Symptoms/coughsputum_cov_matrix_observed_data_v2.csv")),sep=",")
    
    SGACT_regression_coef                 <- (read.csv(paste0(wd,c("/Model - regression coefficients/Physical activity/SGACT_regression_coef_observed_data_v2.csv")),sep=";"))$Value
    SGACT_cov_matrix                      <- read.csv(paste0(wd,c("/Model - regression coefficients/Physical activity/SGACT_cov_matrix_observed_data_v2.csv")),sep=",")
    
    SGTOT_regression_coef                 <- (read.csv(paste0(wd,c("/Model - regression coefficients/Quality of life/SGTOT_regression_coef_observed_data_v2.csv")),sep=";"))$Value
    SGTOT_cov_matrix                      <- read.csv(paste0(wd,c("/Model - regression coefficients/Quality of life/SGTOT_cov_matrix_observed_data_v2.csv")),sep=",")
    
    mortality_weibull_regression_coef     <- (read.csv(paste0(wd,c("/Model - regression coefficients/Mortality/mortality_weibull_regression_coef_observed_data.csv")),sep=";"))$Value
    mortality_weibull_cov_matrix          <- read.csv(paste0(wd,c("/Model - regression coefficients/Mortality/mortality_weibull_cov_matrix_observed_data.csv")),sep=";")
    
    exacerbation_weibull_regression_coef  <- (read.csv(paste0(wd,c("/Model - regression coefficients/Exacerbations/exacerbation_weibull_regression_coef_observed_data.csv")),sep=";"))$Value
    exacerbation_weibull_cov_matrix       <- read.csv(paste0(wd,c("/Model - regression coefficients/Exacerbations/exacerbation_weibull_cov_matrix_observed_data.csv")),sep=",")
    
    exacerbation_severity_regression_coef <- (read.csv(paste0(wd,c("/Model - regression coefficients/Exacerbations/exacerbation_severity_regression_coef_observed_data.csv")),sep=";"))$Estimate
    exacerbation_severity_cov_matrix      <- read.csv(paste0(wd,c("/Model - regression coefficients/Exacerbations/exacerbation_severity_cov_matrix_observed_data.csv")),sep=";")
    
    pneumonia_weibull_regression_coef     <- (read.csv(paste0(wd,c("/Model - regression coefficients/Pneumonia/pneu_weibull_regression_coef_observed_data.csv")),sep=";"))$Value
    pneumonia_weibull_cov_matrix          <- read.csv(paste0(wd,c("/Model - regression coefficients/Pneumonia/pneu_weibull_cov_matrix_observed_data.csv")),sep=",")
    
    pneumonia_hosp_regression_coef        <- (read.csv(paste0(wd,c("/Model - regression coefficients/Pneumonia/pneu_hosp_regression_coef_observed_data.csv")),sep=";"))$Estimate
    pneumonia_hosp_cov_matrix             <- read.csv(paste0(wd,c("/Model - regression coefficients/Pneumonia/pneu_hosp_cov_matrix_observed_data.csv")),sep=",")
    
    
    ### Step 2: assign predictors
    
    ### MIND THE PREDICTORS: WHEN DATA IS OBSERVED WE DONT USED CWE_TOT
    
    ### Note: I thouhgt of making the lists below more efficient but the order is important. This could be an improvement for the future.
    
    fev1_predictors <- c("ANLYEAR","female","AGE","BMIclass2","BMIclass3","SMOKER","SMPKY","other_CVD"
                         ,"Reversibility","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS"
                         ,"FEVA_BL","modexac_yn","sevexac_yn")
    
    cwe_tot_predictors <- c("female","BMIclass2","BMIclass3","SMOKER","SMPKY","other_CVD"
                            ,"MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS"
                            ,"age_time","fevppa","lag_SGACT","lag_CWE_TOT","modexac_yn","sevexac_yn")
    
    breathless_predictors <- c("female","BMIclass2","BMIclass3","SMOKER","SMPKY_SCALED","other_CVD"
                               ,"Reversibility_SCALED","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS"
                               ,"ANLYEAR_SCALED","AGE_SCALED", "fevppa_SCALED","modexac_yn","sevexac_yn",
                               "SGACT_SCALED","lag_breathlessyn")
    
    coughsputum_predictors <- c("female","BMIclass2","BMIclass3","SMOKER","SMPKY_SCALED","other_CVD"
                                ,"Reversibility_SCALED","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS"
                                ,"ANLYEAR_SCALED","AGE_SCALED","fevppa_SCALED","modexac_yn","sevexac_yn"
                                ,"SGACT_SCALED","lag_coughsputumyn")
    
    SGACT_predictors <- c("female","BMIclass2","BMIclass3","SMOKER","SMPKY","other_CVD"
                          ,"Reversibility","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS"
                          ,"ANLYEAR","AGE","fevppa","lag_SGACT", "lag_breathlessyn"
                          ,"lag_coughsputumyn","lag_SGTOT")
    
    SGTOT_predictors <- c("female","BMIclass2","BMIclass3","SMOKER","SMPKY","other_CVD"
                          ,"Reversibility","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS"
                          ,"ANLYEAR","AGE","fevppa","lag_SGTOT","modexac_yn","sevexac_yn","SGACT"
                          ,"breathlessyn","coughsputumyn","pneu_yn") #No adverse events at baseline so all =0
    
    mortality_predictors <- c("female","AGE","BMIclass2","BMIclass3","SMOKER","SMPKY","other_CVD"
                              ,"Reversibility","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS"
                              ,"fevppa","prevsevexacyn","SGACT"
                              ,"breathlessyn","coughsputumyn","SGTOT")
    
    exacerbation_predictors <- c("female","BMIclass2","BMIclass3","SMOKER","SMPKY","other_CVD"
                                 ,"Reversibility","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS"
                                 ,"lag_fevppa","lag_SGACT","prevtotexacyn","prevsevexacyn"
                                 ,"lag_SGTOT","age_time")
    
    exacerbation_severity_predictors <- c("female","BMIclass2","BMIclass3","SMOKER","SMPKY_SCALED","other_CVD"
                                          ,"Reversibility_SCALED","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS"
                                          , "lag_fevppa_SCALED","lag_SGACT_SCALED","prevtotexacyn", "prevsevexacyn"
                                          , "lag_SGTOT_SCALED","age_time_SCALED")
    
    pneumonia_predictors <- c("female", "BMIclass2","BMIclass3","SMOKER","SMPKY","other_CVD"
                              ,"Reversibility","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS"
                              ,"age_time")
    
    pneumonia_hosp_predictors <- c("female", "BMIclass2","BMIclass3","SMOKER","SMPKY_SCALED","other_CVD"
                                   ,"Reversibility_SCALED","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS"
                                   ,"age_time_SCALED")
    
    ### Step 3: Read the data with the patient characteristics and select the complete cases (observations with all patient characteristics)
    baseline_characteristics_run <- read.csv(paste0(wd,c("/Model - datasets/baseline_characteristics_observed_data.csv")),sep=",")
    complete_cases <- baseline_characteristics_run[colnames(baseline_characteristics_run)[-c(17,18,51)]][complete.cases(baseline_characteristics_run[colnames(baseline_characteristics_run)[-c(17,18,51)]]),]
    
    
    ### Step 4: indicate the patient characteristics that we will save during the simulation. Note CWE_TOT is not included 
    ###         as predictor when the model is run in "observed" data mode. 
    history_characteristics <- c("SIMID","PTID","ANLYEAR","age_time","FEVA","fevppa","sevexac_yn","modexac_yn","SGACT","SGTOT","coughsputumyn","breathlessyn","pneu_yn","pneu_hosp_yn","dead",
                                   "female","AGE","BMIclass2","BMIclass3","SMOKER","SMPKY","other_CVD",
                                   "Reversibility", "MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS","FEVA_BL","HTSTD",
                                   "lag_SGACT","lag_breathlessyn","lag_coughsputumyn","lag_SGTOT",
                                   "SMPKY_SCALED","Reversibility_SCALED","ANLYEAR_SCALED","AGE_SCALED","fevppa_SCALED","SGACT_SCALED")
    
  }else{
    
    ### Step 1: read regression coefficients
    
    lung_function_regression_coef         <- (read.csv(paste0(wd,c("/Model - regression coefficients/Lung function/lung_function_regression_coef_predicted_data2.csv")),sep=","))$Value
    lung_function_cov_matrix              <- read.csv(paste0(wd,c("/Model - regression coefficients/Lung function/lung_function_cov_matrix_predicted_data2.csv")),sep=";")
    
    cwe_tot_regression_coef               <- (read.csv(paste0(wd,c("/Model - regression coefficients/Exercise capacity/cwe_tot_regression_coef_predicted_data_v3.csv")),sep=","))$Value
    cwe_tot_cov_matrix                    <- read.csv(paste0(wd,c("/Model - regression coefficients/Exercise capacity/cwe_tot_cov_matrix_predicted_data_v3.csv")),sep=";")
    
    breathless_regression_coef            <- (read.csv(paste0(wd,c("/Model - regression coefficients/Symptoms/breathless_regression_coef_predicted_data_v2.csv")),sep=";"))$Estimate
    breathless_cov_matrix                 <- read.csv(paste0(wd,c("/Model - regression coefficients/Symptoms/breathless_cov_matrix_predicted_data_v2.csv")),sep=";")
    
    coughsputum_regression_coef           <- (read.csv(paste0(wd,c("/Model - regression coefficients/Symptoms/coughsputum_regression_coef_predicted_data_v2.csv")),sep=";"))$Estimate
    coughsputum_cov_matrix                <- read.csv(paste0(wd,c("/Model - regression coefficients/Symptoms/coughsputum_cov_matrix_predicted_data_v2.csv")),sep=",")
    
    SGACT_regression_coef                 <- (read.csv(paste0(wd,c("/Model - regression coefficients/Physical activity/SGACT_regression_coef_predicted_data_v2.csv")),sep=";"))$Value
    SGACT_cov_matrix                      <- read.csv(paste0(wd,c("/Model - regression coefficients/Physical activity/SGACT_cov_matrix_predicted_data_v2.csv")),sep=",")
    
    SGTOT_regression_coef                 <- (read.csv(paste0(wd,c("/Model - regression coefficients/Quality of life/SGTOT_regression_coef_predicted_data_v2.csv")),sep=";"))$Value
    SGTOT_cov_matrix                      <- read.csv(paste0(wd,c("/Model - regression coefficients/Quality of life/SGTOT_cov_matrix_predicted_data_v2.csv")),sep=",")
    
    mortality_weibull_regression_coef     <- (read.csv(paste0(wd,c("/Model - regression coefficients/Mortality/mortality_weibull_regression_coef_predicted_data.csv")),sep=";"))$Value
    mortality_weibull_cov_matrix          <- read.csv(paste0(wd,c("/Model - regression coefficients/Mortality/mortality_weibull_cov_matrix_predicted_data.csv")),sep=",")
    
    exacerbation_weibull_regression_coef  <- (read.csv(paste0(wd,c("/Model - regression coefficients/Exacerbations/exacerbation_weibull_regression_coef_predicted_data.csv")),sep=";"))$Value
    exacerbation_weibull_cov_matrix       <- read.csv(paste0(wd,c("/Model - regression coefficients/Exacerbations/exacerbation_weibull_cov_matrix_predicted_data.csv")),sep=",")
    
    exacerbation_severity_regression_coef <- (read.csv(paste0(wd,c("/Model - regression coefficients/Exacerbations/exacerbation_severity_regression_coef_predicted_data.csv")),sep=";"))$Estimate
    exacerbation_severity_cov_matrix      <- read.csv(paste0(wd,c("/Model - regression coefficients/Exacerbations/exacerbation_severity_cov_matrix_predicted_data.csv")),sep=",")
    
    pneumonia_weibull_regression_coef     <- (read.csv(paste0(wd,c("/Model - regression coefficients/Pneumonia/pneu_weibull_regression_coef_predicted_data.csv")),sep=";"))$Value
    pneumonia_weibull_cov_matrix          <- read.csv(paste0(wd,c("/Model - regression coefficients/Pneumonia/pneu_weibull_cov_matrix_predicted_data.csv")),sep=",")
    
    pneumonia_hosp_regression_coef        <- (read.csv(paste0(wd,c("/Model - regression coefficients/Pneumonia/pneu_hosp_regression_coef_predicted_data.csv")),sep=";"))$Estimate
    pneumonia_hosp_cov_matrix             <- read.csv(paste0(wd,c("/Model - regression coefficients/Pneumonia/pneu_hosp_cov_matrix_predicted_data.csv")),sep=",")
    
    ### Step 2: assign predictors
    
    fev1_predictors       <- c("ANLYEAR","female","AGE","BMIclass2","BMIclass3","SMOKER","SMPKY","other_CVD",
                                "Reversibility","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","EOS_yn","ICS",
                                "FEVA_BL","modexac_yn","sevexac_yn")
    
    cwe_tot_predictors    <- c("female","BMIclass2","BMIclass3","SMOKER","SMPKY","other_CVD",
                               "Reversibility","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS", "EOS_yn",
                               "age_time","fevppa","lag_SGACT","lag_CWE_TOT","modexac_yn","sevexac_yn")
    
    breathless_predictors <- c("female","BMIclass2","BMIclass3","SMOKER","SMPKY_SCALED","other_CVD",
                               "Reversibility_SCALED","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS","EOS_yn",
                               "ANLYEAR_SCALED","AGE_SCALED","fevppa_SCALED","modexac_yn","sevexac_yn","SGACT_SCALED","CWE_TOT_SCALED","lag_breathlessyn") 
    
    coughsputum_predictors <- c("female","BMIclass2","BMIclass3","SMOKER","SMPKY_SCALED","other_CVD",
                                "Reversibility_SCALED","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS", "EOS_yn",
                                "ANLYEAR_SCALED","AGE_SCALED","fevppa_SCALED","modexac_yn","sevexac_yn","SGACT_SCALED","CWE_TOT_SCALED","lag_coughsputumyn")
    
    SGACT_predictors       <- c("female","BMIclass2","BMIclass3","SMOKER","SMPKY","other_CVD",
                                "Reversibility","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS","EOS_yn",
                                "ANLYEAR", "AGE","fevppa","lag_SGACT","CWE_TOT","lag_breathlessyn","lag_coughsputumyn","lag_SGTOT")
    
    SGTOT_predictors       <- c("female","BMIclass2","BMIclass3","SMOKER","SMPKY","other_CVD",
                                "Reversibility","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS","EOS_yn",
                                "ANLYEAR", "AGE","fevppa","lag_SGTOT","modexac_yn","sevexac_yn","SGACT","CWE_TOT",
                                "breathlessyn","coughsputumyn","pneu_yn") #No adverse events at baseline so all =0
    
    mortality_predictors   <- c("female","AGE","BMIclass2","BMIclass3","SMOKER","SMPKY","other_CVD",
                                "Reversibility","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS","EOS_yn",
                                "fevppa","prevsevexacyn","SGACT","CWE_TOT","breathlessyn","coughsputumyn","SGTOT")
    
    exacerbation_predictors <- c("female","BMIclass2","BMIclass3","SMOKER","SMPKY","other_CVD",
                                 "Reversibility","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS","EOS_yn",
                                 "lag_fevppa","lag_SGACT","prevtotexacyn","prevsevexacyn","lag_SGTOT","age_time")
    
    exacerbation_severity_predictors <- c("female","BMIclass2","BMIclass3","SMOKER","SMPKY_SCALED","other_CVD",
                                          "Reversibility_SCALED","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS","EOS_yn",
                                          "lag_fevppa_SCALED","lag_SGACT_SCALED","prevtotexacyn", "prevsevexacyn", "lag_SGTOT_SCALED","age_time_SCALED")
    
    pneumonia_predictors <- c("female","BMIclass2","BMIclass3","SMOKER","SMPKY","other_CVD",
                              "Reversibility","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS","EOS_yn","age_time")
    
    pneumonia_hosp_predictors <- c("female", "BMIclass2","BMIclass3","SMOKER","SMPKY_SCALED","other_CVD",
                                   "Reversibility_SCALED","MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS","EOS_yn","age_time_SCALED")
    
    ### Step 3: Read the data with the patient characteristics and select the complete cases
    baseline_characteristics_run <- read.csv(paste0(wd,c("/Model - datasets/baseline_characteristics_predicted_data.csv")),sep=",")
    complete_cases <- baseline_characteristics_run[colnames(baseline_characteristics_run)][complete.cases(baseline_characteristics_run[colnames(baseline_characteristics_run)]),]
    
    ### Step 4: indicate the patient characteristics that we will save during the simulation. Note CWE_TOT is not included 
    ###         as predictor when the model is run in "observed" data mode. 
    history_characteristics <- c("SIMID","PTID","ANLYEAR","age_time","FEVA","fevppa","sevexac_yn","modexac_yn","CWE_TOT","SGACT","SGTOT","coughsputumyn","breathlessyn","pneu_yn","pneu_hosp_yn","dead",
                                 "female","AGE","BMIclass2","BMIclass3","SMOKER","SMPKY","other_CVD",
                                 "Reversibility", "MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","EOS_yn","ICS","FEVA_BL","HTSTD",
                                 "lag_SGACT","lag_CWE_TOT","lag_breathlessyn","lag_coughsputumyn","lag_SGTOT",
                                 "SMPKY_SCALED","Reversibility_SCALED","ANLYEAR_SCALED","AGE_SCALED","fevppa_SCALED","SGACT_SCALED","CWE_TOT_SCALED")
    
  }
  
  
  ### The line below was used to run the model for the "average" patient mode. Not sure if we will do this in the future.
  ### baseline_characteristics_run_average <- data.frame(t(colMeans(baseline_characteristics_run[,-(1:2)],na.rm=TRUE)))
  
  
  ### Step 5: select subgroups 
  
  complete_cases <- switch(subgroup_input,
                           "low BMI"               = complete_cases[which(complete_cases$BMIclass==1),],
                           "high BMI"              = complete_cases[which(complete_cases$BMIclass==3),],
                           "frequent exacerbators" = complete_cases[which(complete_cases$prevsevexacyn==1),],
                           "emphysema"             = complete_cases[which(complete_cases$EMPHDIA==1),],
                           "Dynagito"              = complete_cases[which(complete_cases$prevtotexacyn==1 & complete_cases$fevppa<60),],
                           complete_cases) # last item is default = all population
  

  ##################################################
  ########## MAIN PART I: simulate events ##########
  ##################################################
  
  
  ### Regression coefficients for PSA -- once per PSA
  
  if(run_PSA_input == 1){
    
    lung_function_regression_coef <- mvrnorm(1,lung_function_regression_coef,lung_function_cov_matrix[,-1])
    
    if(run_obs_input==0){
      cwe_tot_regression_coef <- mvrnorm(1,cwe_tot_regression_coef,cwe_tot_cov_matrix[,-1])
    }
    
    SGACT_regression_coef                 <- mvrnorm(1,SGACT_regression_coef,SGACT_cov_matrix[,-1])
    SGTOT_regression_coef                 <- mvrnorm(1,SGTOT_regression_coef,SGTOT_cov_matrix[,-1])
    exacerbation_weibull_regression_coef  <- mvrnorm(1,exacerbation_weibull_regression_coef,exacerbation_weibull_cov_matrix[,-1])
    pneumonia_weibull_regression_coef     <- mvrnorm(1,pneumonia_weibull_regression_coef,pneumonia_weibull_cov_matrix[,-1])
    exacerbation_severity_regression_coef <- mvrnorm(1,exacerbation_severity_regression_coef,exacerbation_severity_cov_matrix[,-1])
    pneumonia_hosp_regression_coef        <- mvrnorm(1,pneumonia_hosp_regression_coef,pneumonia_hosp_cov_matrix[,-1])
    breathless_regression_coef            <- mvrnorm(1,breathless_regression_coef,breathless_cov_matrix[,-1])
    coughsputum_regression_coef           <- mvrnorm(1,coughsputum_regression_coef,coughsputum_cov_matrix[,-1])
    mortality_weibull_regression_coef     <- mvrnorm(1,mortality_weibull_regression_coef,mortality_weibull_cov_matrix[,-1]) 
    
    ### Treatment effect parameters must come at the end to ensure the same random seed in treatment arm  
    ### Consider changing the names of thees 3 calculated variables
    exac_treatment_effect_tte_input        <- runif(1,min(c(exac_treatment_effect_tte_input-(exac_treatment_effect_tte_input-1)*0.1,exac_treatment_effect_tte_input+(exac_treatment_effect_tte_input-1)*0.1)),max(c(exac_treatment_effect_tte_input-(exac_treatment_effect_tte_input-1)*0.1,exac_treatment_effect_tte_input+(exac_treatment_effect_tte_input-1)*0.1)))
    exac_treatment_effect_sevexaprob_input <- runif(1,min(c(exac_treatment_effect_sevexaprob_input-(exac_treatment_effect_sevexaprob_input-1)*0.1,exac_treatment_effect_sevexaprob_input+(exac_treatment_effect_sevexaprob_input-1)*0.1)),max(c(exac_treatment_effect_sevexaprob_input-(exac_treatment_effect_sevexaprob_input-1)*0.1,exac_treatment_effect_sevexaprob_input+(exac_treatment_effect_sevexaprob_input-1)*0.1)))
    fev1_treatment_effect_input            <- runif(1,min(c(fev1_treatment_effect_input-(fev1_treatment_effect_input-1)*0.1,fev1_treatment_effect_input+(fev1_treatment_effect_input-1)*0.1)),max(c(fev1_treatment_effect_input-(fev1_treatment_effect_input-1)*0.1,fev1_treatment_effect_input+(fev1_treatment_effect_input-1)*0.1)))
    
  } # end if regression coef PSA
  
  
  ### Step 6: sample with replacement from the pool of patients
  
  ### IMPORTANT: Set a random seed to be able to replicate the results.
  ### This first seed is used to draw the same pool of patients when the function is called multiple times.
  ### Decide later if this has to go inside or outside the function. It depends on whether this seed is going to be 
  ### always fixed (inside) or if it can change (outside - and make it random or user defined)
  
  set.seed(seed_input) 
  
  ### Select the sample used in the simulation: sample with replacement from the dataset
  #simulation_baseline_patients <-  complete_cases ##Average
  simulation_baseline_patients <-  complete_cases[sample(nrow(complete_cases), patient_size_input, replace = TRUE), ]
  
  ### Add a simulation ID variable to summarize results after the simulation is finished.
  ### This is needed because patients might be repeated in the simulation, and then aggregating 
  ### results by patient ID is not correct.
  SIMID <- rep("NA",patient_size_input)
  simulation_baseline_patients <- cbind(simulation_baseline_patients,SIMID)
  
  ### Create the simulation patient history table (for now is just empty)
  simulation_patients_history <- simulation_baseline_patients[FALSE,c(history_characteristics)]
  
  ### Choose the 1st patient (later make a loop)
  patient_index <- 1
  
  

  
  
  
  ### Begin the loop on the simulation size (i.e. the number of patients you want to simulate)
  for(patient_index in 1:patient_size_input){
    
    # Print the patient index to know how advanced is the simulation.
    # Try to show this in the interface.
    if(run_PSA_input == 0){print(patient_index)}

    # Pick the current patient from those selected fomr baseline
    current_patient <- simulation_baseline_patients[patient_index,]
    
    #####################################################################################
    # STEP1 : Predict the continuous variables at baseline since observed != predicted. #
    # Otherwise, we get "strange" values after the first event.                         #
    #####################################################################################
    
    current_patient$SIMID <- patient_index
    
    # Baseline predicted FEV1
    current_patient$FEVA          <- max(0,predicted_fev1(run_obs_input,lung_function_regression_coef,current_patient[fev1_predictors],fev1_treatment_effect_input)$fev_1) 
    baseline_fevppa_calc          <- fevppa_calc(current_patient$female,current_patient$HTSTD,current_patient$age_time,current_patient$FEVA)
    current_patient$fevppa        <- baseline_fevppa_calc$fevppa
    current_patient$FEV1pred      <- baseline_fevppa_calc$FEV1_pred
    current_patient$fevppa_SCALED <- (current_patient$fevppa - attr(scale(baseline_characteristics_run$fevppa),"scaled:center"))/attr(scale(baseline_characteristics_run$fevppa),"scaled:scale")
    
    
    # Baseline predicted CWE (min observed is 42. here for now we truncate at 0)
    # Note that this is not used then the model is run in observed data mode.
    if(run_obs_input==0){
      current_patient$CWE_TOT        <- cwe_treatment_effect_input*max(0,predicted_cwe_tot(cwe_tot_regression_coef,current_patient[cwe_tot_predictors])$cwe_tot) #max(0,predicted_cwe_tot(cwe_tot_regression_coef$Value,current_patient[cwe_tot_predictors])$cwe_tot)
      current_patient$CWE_TOT_SCALED <- (current_patient$CWE_TOT - attr(scale(baseline_characteristics_run$CWE_TOT),"scaled:center"))/attr(scale(baseline_characteristics_run$CWE_TOT),"scaled:scale")
    }
    
    # Baseline predicted SGACT
    current_patient$SGACT        <- min(100,max(0,sgact_treatment_effect_input+predicted_SGACT(SGACT_regression_coef,current_patient[SGACT_predictors])$SGACT))
    current_patient$SGACT_SCALED <- (current_patient$SGACT - attr(scale(baseline_characteristics_run$SGACT),"scaled:center"))/attr(scale(baseline_characteristics_run$SGACT),"scaled:scale")
    
    
    # Baseline predicted breathlesness and coughsputum: NOT CHANGED because these are not continuous variables
    
    # Baseline predicted SGTOT
    current_patient$SGTOT        <- min(100,max(0,sgtot_treatment_effect_input+predicted_SGTOT(SGTOT_regression_coef,current_patient[SGTOT_predictors])$SGTOT))
    current_patient$SGTOT_SCALED <- (current_patient$SGTOT - attr(scale(baseline_characteristics_run$SGTOT),"scaled:center"))/attr(scale(baseline_characteristics_run$SGTOT),"scaled:scale")
    
    
    # Save the characteristics to be used in the simulation history (not those that are stable, only those changing) 
    simulation_patients_history <- rbind(simulation_patients_history,current_patient[history_characteristics])
    
    
    #############################################
    # STEP2: Sample from mortality distribution #
    #############################################
    
    ### This seed will be used to draw the life expectancy
    if(run_PSA_input == 0){set.seed(current_patient$SIMID)}else{set.seed((seed_input*100)+current_patient$SIMID)} # loop size is 500x100 so 100 is number of patients
    
    # Mortality at baseline: Weibull
    baseline_remaining_life_exp_parameters <- predicted_mortality_weibull(mortality_weibull_regression_coef,current_patient[mortality_predictors])
    baseline_remaining_life_exp            <- rweibull(1,baseline_remaining_life_exp_parameters$shape_mortality_weibull,baseline_remaining_life_exp_parameters$scale_mortality_weibull)/365
    baseline_remaining_life_exp_mean       <- baseline_remaining_life_exp_parameters$scale_mortality_weibull*gamma(1+1/baseline_remaining_life_exp_parameters$shape_mortality_weibull)/365
    
    #print(baseline_remaining_life_exp)
    
    # Save the sampled life expectancy. This will be used as reference.
    current_remaining_life_exp  <- baseline_remaining_life_exp
    lag_current_mortality       <- baseline_remaining_life_exp_mean 
    
    
    
    ###############################################################
    # STEP 3: Start the "timed" simulation (while loop = clock)   #
    ###############################################################
    
    current_event   <- 1
    factor_for_seed <- 100 #test = 1 to reporduce results. I normally used 100 but can be anything large enough
    
    while(current_patient$dead==0){
      
      ### These seeds will be used to draw the time to events while the patient is alive 
      if(run_PSA_input == 0){
        set.seed(factor_for_seed*current_patient$SIMID + current_event)
      }else{
          set.seed(factor_for_seed*seed_input*current_patient$SIMID + current_event)
        }
      
      #####################################################
      # STEP 3.1: Sample exacerbation and pneumonia time  #
      #####################################################
      
      # Exacerbation is Weibull
      current_exacerbation_parameters <- predicted_exacerbation_weibull(exacerbation_weibull_regression_coef,current_patient[exacerbation_predictors])
      current_exacerbation_time       <- rweibull(1,current_exacerbation_parameters$shape_exacerbation_weibull,current_exacerbation_parameters$scale_exacerbation_weibull)
      current_exacerbation_time       <- exac_treatment_effect_tte_input*current_exacerbation_time
      
      # pneumonia is Weibull
      current_pneumonia_parameters <- predicted_pneumonia_weibull(pneumonia_weibull_regression_coef,current_patient[pneumonia_predictors])
      current_pneumonia_time       <- rweibull(1,current_pneumonia_parameters$shape_pneumonia_weibull,current_pneumonia_parameters$scale_pneumonia_weibull)
      
      ##############################################
      # STEP 3.2: Update patient characteristics   #
      ##############################################
      
      # We first copy all the previous characteristics
      current_patient_update <- current_patient
      
      #######################################################
      # STEP 3.2.1: Update first the lagged characteristics #
      #######################################################
      
      current_patient_update$lag_fevppa        <- current_patient$fevppa
      current_patient_update$lag_fevppa_SCALED <- (current_patient$fevppa - attr(scale(baseline_characteristics_run$fevppa),"scaled:center"))/attr(scale(baseline_characteristics_run$fevppa),"scaled:scale")
      
      current_patient_update$lag_SGACT         <- current_patient$SGACT
      current_patient_update$lag_SGACT_SCALED  <- (current_patient$SGACT - attr(scale(baseline_characteristics_run$SGACT),"scaled:center"))/attr(scale(baseline_characteristics_run$SGACT),"scaled:scale")
      
      current_patient_update$lag_CWE_TOT       <- current_patient$CWE_TOT
      current_patient_update$lag_breathlessyn  <- current_patient$breathlessyn
      current_patient_update$lag_coughsputumyn <- current_patient$coughsputumyn
      
      current_patient_update$lag_SGTOT         <- current_patient$SGTOT
      current_patient_update$lag_SGTOT_SCALED  <- (current_patient$SGTOT - attr(scale(baseline_characteristics_run$SGTOT),"scaled:center"))/attr(scale(baseline_characteristics_run$SGTOT),"scaled:scale")
      
      current_patient_update$prevsevexacyn <- current_patient$sevexac_yn
      current_patient_update$prevtotexacyn <- ifelse(current_patient$sevexac_yn==1 | current_patient$modexac_yn==1,1,0)
      
      
      #######################################################
      # STEP 3.2.2: Update depending on the event occurred  #
      #######################################################
      
      # If exacerbation happened, then update at exacerbation time
      if(min(current_remaining_life_exp,current_exacerbation_time,current_pneumonia_time)==current_exacerbation_time){
        
        # Patient is still alive
        current_patient_update$dead <- 0
        
        # Update time-relatedcharacteristics
        current_patient_update$ANLYEAR         <- current_patient$ANLYEAR  + current_exacerbation_time
        current_patient_update$ANLYEAR_SCALED  <- (current_patient_update$ANLYEAR - if(run_obs_input==1){0.962629}else{1.529411})/if(run_obs_input==1){0.9436266}else{1.376549} # hardcoded based on dataset
        current_patient_update$lag_ANLYEAR     <- current_exacerbation_time
        current_patient_update$age_time        <- current_patient$age_time + current_exacerbation_time
        current_patient_update$age_time_SCALED <- (current_patient_update$age_time - attr(scale(baseline_characteristics_run$age_time),"scaled:center"))/attr(scale(baseline_characteristics_run$age_time),"scaled:scale")
        
        # If age + current life expectancy is > 100 years then we force death
        if(current_patient_update$age_time>100){current_patient_update$dead <- 1}
        
        # Update exacerbation status. Decide first whether the exacerbation is moderate or severe.
        
        current_severity_prob   <- exac_treatment_effect_sevexaprob_input*predicted_exacerbation_severity(exacerbation_severity_regression_coef,current_patient_update[exacerbation_severity_predictors])$p.exacerbation_severity
        current_severity_sample <- rbinom(1,1,current_severity_prob)
        if(current_severity_sample==1){
          current_patient_update$modexac_yn <- 0
          current_patient_update$sevexac_yn <- 1
          
          # Additional death risk because of severe exacerbation 
          current_sevexa_death_prob   <- 0.063 ### hardcoded for now!
          current_sevexa_death_sample <- rbinom(1,1,current_sevexa_death_prob)
          
          if(current_sevexa_death_sample==1){
            current_patient_update$dead <- 1
          }else{
            current_patient_update$dead <- 0
          }
          
          
        }else{
          current_patient_update$modexac_yn <- 1
          current_patient_update$sevexac_yn <- 0
        }
        
        # Update pneumonia status 
        current_patient_update$pneu_yn      <- 0
        current_patient_update$pneu_hosp_yn <- 0
        
      }
      
      
      # If pneumonia happened, then update at pneumonia time
      if(min(current_remaining_life_exp,current_exacerbation_time,current_pneumonia_time)==current_pneumonia_time){
        
        # Update time-relatedcharacteristics
        current_patient_update$ANLYEAR         <- current_patient$ANLYEAR  + current_pneumonia_time
        current_patient_update$ANLYEAR_SCALED  <- (current_patient_update$ANLYEAR - if(run_obs_input==1){0.962629}else{1.529411})/if(run_obs_input==1){0.9436266}else{1.376549} # hardcoded based on dataset
        current_patient_update$lag_ANLYEAR     <- current_pneumonia_time
        current_patient_update$age_time        <- current_patient$age_time + current_pneumonia_time
        current_patient_update$age_time_SCALED <- (current_patient_update$age_time - attr(scale(baseline_characteristics_run$age_time),"scaled:center"))/attr(scale(baseline_characteristics_run$age_time),"scaled:scale")
        
        
        # Force death at 100 years
        if(current_patient_update$age_time>100){current_patient_update$dead <- 1}
        
        # Update pneumonia status. 
        current_patient_update$pneu_yn      <- 1
        current_patient_update$pneu_hosp_yn <- rbinom(1,1,predicted_pneumonia_hosp(pneumonia_hosp_regression_coef,current_patient_update[pneumonia_hosp_predictors])$p.pneumonia.hosp)
        
        # Patient might die because of pneumonia after hospitalisation
        if(current_patient_update$pneu_hosp_yn==1){
          
          current_pneu_death_prob   <- 1607/19786 ### hardcoded based on dataset
          current_pneu_death_sample <- rbinom(1,1,current_pneu_death_prob)
          
          if(current_pneu_death_sample==1){
            current_patient_update$dead <- 1
          }else{
            current_patient_update$dead <- 0
          }
        }
        
      }
      
      
      # If death  happened then update age and finish the simulation for this patient
      if(min(current_remaining_life_exp,current_exacerbation_time)==current_remaining_life_exp){
        
        # Patient is dead
        current_patient_update$dead <- 1
        
        # Update time-relatedcharacteristics
        current_patient_update$ANLYEAR          <- current_patient$ANLYEAR  + current_remaining_life_exp
        current_patient_update$ANLYEAR_SCALED   <- (current_patient_update$ANLYEAR - if(run_obs_input==1){0.962629}else{1.529411})/if(run_obs_input==1){0.9436266}else{1.376549} # hardcoded based on dataset
        current_patient_update$lag_ANLYEAR      <- current_remaining_life_exp
        current_patient_update$age_time         <- current_patient$age_time + current_remaining_life_exp
        current_patient_update$age_time_SCALED  <- (current_patient_update$age_time - attr(scale(baseline_characteristics_run$age_time),"scaled:center"))/attr(scale(baseline_characteristics_run$age_time),"scaled:scale")
        
        # Update exacerbation status. 
        current_patient_update$modexac_yn <- 0
        current_patient_update$sevexac_yn <- 0
        
        # Update pneumonia status. 
        current_patient_update$pneu_yn      <- 0
        current_patient_update$pneu_hosp_yn <- 0
      }
      
      
      ############################################################################
      # STEP 3.2.3: Update continuous variables depending on the event occurred  #
      ############################################################################
      
      # Update FEV1, CWE, SGACT, symptoms and SGTOT. The order is important here:
      
      # Update FEV1
      current_patient_update$FEVA          <- max(0,predicted_fev1(run_obs_input,lung_function_regression_coef,current_patient_update[fev1_predictors],fev1_treatment_effect_input)$fev_1) 
      current_fevppa_calc                  <- fevppa_calc(current_patient_update$female,current_patient_update$HTSTD,current_patient_update$age_time,current_patient_update$FEVA)
      current_patient_update$fevppa        <- current_fevppa_calc$fevppa
      current_patient_update$FEV1pred      <- current_fevppa_calc$FEV1_pred
      current_patient_update$fevppa_SCALED <- (current_patient_update$fevppa - attr(scale(baseline_characteristics_run$fevppa),"scaled:center"))/attr(scale(baseline_characteristics_run$fevppa),"scaled:scale")
      
      # Force death if FEV1 < 0.2
      if(current_patient_update$FEVA < 0.2){current_patient_update$dead <- 1}
      
      
      # Update CWE (min observed is 42. here for now we truncate at 0)
      if(run_obs_input==0){
        current_patient_update$CWE_TOT        <- cwe_treatment_effect_input*max(0,predicted_cwe_tot(cwe_tot_regression_coef,current_patient_update[cwe_tot_predictors])$cwe_tot) #max(0,predicted_cwe_tot(cwe_tot_regression_coef$Value,current_patient_update[cwe_tot_predictors])$cwe_tot)
        current_patient_update$CWE_TOT_SCALED <- (current_patient_update$CWE_TOT - attr(scale(baseline_characteristics_run$CWE_TOT),"scaled:center"))/attr(scale(baseline_characteristics_run$CWE_TOT),"scaled:scale")
      }
      
      # Update SGACT
      current_patient_update$SGACT        <- min(100,max(0,sgact_treatment_effect_input+predicted_SGACT(SGACT_regression_coef,current_patient_update[SGACT_predictors])$SGACT))
      current_patient_update$SGACT_SCALED <- (current_patient_update$SGACT - attr(scale(baseline_characteristics_run$SGACT),"scaled:center"))/attr(scale(baseline_characteristics_run$SGACT),"scaled:scale")
      
      # Update breathlesness and coughsputum
      current_patient_update$breathlessyn  <- rbinom(1,1,breathless_treatment_effect_input*predicted_breathless(breathless_regression_coef,current_patient_update[breathless_predictors])$p.breathless)
      current_patient_update$coughsputumyn <- rbinom(1,1,coughsputum_treatment_effect_input*predicted_coughsputum(coughsputum_regression_coef,current_patient_update[coughsputum_predictors])$p.coughsputum)
      
      # Update SGTOT
      current_patient_update$SGTOT        <- min(100,max(0,sgtot_treatment_effect_input+predicted_SGTOT(SGTOT_regression_coef,current_patient_update[SGTOT_predictors])$SGTOT))
      current_patient_update$SGTOT_SCALED <- (current_patient_update$SGTOT - attr(scale(baseline_characteristics_run$SGTOT),"scaled:center"))/attr(scale(baseline_characteristics_run$SGTOT),"scaled:scale")
      
      ### When all characteristics have been updated, we add these to the patient history
      simulation_patients_history <- rbind(simulation_patients_history,current_patient_update[history_characteristics])
      
      ### And update current patient and go up to while loop
      current_patient <- current_patient_update
      
      
      ######################################################################
      # STEP 3.2.4: Adjust Mortality according to updated characteristics  #
      ######################################################################
      
      ### Fixed to Weibull
      current_mortality_parameters <- predicted_mortality_weibull(mortality_weibull_regression_coef,current_patient[mortality_predictors])
      current_mortality            <- current_mortality_parameters$scale_mortality_weibull*gamma(1+1/current_mortality_parameters$shape_mortality_weibull)/365
      
      ### Adjust remaining life expectancy according to improvement or worsened in condition wrt baseline or previous period
      current_remaining_life_exp <- max(0,(current_remaining_life_exp - current_patient$lag_ANLYEAR)*(current_mortality/lag_current_mortality))
      
      ### Update mortality factor and event index
      lag_current_mortality    <- current_mortality 
      current_event            <- current_event + 1
      
    } #end while loop and move to another patient
    
    patient_index <- patient_index + 1
    
  } #end for loop in number of patients
  
  # Output from PART I
  # View(simulation_patients_history)
  
  #################################################################################
  ########## MAIN PART II: Update intermediate outcomes after every year ##########
  #################################################################################
  
  patient_characteristics_saved <- if(run_obs_input==0){c("SIMID","PTID","ANLYEAR","age_time","FEVA","fevppa","sevexac_yn","modexac_yn","CWE_TOT",
                                                          "SGACT","SGTOT","coughsputumyn","breathlessyn","pneu_yn","pneu_hosp_yn","dead")}
                                                   else{c("SIMID","PTID","ANLYEAR","age_time","FEVA","fevppa","sevexac_yn","modexac_yn",
                                                          "SGACT","SGTOT","coughsputumyn","breathlessyn","pneu_yn","pneu_hosp_yn","dead")}
  
  
  #patient_event_history_update <- simulation_baseline_patients[FALSE,c(patient_characteristics_saved)]
  patient_event_history_update <- simulation_patients_history[FALSE,c(patient_characteristics_saved)]
  
  
  for(i in 1:max(simulation_patients_history$SIMID)){
    
    if(run_PSA_input == 0){print(i+patient_size_input)}
    
    current_patient_event_history <- simulation_patients_history[which(simulation_patients_history$SIMID == i),]
    where    <- tail(1:nrow(current_patient_event_history),-1) 
    how_many <- floor(diff(current_patient_event_history$ANLYEAR))### Not correct either: this works when the time to event is larger than 1 year! So floor has to be larger than 0
    
    current_patient_event_history_update <- data.frame(matrix(,ncol=ncol(current_patient_event_history),nrow=max(where)+sum(how_many)))
    colnames(current_patient_event_history_update) <- c(history_characteristics)
    
    current_patient_event_history_update[c(1,cumsum(how_many)+where),] <- current_patient_event_history[1:nrow(current_patient_event_history),]
    current_patient_event_history_update$SIMID <- current_patient_event_history$SIMID[1]
    current_patient_event_history_update$PTID  <- current_patient_event_history$PTID[1]
    
    for(j in 2:(nrow(current_patient_event_history_update)-1)){
      
      if(is.na(current_patient_event_history_update[j,]$ANLYEAR)==TRUE){
        
        current_patient_event_history_update[j,]$female        <- current_patient_event_history_update[j-1,]$female
        current_patient_event_history_update[j,]$AGE           <- current_patient_event_history_update[j-1,]$AGE
        current_patient_event_history_update[j,]$AGE_SCALED    <- current_patient_event_history_update[j-1,]$AGE_SCALED
        current_patient_event_history_update[j,]$BMIclass2     <- current_patient_event_history_update[j-1,]$BMIclass2  
        current_patient_event_history_update[j,]$BMIclass3     <- current_patient_event_history_update[j-1,]$BMIclass3
        current_patient_event_history_update[j,]$SMOKER        <- current_patient_event_history_update[j-1,]$SMOKER
        current_patient_event_history_update[j,]$SMPKY         <- current_patient_event_history_update[j-1,]$SMPKY
        current_patient_event_history_update[j,]$SMPKY_SCALED  <- current_patient_event_history_update[j-1,]$SMPKY_SCALED
        current_patient_event_history_update[j,]$other_CVD     <- current_patient_event_history_update[j-1,]$other_CVD  
        current_patient_event_history_update[j,]$Reversibility        <- current_patient_event_history_update[j-1,]$Reversibility
        current_patient_event_history_update[j,]$Reversibility_SCALED <- current_patient_event_history_update[j-1,]$Reversibility_SCALED 
        current_patient_event_history_update[j,]$MH_DI                <- current_patient_event_history_update[j-1,]$MH_DI
        current_patient_event_history_update[j,]$MH_DE                <- current_patient_event_history_update[j-1,]$MH_DE
        current_patient_event_history_update[j,]$MH_CF                <- current_patient_event_history_update[j-1,]$MH_CF
        current_patient_event_history_update[j,]$MH_RR                <- current_patient_event_history_update[j-1,]$MH_RR
        current_patient_event_history_update[j,]$EMPHDIA              <- current_patient_event_history_update[j-1,]$EMPHDIA
        if(run_obs_input==0){
          current_patient_event_history_update[j,]$EOS_yn               <- current_patient_event_history_update[j-1,]$EOS_yn
        }
        current_patient_event_history_update[j,]$ICS                  <- current_patient_event_history_update[j-1,]$ICS
        current_patient_event_history_update[j,]$FEVA_BL              <- current_patient_event_history_update[j-1,]$FEVA_BL
        current_patient_event_history_update[j,]$HTSTD                <- current_patient_event_history_update[j-1,]$HTSTD
        
        
        current_patient_event_history_update[j,]$ANLYEAR         <- current_patient_event_history_update[j-1,]$ANLYEAR + 1 
        current_patient_event_history_update[j,]$ANLYEAR_SCALED  <- (current_patient_event_history_update[j,]$ANLYEAR - if(run_obs_input==1){0.962629}else{1.529411})/if(run_obs_input==1){0.9436266}else{1.376549}
        
        current_patient_event_history_update[j,]$age_time     <- current_patient_event_history_update[j-1,]$age_time + 1 
        current_patient_event_history_update[j,]$sevexac_yn   <- 0
        current_patient_event_history_update[j,]$modexac_yn   <- 0
        current_patient_event_history_update[j,]$pneu_yn      <- 0
        current_patient_event_history_update[j,]$pneu_hosp_yn <- 0
        current_patient_event_history_update[j,]$dead         <- 0
        
        
        current_patient_event_history_update[j,]$FEVA         <- max(0,predicted_fev1(run_obs_input,lung_function_regression_coef,current_patient_event_history_update[j,][fev1_predictors],fev1_treatment_effect_input)$fev_1) 
        current_patient_event_history_update_fevppa_calc      <- fevppa_calc(current_patient_event_history_update[j,]$female,current_patient_event_history_update[j,]$HTSTD,current_patient_event_history_update[j,]$age_time,current_patient_event_history_update[j,]$FEVA)
        current_patient_event_history_update[j,]$fevppa       <- current_patient_event_history_update_fevppa_calc$fevppa
        #current_patient_event_history_update[j,]$FEV1pred     <- current_patient_event_history_update_fevppa_calc$FEV1_pred
        current_patient_event_history_update[j,]$fevppa_SCALED <- (current_patient_event_history_update[j,]$fevppa - attr(scale(baseline_characteristics_run$fevppa),"scaled:center"))/attr(scale(baseline_characteristics_run$fevppa),"scaled:scale")
        
        current_patient_event_history_update[j,]$lag_SGACT    <- current_patient_event_history_update[j-1,]$SGACT 
        
        if(run_obs_input==0){
          current_patient_event_history_update[j,]$lag_CWE_TOT    <- current_patient_event_history_update[j-1,]$CWE_TOT
          current_patient_event_history_update[j,]$CWE_TOT        <- max(0,predicted_cwe_tot(cwe_tot_regression_coef,current_patient_event_history_update[j,][cwe_tot_predictors])$cwe_tot) #max(0,predicted_cwe_tot(cwe_tot_regression_coef$Value,current_patient_event_history_update[j,][cwe_tot_predictors])$cwe_tot)    
          current_patient_event_history_update[j,]$CWE_TOT_SCALED <- (current_patient_event_history_update[j,]$CWE_TOT - attr(scale(baseline_characteristics_run$CWE_TOT),"scaled:center"))/attr(scale(baseline_characteristics_run$CWE_TOT),"scaled:scale")
        }
        
        
        current_patient_event_history_update[j,]$lag_breathlessyn  <- current_patient_event_history_update[j-1,]$breathlessyn
        current_patient_event_history_update[j,]$lag_coughsputumyn <- current_patient_event_history_update[j-1,]$coughsputumyn
        current_patient_event_history_update[j,]$lag_SGTOT         <- current_patient_event_history_update[j-1,]$SGTOT
        
        
        current_patient_event_history_update[j,]$SGACT             <- min(100,max(0,predicted_SGACT(SGACT_regression_coef,current_patient_event_history_update[j,][SGACT_predictors])$SGACT))
        current_patient_event_history_update[j,]$SGACT_SCALED      <- (current_patient_event_history_update[j,]$SGACT - attr(scale(baseline_characteristics_run$SGACT),"scaled:center"))/attr(scale(baseline_characteristics_run$SGACT),"scaled:scale")
        
        
        current_patient_event_history_update[j,]$breathlessyn  <- rbinom(1,1,rbinom(1,1,predicted_breathless(breathless_regression_coef,current_patient_event_history_update[j,][breathless_predictors])$p.breathless))
        current_patient_event_history_update[j,]$coughsputumyn <- rbinom(1,1,rbinom(1,1,predicted_coughsputum(coughsputum_regression_coef,current_patient_event_history_update[j,][coughsputum_predictors])$p.coughsputum))
        
        
        current_patient_event_history_update[j,]$SGTOT         <- min(100,max(0,predicted_SGTOT(SGTOT_regression_coef,current_patient_event_history_update[j,][SGTOT_predictors])$SGTOT))
        #current_patient_update$SGTOT_SCALED <- (current_patient_update$SGTOT - attr(scale(baseline_characteristics_run$SGTOT),"scaled:center"))/attr(scale(baseline_characteristics_run$SGTOT),"scaled:scale")
        
      }
      
    } #end for per patient
    
    patient_event_history_update <- rbind(patient_event_history_update,current_patient_event_history_update[,c(patient_characteristics_saved)])
    
  }
  
  # Output from PART II
  # View(patient_event_history_update)
  
  #################################################################
  ########## MAIN PART III: Calculate aggregated results ##########
  #################################################################
  
  ### We first need to make additional columns for the "diff" variables
  patient_event_history_update$diff_ANLYEAR <- "NA"
  patient_event_history_update$diff_FEVA    <- "NA"
  patient_event_history_update$diff_SGTOT   <- "NA"
  patient_event_history_update$diff_SGACT   <- "NA"
  if(run_obs_input==0){patient_event_history_update$diff_CWE_TOT <- "NA"}
  
  ### Calculate the "diff" variables to calculate the change in outcomes per year
  diff_ANLYEAR <- ddply(patient_event_history_update, "SIMID", summarize, diff_ANLYEAR = c(0,diff(ANLYEAR)))
  diff_FEVA    <- ddply(patient_event_history_update, "SIMID", summarize, diff_FEVA    = c(0,diff(FEVA)))
  diff_SGTOT   <- ddply(patient_event_history_update, "SIMID", summarize, diff_SGTOT   = c(0,diff(SGTOT)))
  diff_SGACT   <- ddply(patient_event_history_update, "SIMID", summarize, diff_SGACT   = c(0,diff(SGACT)))
  if(run_obs_input==0){diff_CWE_TOT <- ddply(patient_event_history_update, "SIMID", summarize, diff_CWE_TOT = c(0,diff(CWE_TOT)))}
  
  patient_event_history_update$diff_ANLYEAR <- diff_ANLYEAR$diff_ANLYEAR
  patient_event_history_update$diff_FEVA    <- diff_FEVA$diff_FEVA
  patient_event_history_update$diff_SGTOT   <- diff_SGTOT$diff_SGTOT
  patient_event_history_update$diff_SGACT   <- diff_SGACT$diff_SGACT
  if(run_obs_input==0){patient_event_history_update$diff_CWE_TOT <- diff_CWE_TOT$diff_CWE_TOT}
  
  ### Calculate model outcomes ###
  
  ### Adjust the number of exacerbations (no exacerbation when pneumonia happened)
  patient_event_history_update[which(patient_event_history_update$pneu_yn==1),"modexac_yn"] <- 0
  patient_event_history_update[which(patient_event_history_update$pneu_yn==1),"sevexac_yn"] <- 0
  
  
  ### Lung function: Mean FEV1 decline per year
  mean_annual_fev1_decline <- round(mean(patient_event_history_update$diff_FEVA)/mean(patient_event_history_update$diff_ANLYEAR),4)
  
  ### Life expectancy
  mean_life_expectancy <- round(mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"]),4)
  
  # Moderate
  modexac_year <-function(low_limit,upper_limit){
    sum(patient_event_history_update[which(patient_event_history_update$ANLYEAR>low_limit & patient_event_history_update$ANLYEAR<=upper_limit),"modexac_yn"])
  }
  mean_mod_exa_rate <- round(modexac_year(0,1000)/(patient_size_input*mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"])),4)
  
  # Severe 
  sevexac_year <-function(low_limit,upper_limit){
    sum(patient_event_history_update[which(patient_event_history_update$ANLYEAR>low_limit & patient_event_history_update$ANLYEAR<=upper_limit),"sevexac_yn"])
  }
  mean_sev_exa_rate <- round(sevexac_year(0,1000)/(patient_size_input*mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"])),4)
  
  
  ### Exercise capacity
  mean_CWE_TOT_change <- if(run_obs_input==0){round(mean(patient_event_history_update$diff_CWE_TOT)/mean(patient_event_history_update$diff_ANLYEAR),4)}else{0}
  
  ### SGRQ activity score 
  mean_SGACT_change <- round(mean(patient_event_history_update$diff_SGACT)/mean(patient_event_history_update$diff_ANLYEAR),4)
  
  ### Cough/sputum
  # Total number of cough/sputum per patient during lifetime
  cum_coughsputum <- aggregate(patient_event_history_update$coughsputumyn, list(Patient = patient_event_history_update$SIMID), sum)
  # Rate per year
  mean_cough_rate <- round(mean(cum_coughsputum$x)/mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"]),4)
  
  # July 2018 added
  # create new time variable
  diff_ANLYEAR_sym <- ddply(patient_event_history_update, "SIMID", summarize, diff_ANLYEAR_sym = c(diff(ANLYEAR),0))
  
  # create data frame just for symptoms
  simulation_clinical_history_sym <- cbind(patient_event_history_update[,c("ANLYEAR", "coughsputumyn", "breathlessyn","dead")], diff_ANLYEAR_sym)
  simulation_dead                 <- simulation_clinical_history_sym[which(simulation_clinical_history_sym$dead == 1),c("SIMID", "ANLYEAR")]
  simulation_cough                <- simulation_clinical_history_sym[,c("SIMID", "diff_ANLYEAR_sym")]
  simulation_cough_time           <- simulation_clinical_history_sym$coughsputumyn * simulation_clinical_history_sym$diff_ANLYEAR_sym
  simulation_cough                <- cbind(simulation_cough, simulation_cough_time)
  mean_time_cough                 <- mean(aggregate(simulation_cough$simulation_cough_time, list(Patient = simulation_cough$SIMID), sum)$x/simulation_dead$ANLYEAR)
  #ci_time_cough         <- quantile(aggregate(simulation_cough$simulation_cough_time, list(Patient = simulation_cough$SIMID), sum)$x/simulation_dead$ANLYEAR, c(0.025,0.975))
  # End - July 2018 added
  
  ### Shortness of breath
  # Total number of cough/sputum per patient during lifetime
  cum_breathless <- aggregate(patient_event_history_update$breathlessyn, list(Patient = patient_event_history_update$SIMID), sum)
  # Rate per year
  mean_breathless_rate <- round(mean(cum_breathless$x)/mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"]),4)
  
  
  # July 2018 added
  simulation_breath      <- simulation_clinical_history_sym[,c("SIMID", "diff_ANLYEAR_sym")]
  simulation_breath_time <- simulation_clinical_history_sym$breathlessyn * simulation_clinical_history_sym$diff_ANLYEAR_sym
  simulation_breath      <- cbind(simulation_breath, simulation_breath_time)
  mean_time_breath       <- mean(aggregate(simulation_breath$simulation_breath_time, list(Patient = simulation_breath$SIMID), sum)$x/simulation_dead$ANLYEAR)
  # End - July 2018 added
  
  
  
  ### Adverse events (pneumonia)
  # Total number of pneumonias per patient during lifetime
  cum_pneu <- aggregate(patient_event_history_update$pneu_yn, list(Patient = patient_event_history_update$SIMID), sum)
  # Rate per year
  mean_pneu_rate <- round(mean(cum_pneu$x)/mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"]),4)
  
  # Total number of pneumonias leading to hospitalisation per patient during lifetime
  cum_pneu_hosp <- aggregate(patient_event_history_update$pneu_hosp_yn, list(Patient = patient_event_history_update$SIMID), sum)
  # Rate per year
  mean_pneu_hosp_rate <- round(mean(cum_pneu_hosp$x)/mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"]),4)
  
  
  ### SGRQ total score 
  mean_SGTOT_change <- round(mean(patient_event_history_update$diff_SGTOT)/mean(patient_event_history_update$diff_ANLYEAR),4)
  
  
  ### QALYs
  # Read baseline characteristics: we need gender to calculate utilities
  ## In the average case it's not needed for the utilities but for the maintenace costs
  baseline_characteristics <- read.csv("Model - datasets/baseline_characteristics_predicted_data.csv",sep=",")
  
  # Add column with gender to the simulated results file
  ##AVERAGE
  patient_event_history_update <- merge(patient_event_history_update,baseline_characteristics[c("PTID","female")],by.x = "PTID",by.y = "PTID")
  #patient_event_history_update$female <- (simulation_patients_history$female)[1]
  
  # Order the dataset by "SIMID". This is very important because after merging the order is lost.
  patient_event_history_update <- patient_event_history_update[order(patient_event_history_update$SIMID,patient_event_history_update$ANLYEAR),]
  

  ### These are the simulated utilities: Changed with respect to the previous version: now it's vectorized and therefore faster.
  utilities <- round(0.9617-0.0013*(patient_event_history_update$SGTOT)-0.0001*((patient_event_history_update$SGTOT)^2) + ifelse(patient_event_history_update$female==0,0.0231,0),4)
  
  
  # Add them to the results dataset 
  patient_event_history_update$utilities <- utilities
  
  # To calculate QALYS we multiply utilities by diff_ANLYEAR, then discount them
  QALYs               <- patient_event_history_update$diff_ANLYEAR*patient_event_history_update$utilities
  
  discount_rate_QALYs <- ifelse(patient_event_history_update$ANLYEAR<=30,0.035,0.035) # UK = 0.035 always #Dynagito 0.04 0.02 #hardcoded
  QALYs_discounted    <- QALYs/(1+discount_rate_QALYs)^patient_event_history_update$ANLYEAR
  
  # Add QALYs to the results table
  patient_event_history_update$QALYs <- QALYs
  patient_event_history_update$QALYs_discounted <- QALYs_discounted
  
  # Finally, aggregate QALYS per patient and compute the average
  QALYs_patient <- aggregate(patient_event_history_update$QALYs,list(SIMID=patient_event_history_update$SIMID),sum)
  mean_qalys <- round(mean(QALYs_patient$x),4)
  
  # discounted QALYs
  QALYs_patient_discounted <- aggregate(patient_event_history_update$QALYs_discounted,list(SIMID=patient_event_history_update$SIMID),sum)
  mean_qalys_disc <- round(mean(QALYs_patient_discounted$x),4)
  
  
  
  ### Costs
  
  # Treatment  costs: in France there's a distinction between HC and Soc.
  # Scenario: mind the price when considering treatment effect scenarios.
  # Be careful when implementing this in the interface!!!
  treatment_price_year_hc       <- 1.12*1*365.25 #if(exac_treatment_effect_tte_input ==1 && fev1_treatment_effect_input==1){1.12*365.25}else{1.12*1.20*365.25} ### Hardcoded here
  treatment_price_year_societal <- treatment_price_year_hc
  
  
  patient_event_history_update$treatment_costs_hc            <- patient_event_history_update$diff_ANLYEAR*treatment_price_year_hc
  discount_rate_costs <- ifelse(patient_event_history_update$ANLYEAR<=30,0.035,0.035) # UK = 0.035 always #Dynagito 0.04 0.02 #hardcoded
  patient_event_history_update$treatment_costs_hc_discounted <- patient_event_history_update$treatment_costs_hc/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  
  patient_event_history_update$treatment_costs_societal            <- patient_event_history_update$diff_ANLYEAR*treatment_price_year_societal
  patient_event_history_update$treatment_costs_societal_discounted <- patient_event_history_update$treatment_costs_societal/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  
  # Exacerbation costs 
  
  # Read costs from Excel file depending on the selected options
  
  # Retirement age
  retirement_age <- 65 # UK = 65 #Dynagito 62 #hardcoded
  
  ### Societal persp. & health care use
  exacerbation_costs_hc_use_row_names <- c("Primary care visits", "Secondary care visits", "Hospital days",
                                           "Ambulance rides", "ER visits", "Course antibiotics","Course oral steroids",
                                           "Work days lost", "Distance primary care clinic", "Distance specialist clinic")
  
  exacerbation_costs_hc_use <- read.csv("Model - datasets/Costs/exacerbation costs hc use.csv",sep=";",row.names = exacerbation_costs_hc_use_row_names)
  
  
  # Societal persp. & average
  exacerbation_costs_average_row_names <- c("Societal", "Societal (retired)", "Health care", "Health care (retired)")
  
  exacerbation_costs_average <- read.csv("Model - datasets/Costs/exacerbation costs average.csv",sep=";",row.names = exacerbation_costs_average_row_names)
  
  ### Dynagito scenario
  #exacerbation_costs_average <- read.csv("Model - datasets/Costs/exacerbation costs average - FR.csv",sep=",",row.names = exacerbation_costs_average_row_names)
  
  mod_exa_costs_societal_average <- exacerbation_costs_average[1,1]   
  sev_exa_costs_societal_average <- exacerbation_costs_average[1,2]   
  
  mod_exa_costs_societal_average_retired <- exacerbation_costs_average[2,1]   
  sev_exa_costs_societal_average_retired <- exacerbation_costs_average[2,2]   
  
  ### Health care persp. & average
  mod_exa_costs_hc_average <- exacerbation_costs_average[3,1]   
  sev_exa_costs_hc_average <- exacerbation_costs_average[3,2]   
  
  mod_exa_costs_hc_average_retired <- exacerbation_costs_average[4,1]   
  sev_exa_costs_hc_average_retired <- exacerbation_costs_average[4,2]   
  
  
  ### Write functions to calculate cost per moderate and severe exacerbation
  ### For the moment I include some items but this can change
  mod_exa_cost <- function(age_input){
    
    if(cost_type=="Health care use"){
      mod_exa_primary_care_cost   <- exacerbation_costs_hc_use[1,1]*exacerbation_costs_hc_use[1,3] 
      mod_exa_secondary_care_cost <- exacerbation_costs_hc_use[2,1]*exacerbation_costs_hc_use[2,3] 
      mod_exa_hospital_days_cost  <- 0 
      mod_exa_ambulance_ride_cost <- exacerbation_costs_hc_use[4,1]*exacerbation_costs_hc_use[4,3] 
      mod_exa_ER_visit_cost       <- exacerbation_costs_hc_use[5,1]*exacerbation_costs_hc_use[5,3] 
      mod_exa_antibiotics_cost    <- exacerbation_costs_hc_use[6,1]*exacerbation_costs_hc_use[6,3]
      mod_exa_steroids_cost       <- exacerbation_costs_hc_use[7,1]*exacerbation_costs_hc_use[7,3]
      mod_exa_work_days_lost_cost <- ifelse(perspective=="Societal" & age_input<=retirement_age, exacerbation_costs_hc_use[8,1]*exacerbation_costs_hc_use[8,3], 0) 
      mod_exa_distance_pcc_cost   <- ifelse(perspective=="Societal", 2*exacerbation_costs_hc_use[1,1]*exacerbation_costs_hc_use[9,1]*exacerbation_costs_hc_use[9,3], 0)
      mod_exa_distance_sc_cost    <- ifelse(perspective=="Societal", 2*exacerbation_costs_hc_use[2,1]*exacerbation_costs_hc_use[10,1]*exacerbation_costs_hc_use[10,3], 0)
      
    }
    
    if(cost_type=="Average"){
      mod_exa_average_cost   <- ifelse(perspective=="Societal", 
                                       ifelse(age_input<=retirement_age,mod_exa_costs_societal_average,mod_exa_costs_societal_average_retired), 
                                       ifelse(age_input<=retirement_age,mod_exa_costs_hc_average,mod_exa_costs_hc_average_retired))
    }
    
    #mod_exa_cost_total <- 
    ifelse(cost_type=="Health care use",
           sum(mod_exa_primary_care_cost,
               mod_exa_secondary_care_cost,
               mod_exa_hospital_days_cost,
               mod_exa_ambulance_ride_cost,
               mod_exa_ER_visit_cost,
               mod_exa_antibiotics_cost,
               mod_exa_steroids_cost,
               mod_exa_work_days_lost_cost,
               mod_exa_distance_pcc_cost,
               mod_exa_distance_sc_cost),
               mod_exa_average_cost)
    
  }
  
  
  
  # Severe exacerbations
  sev_exa_cost <- function(age_input){
    
    if(cost_type=="Health care use"){
      sev_exa_primary_care_cost   <- exacerbation_costs_hc_use[1,2]*exacerbation_costs_hc_use[1,3] 
      sev_exa_secondary_care_cost <- exacerbation_costs_hc_use[2,2]*exacerbation_costs_hc_use[2,3] 
      sev_exa_hospital_days_cost  <- exacerbation_costs_hc_use[3,2]*exacerbation_costs_hc_use[3,3]  
      sev_exa_ambulance_ride_cost <- exacerbation_costs_hc_use[4,2]*exacerbation_costs_hc_use[4,3] 
      sev_exa_ER_visit_cost       <- exacerbation_costs_hc_use[5,2]*exacerbation_costs_hc_use[5,3] 
      sev_exa_antibiotics_cost    <- exacerbation_costs_hc_use[6,2]*exacerbation_costs_hc_use[6,3]
      sev_exa_steroids_cost       <- exacerbation_costs_hc_use[7,2]*exacerbation_costs_hc_use[7,3]
      sev_exa_work_days_lost_cost <- ifelse(perspective=="Societal" & age_input<=retirement_age, exacerbation_costs_hc_use[8,2]*exacerbation_costs_hc_use[8,3], 0) 
      sev_exa_distance_pcc_cost   <- ifelse(perspective=="Societal", 2*exacerbation_costs_hc_use[1,2]*exacerbation_costs_hc_use[9,2]*exacerbation_costs_hc_use[9,3], 0)
      sev_exa_distance_sc_cost    <- ifelse(perspective=="Societal", 2*exacerbation_costs_hc_use[2,2]*exacerbation_costs_hc_use[10,2]*exacerbation_costs_hc_use[10,3], 0)
      
    }
    
    if(cost_type=="Average"){
      sev_exa_average_cost   <- ifelse(perspective=="Societal", 
                                       ifelse(age_input<=retirement_age,sev_exa_costs_societal_average,sev_exa_costs_societal_average_retired), 
                                       ifelse(age_input<=retirement_age,sev_exa_costs_hc_average,sev_exa_costs_hc_average_retired))
    }
    
    #sev_exa_cost_total <- 
    ifelse(cost_type=="Health care use",
           sum(sev_exa_primary_care_cost,
               sev_exa_secondary_care_cost,
               sev_exa_hospital_days_cost,
               sev_exa_ambulance_ride_cost,
               sev_exa_ER_visit_cost,
               sev_exa_antibiotics_cost,
               sev_exa_steroids_cost,
               sev_exa_work_days_lost_cost,
               sev_exa_distance_pcc_cost,
               sev_exa_distance_sc_cost),
           sev_exa_average_cost)
    
  }
  
  
  
  ### Now calculate the exacerbation costs associated to exacerbations in the simulation
  ### This has to be age-dependent when the perspective is societal (productivity lossses)
  
  
  ### And apply it to the model results
  sev_exacerbation_costs_calc_sim <- function(index){
    if(patient_event_history_update[index,"ANLYEAR"]>0){
      sev_exa_cost(patient_event_history_update[index,"age_time"])*patient_event_history_update[index,"sevexac_yn"]
    }else{0}
    
  }
  
  
  mod_exacerbation_costs_calc_sim <- function(index){
    if(patient_event_history_update[index,"ANLYEAR"]>0){
      mod_exa_cost(patient_event_history_update[index,"age_time"])*patient_event_history_update[index,"modexac_yn"]
    }else{0}
  }
  
  
  ### Choose first perspective and type of cost data: these should be input parameters!!!
  cost_type   <- "Health care use" #"Health care use" ### Choose between "Health care use" or "Average" # Dynagito uses Average
  perspective <- "Societal" ### Choose between "Societal" or "Health care"
  
  
  ### These are the simulated exacerbation costs
  sev_exa_costs_societal            <- sapply(1:nrow(patient_event_history_update),sev_exacerbation_costs_calc_sim)
  sev_exa_costs_societal_discounted <- sev_exa_costs_societal/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  mod_exa_costs_societal            <- sapply(1:nrow(patient_event_history_update),mod_exacerbation_costs_calc_sim)
  mod_exa_costs_societal_discounted <- mod_exa_costs_societal/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  
  patient_event_history_update$sev_exa_costs_societal            <- sev_exa_costs_societal
  patient_event_history_update$sev_exa_costs_societal_discounted <- sev_exa_costs_societal_discounted
  patient_event_history_update$mod_exa_costs_societal            <- mod_exa_costs_societal
  patient_event_history_update$mod_exa_costs_societal_discounted <- mod_exa_costs_societal_discounted
  
  patient_event_history_update$exacerbation_costs_societal            <- patient_event_history_update$sev_exa_costs_societal + patient_event_history_update$mod_exa_costs_societal
  patient_event_history_update$exacerbation_costs_societal_discounted <- patient_event_history_update$sev_exa_costs_societal_discounted + patient_event_history_update$mod_exa_costs_societal_discounted
  
  # Exacerbation costs aggregated per patient and average 
  mod_exa_cost_societal_patient <- aggregate(patient_event_history_update$mod_exa_costs_societal,list(SIMID=patient_event_history_update$SIMID),sum)
  sev_exa_cost_societal_patient <- aggregate(patient_event_history_update$sev_exa_costs_societal,list(SIMID=patient_event_history_update$SIMID),sum)
  
  
  ### Discounted
  mod_exa_cost_societal_patient_discounted <- aggregate(patient_event_history_update$mod_exa_costs_societal_discounted,list(SIMID=patient_event_history_update$SIMID),sum)
  sev_exa_cost_societal_patient_discounted <- aggregate(patient_event_history_update$sev_exa_costs_societal_discounted,list(SIMID=patient_event_history_update$SIMID),sum)
  
  
  ### Choose first perspective and type of cost data
  perspective <- "Health care" ### Choose between "Societal" or "Health care"
  
  
  ### These are the simulated exacerbation costs
  sev_exa_costs_hc            <- sapply(1:nrow(patient_event_history_update),sev_exacerbation_costs_calc_sim)
  sev_exa_costs_hc_discounted <- sev_exa_costs_hc/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  mod_exa_costs_hc            <- sapply(1:nrow(patient_event_history_update),mod_exacerbation_costs_calc_sim)
  mod_exa_costs_hc_discounted <- mod_exa_costs_hc/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  
  patient_event_history_update$sev_exa_costs_hc            <- sev_exa_costs_hc
  patient_event_history_update$sev_exa_costs_hc_discounted <- sev_exa_costs_hc_discounted
  patient_event_history_update$mod_exa_costs_hc            <- mod_exa_costs_hc
  patient_event_history_update$mod_exa_costs_hc_discounted <- mod_exa_costs_hc_discounted
  
  patient_event_history_update$exacerbation_costs_hc            <- patient_event_history_update$sev_exa_costs_hc + patient_event_history_update$mod_exa_costs_hc
  patient_event_history_update$exacerbation_costs_hc_discounted <- patient_event_history_update$sev_exa_costs_hc_discounted + patient_event_history_update$mod_exa_costs_hc_discounted
  
  ### Exacerbation costs aggregated per patient and average 
  mod_exa_cost_hc_patient <- aggregate(patient_event_history_update$mod_exa_costs_hc,list(SIMID=patient_event_history_update$SIMID),sum)
  sev_exa_cost_hc_patient <- aggregate(patient_event_history_update$sev_exa_costs_hc,list(SIMID=patient_event_history_update$SIMID),sum)
  
  # Exacerbation costs discounted
  mod_exa_cost_hc_patient_discounted <- aggregate(patient_event_history_update$mod_exa_costs_hc_discounted,list(SIMID=patient_event_history_update$SIMID),sum)
  sev_exa_cost_hc_patient_discounted <- aggregate(patient_event_history_update$sev_exa_costs_hc_discounted,list(SIMID=patient_event_history_update$SIMID),sum)
  
  
  ## Maintenance costs 
  
  
  ### GP and specialist visits
  
  ### Need to define predictive functions like in the simulation core of the model
  ### Then apply these functions per row of the simulation results file as above (after merging baseline patient characteristics)
  
  # GP visits
  gpvisits_regression_coef <- read.csv("Model - regression coefficients/Costs/gpvisits_regression_coef_observed_data.csv",sep=";")
  
  # Predict GP visits
  gpvisits_stable_predictors <- c("female","BMIclass2","BMIclass3","SMOKER","SMPKY","other_CVD"
                                  ,"MH_DI","MH_DE","MH_CF","MH_RR","EMPHDIA","ICS","EOS_yn")
  
  gpvisits_predictors <- c(gpvisits_stable_predictors,"age_time","fevppa","SGACT","coughsputumyn","breathlessyn","SGTOT","modexac_yn","sevexac_yn")
  
  
  ### Merge the two datsets: I'm creating another file here because I don't want to save all the patient
  ### characteristics in the simulated results file
  #AVERAGE
  patient_event_history_update_maintenance_costs <- merge(patient_event_history_update[,1:19],baseline_characteristics[c("PTID",gpvisits_stable_predictors)],by.x = "PTID",by.y = "PTID")
  
  #patient_event_history_update_maintenance_costs <- merge(patient_event_history_update[,1:19],baseline_characteristics_run_average[,c(gpvisits_stable_predictors)])
  
  ### Order the dataset by "SIMID". This is very important because after merging the order is lost.
  patient_event_history_update_maintenance_costs <- patient_event_history_update_maintenance_costs[order(patient_event_history_update_maintenance_costs$SIMID),]
  
  predicted_gpvisits <- function(regression_coefficents_gpvisits_input, patient_characteristics_gpvisits_input){
    patient_characteristics_gpvisits_input <- as.numeric(patient_characteristics_gpvisits_input)
    gpvisits <- sum(regression_coefficents_gpvisits_input*c(1,patient_characteristics_gpvisits_input))
    return(list(gpvisits=gpvisits))
    
  }
  
  # Specialist visits
  specvisits_regression_coef <- read.csv("Model - regression coefficients/Costs/specvisits_regression_coef_observed_data.csv",sep=",")
  
  # Predict spec visits
  specvisits_predictors <- gpvisits_predictors
  
  predicted_specvisits <- function(regression_coefficents_specvisits_input, patient_characteristics_specvisits_input){
    patient_characteristics_specvisits_input <- as.numeric(patient_characteristics_specvisits_input)
    specvisits <- sum(regression_coefficents_specvisits_input*c(1,patient_characteristics_specvisits_input))
    return(list(specvisits=specvisits))
    
  }
  
  
  
  ### And finally predict for the simulation results
  predicted_gpvisits_calc <-function(index){
    exp(predicted_gpvisits(gpvisits_regression_coef$Estimate,patient_event_history_update_maintenance_costs[index,gpvisits_predictors])$gpvisits)
  }
  
  num_gpvisits <- sapply(1:nrow(patient_event_history_update),predicted_gpvisits_calc) 
  
  
  predicted_specvisits_calc <-function(index){
    exp(predicted_specvisits(specvisits_regression_coef$Estimate,patient_event_history_update_maintenance_costs[index,specvisits_predictors])$specvisits)
  }
  
  num_specvisits <- sapply(1:nrow(patient_event_history_update),predicted_specvisits_calc)
  
  ### Add these to the simulated results file
  patient_event_history_update$num_gpvisits   <- num_gpvisits
  patient_event_history_update$num_specvisits <- num_specvisits
  
  
  ### Read first the maintenance costs data
  maintenance_costs_row_names <- c("Number primary care visits per year",
                                   "Number specialist visits per year",
                                   "Spirometries per year",
                                   "Influenza vaccination per year",
                                   "ICS costs/year",
                                   "Distance primary clinic",
                                   "Distance specialist clinic")
  
  # Base case
  maintenance_costs <- read.csv("Model - datasets/Costs/maintenance costs.csv",sep=";",row.names = maintenance_costs_row_names)
  
  #Dynagito scenario
  #maintenance_costs <- read.csv("Model - datasets/Costs/maintenance costs - FR.csv",sep=",",row.names = maintenance_costs_row_names)
  
  #maintenance_costs
  
  ### Adjust then the number of gp and spec visits and multiply by the time interval where occurs
  num_gpvisits_observed <- 0.422
  num_gpvisits_country  <- maintenance_costs[1,1]
  patient_event_history_update$adjusted_num_gpvisits <- (patient_event_history_update$num_gpvisits*num_gpvisits_country/num_gpvisits_observed)*patient_event_history_update$diff_ANLYEAR
  
  
  num_specvisits_observed <- 0.545
  num_specvisits_country  <- maintenance_costs[2,1]
  patient_event_history_update$adjusted_num_specvisits <- (patient_event_history_update$num_specvisits*num_specvisits_country/num_specvisits_observed)*patient_event_history_update$diff_ANLYEAR
  
  
  ### Add the costs per visit
  
  #Dynagito and base case as well
  patient_event_history_update$maintenance_gpvisits_cost_hc              <- patient_event_history_update$adjusted_num_gpvisits*maintenance_costs[1,2]
  patient_event_history_update$maintenance_gpvisits_cost_hc_discounted   <- patient_event_history_update$maintenance_gpvisits_cost_hc/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  patient_event_history_update$maintenance_specvisits_cost_hc            <- patient_event_history_update$adjusted_num_specvisits*maintenance_costs[2,2]
  patient_event_history_update$maintenance_specvisits_cost_hc_discounted <- patient_event_history_update$maintenance_specvisits_cost_hc/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  
  patient_event_history_update$maintenance_gpvisits_cost_societal              <- patient_event_history_update$adjusted_num_gpvisits*maintenance_costs[1,3]
  patient_event_history_update$maintenance_gpvisits_cost_societal_discounted   <- patient_event_history_update$maintenance_gpvisits_cost_societal/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  patient_event_history_update$maintenance_specvisits_cost_societal            <- patient_event_history_update$adjusted_num_specvisits*maintenance_costs[2,3]
  patient_event_history_update$maintenance_specvisits_cost_societal_discounted <- patient_event_history_update$maintenance_specvisits_cost_societal/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  
  
  ### Now aggregate per SIMID: so what we get is the total number of visits during the lifetime
  
  ### Number of gpvisits
  num_gpvisits_patient <- aggregate(patient_event_history_update$adjusted_num_gpvisits,list(SIMID=patient_event_history_update$SIMID),sum)
  
  ### Number of specialist visits
  num_specvisits_patient <- aggregate(patient_event_history_update$adjusted_num_specvisits,list(SIMID=patient_event_history_update$SIMID),sum)
  
  
  ### Maintenance costs gpvisits
  gpvisits_costs_patient_hc            <- aggregate(patient_event_history_update$maintenance_gpvisits_cost_hc,list(SIMID=patient_event_history_update$SIMID),sum)
  gpvisits_costs_patient_hc_discounted <- aggregate(patient_event_history_update$maintenance_gpvisits_cost_hc_discounted,list(SIMID=patient_event_history_update$SIMID),sum)
  
  gpvisits_costs_patient_societal            <- aggregate(patient_event_history_update$maintenance_gpvisits_cost_societal,list(SIMID=patient_event_history_update$SIMID),sum)
  gpvisits_costs_patient_societal_discounted <- aggregate(patient_event_history_update$maintenance_gpvisits_cost_societal_discounted,list(SIMID=patient_event_history_update$SIMID),sum)
  
  ### Maintenance costs spec visits
  specvisits_costs_patient_hc            <- aggregate(patient_event_history_update$maintenance_specvisits_cost_hc,list(SIMID=patient_event_history_update$SIMID),sum)
  specvisits_costs_patient_hc_discounted <- aggregate(patient_event_history_update$maintenance_specvisits_cost_hc_discounted,list(SIMID=patient_event_history_update$SIMID),sum)
  
  specvisits_costs_patient_societal            <- aggregate(patient_event_history_update$maintenance_specvisits_cost_societal,list(SIMID=patient_event_history_update$SIMID),sum)
  specvisits_costs_patient_societal_discounted <- aggregate(patient_event_history_update$maintenance_specvisits_cost_societal_discounted,list(SIMID=patient_event_history_update$SIMID),sum)
  
  
  ### Spirometries costs
  patient_event_history_update$maintenance_spirometries_cost_hc            <- patient_event_history_update$diff_ANLYEAR*maintenance_costs[3,1]*maintenance_costs[3,2]
  patient_event_history_update$maintenance_spirometries_cost_hc_discounted <- patient_event_history_update$maintenance_spirometries_cost_hc/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR 
  
  patient_event_history_update$maintenance_spirometries_cost_societal            <- patient_event_history_update$diff_ANLYEAR*maintenance_costs[3,1]*maintenance_costs[3,3]
  patient_event_history_update$maintenance_spirometries_cost_societal_discounted <- patient_event_history_update$maintenance_spirometries_cost_societal/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR 
  
  
  ### influenza vaccination costs
  patient_event_history_update$maintenance_influenza_vaccination_cost_hc            <- patient_event_history_update$diff_ANLYEAR*maintenance_costs[4,1]*maintenance_costs[4,2]
  patient_event_history_update$maintenance_influenza_vaccination_cost_hc_discounted <- patient_event_history_update$maintenance_influenza_vaccination_cost_hc/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR  
  
  patient_event_history_update$maintenance_influenza_vaccination_cost_societal            <- patient_event_history_update$diff_ANLYEAR*maintenance_costs[4,1]*maintenance_costs[4,3]
  patient_event_history_update$maintenance_influenza_vaccination_cost_societal_discounted <- patient_event_history_update$maintenance_influenza_vaccination_cost_societal/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR  
  
  
  ### For ICS costs we need to check the value at baseline
  patient_event_history_update$maintenance_ics_cost_hc            <- patient_event_history_update_maintenance_costs$ICS*patient_event_history_update$diff_ANLYEAR*maintenance_costs[5,1]*maintenance_costs[5,2]
  patient_event_history_update$maintenance_ics_cost_hc_discounted <- patient_event_history_update$maintenance_ics_cost_hc/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  
  patient_event_history_update$maintenance_ics_cost_societal            <- patient_event_history_update_maintenance_costs$ICS*patient_event_history_update$diff_ANLYEAR*maintenance_costs[5,1]*maintenance_costs[5,3]
  patient_event_history_update$maintenance_ics_cost_societal_discounted <- patient_event_history_update$maintenance_ics_cost_societal/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  
  
  ### Travel costs depend on the number of visits calculated above
  perspective <- "Societal"
  
  patient_event_history_update$maintenance_distance_pcc_cost_societal            <- if(perspective=="Societal"){2*maintenance_costs[6,1]*maintenance_costs[6,2]*patient_event_history_update$adjusted_num_gpvisits}else{rep(0,row(patient_event_history_update))}
  patient_event_history_update$maintenance_distance_pcc_cost_societal_discounted <- patient_event_history_update$maintenance_distance_pcc_cost_societal/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  patient_event_history_update$maintenance_distance_sc_cost_societal             <- if(perspective=="Societal"){2*maintenance_costs[7,1]*maintenance_costs[7,2]*patient_event_history_update$adjusted_num_specvisits}else{rep(0,row(patient_event_history_update))}
  patient_event_history_update$maintenance_distance_sc_cost_societal_discounted  <- patient_event_history_update$maintenance_distance_sc_cost_societal/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  
  
  # health care perspective
  perspective <- "Health care"
  
  patient_event_history_update$maintenance_distance_pcc_cost_hc            <- if(perspective=="Societal"){2*maintenance_costs[6,1]*maintenance_costs[6,2]*patient_event_history_update$adjusted_num_gpvisits}else{rep(0,nrow(patient_event_history_update))}
  patient_event_history_update$maintenance_distance_pcc_cost_hc_discounted <- patient_event_history_update$maintenance_distance_pcc_cost_hc/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  patient_event_history_update$maintenance_distance_sc_cost_hc             <- if(perspective=="Societal"){2*maintenance_costs[7,1]*maintenance_costs[7,2]*patient_event_history_update$adjusted_num_specvisits}else{rep(0,nrow(patient_event_history_update))}
  patient_event_history_update$maintenance_distance_sc_cost_hc_discounted  <- patient_event_history_update$maintenance_distance_sc_cost_hc/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  
  
  ### Total maintenance costs 
  
  #Total maintenance costs: societal 
  patient_event_history_update$maintenance_costs_total_societal            <- patient_event_history_update$maintenance_gpvisits_cost_societal + patient_event_history_update$maintenance_specvisits_cost_societal + patient_event_history_update$maintenance_spirometries_cost_societal + patient_event_history_update$maintenance_influenza_vaccination_cost_societal + patient_event_history_update$maintenance_ics_cost_societal + patient_event_history_update$maintenance_distance_pcc_cost_societal + patient_event_history_update$maintenance_distance_sc_cost_societal
  patient_event_history_update$maintenance_costs_total_societal_discounted <- patient_event_history_update$maintenance_gpvisits_cost_societal_discounted + patient_event_history_update$maintenance_specvisits_cost_societal_discounted + patient_event_history_update$maintenance_spirometries_cost_societal_discounted + patient_event_history_update$maintenance_influenza_vaccination_cost_societal_discounted + patient_event_history_update$maintenance_ics_cost_societal_discounted + patient_event_history_update$maintenance_distance_pcc_cost_societal_discounted + patient_event_history_update$maintenance_distance_sc_cost_societal_discounted
  
  # Maintenance costs aggregated per patient and average 
  maintenance_costs_total_societal_patient            <- aggregate(patient_event_history_update$maintenance_costs_total_societal,list(SIMID=patient_event_history_update$SIMID),sum)
  maintenance_costs_total_societal_discounted_patient <- aggregate(patient_event_history_update$maintenance_costs_total_societal_discounted,list(SIMID=patient_event_history_update$SIMID),sum)
  
  
  #Total maintenance costs: hc 
  patient_event_history_update$maintenance_costs_total_hc            <- patient_event_history_update$maintenance_gpvisits_cost_hc + patient_event_history_update$maintenance_specvisits_cost_hc + patient_event_history_update$maintenance_spirometries_cost_hc + patient_event_history_update$maintenance_influenza_vaccination_cost_hc + patient_event_history_update$maintenance_ics_cost_hc + patient_event_history_update$maintenance_distance_pcc_cost_hc + patient_event_history_update$maintenance_distance_sc_cost_hc
  patient_event_history_update$maintenance_costs_total_hc_discounted <- patient_event_history_update$maintenance_gpvisits_cost_hc_discounted + patient_event_history_update$maintenance_specvisits_cost_hc_discounted + patient_event_history_update$maintenance_spirometries_cost_hc_discounted + patient_event_history_update$maintenance_influenza_vaccination_cost_hc_discounted + patient_event_history_update$maintenance_ics_cost_hc_discounted + patient_event_history_update$maintenance_distance_pcc_cost_hc_discounted + patient_event_history_update$maintenance_distance_sc_cost_hc_discounted
  # Maintenance costs aggregated per patient and average 
  maintenance_costs_total_hc_patient            <- aggregate(patient_event_history_update$maintenance_costs_total_hc,list(SIMID=patient_event_history_update$SIMID),sum)
  maintenance_costs_total_hc_discounted_patient <- aggregate(patient_event_history_update$maintenance_costs_total_hc_discounted,list(SIMID=patient_event_history_update$SIMID),sum)
  
  
  ## Side effects costs 
  
  #The assumption here is that the costs associated to pneumonia + hospitalization = cost severe exacerbation.
  #Cost pneumonia - hospitalization = cost moderate exacerbation.
  
  pneu_costs_calc_sim <- function(index){
    
    if(patient_event_history_update[index,"pneu_yn"]==1){
      if(patient_event_history_update[index,"pneu_hosp_yn"]==1){
        sev_exa_cost(patient_event_history_update[index,"age_time"])
      }else{
        mod_exa_cost(patient_event_history_update[index,"age_time"])
      }
    }else{0}
    
  }
  
  perspective <- "Societal"
  
  ### These are the simulated exacerbation costs
  patient_event_history_update$pneu_costs_societal            <- sapply(1:nrow(patient_event_history_update),pneu_costs_calc_sim)
  patient_event_history_update$pneu_costs_societal_discounted <- patient_event_history_update$pneu_costs_societal/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  
  ### Total pneumonia costs per patient (aggregate and then average) societal
  pneu_costs_societal_total_patient            <- aggregate(patient_event_history_update$pneu_costs_societal,list(SIMID=patient_event_history_update$SIMID),sum)
  pneu_costs_societal_total_patient_discounted <- aggregate(patient_event_history_update$pneu_costs_societal_discounted,list(SIMID=patient_event_history_update$SIMID),sum)
  
  
  
  perspective <- "Health care"
  ### These are the simulated exacerbation costs
  patient_event_history_update$pneu_costs_hc            <- sapply(1:nrow(patient_event_history_update),pneu_costs_calc_sim)
  patient_event_history_update$pneu_costs_hc_discounted <- patient_event_history_update$pneu_costs_hc/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  
  ### Total pneumonia costs per patient (aggregate and then average) hc
  pneu_costs_hc_total_patient            <- aggregate(patient_event_history_update$pneu_costs_hc,list(SIMID=patient_event_history_update$SIMID),sum)
  pneu_costs_hc_total_patient_discounted <- aggregate(patient_event_history_update$pneu_costs_hc_discounted,list(SIMID=patient_event_history_update$SIMID),sum)
  
  
  ### Total costs
  
  ### total societal undiscounted
  patient_event_history_update$total_costs_societal <- patient_event_history_update$treatment_costs_societal + patient_event_history_update$exacerbation_costs_societal + patient_event_history_update$maintenance_costs_total_societal + patient_event_history_update$pneu_costs_societal            
  
  ### Total costs per patient (aggregate and then average) societal
  total_costs_societal_patient <- aggregate(patient_event_history_update$total_costs_societal,list(SIMID=patient_event_history_update$SIMID),sum)
  
  ### total societal discounted
  patient_event_history_update$total_costs_societal_discounted <- patient_event_history_update$treatment_costs_societal_discounted + patient_event_history_update$exacerbation_costs_societal_discounted + patient_event_history_update$maintenance_costs_total_societal_discounted + patient_event_history_update$pneu_costs_societal_discounted            
  
  ### Total costs per patient (aggregate and then average) societal
  total_costs_societal_discounted_patient <- aggregate(patient_event_history_update$total_costs_societal_discounted,list(SIMID=patient_event_history_update$SIMID),sum)
  
  
  ### total HC undiscounted
  patient_event_history_update$total_costs_hc <- patient_event_history_update$treatment_costs_hc + patient_event_history_update$exacerbation_costs_hc + patient_event_history_update$maintenance_costs_total_hc + patient_event_history_update$pneu_costs_hc            
  
  ### Total costs per patient (aggregate and then average) hc
  total_costs_hc_patient <- aggregate(patient_event_history_update$total_costs_hc,list(SIMID=patient_event_history_update$SIMID),sum)
  
  ### total HC discounted
  patient_event_history_update$total_costs_hc_discounted <- patient_event_history_update$treatment_costs_hc_discounted + patient_event_history_update$exacerbation_costs_hc_discounted + patient_event_history_update$maintenance_costs_total_hc_discounted + patient_event_history_update$pneu_costs_hc_discounted            
  
  ### Total costs per patient (aggregate and then average) hc
  total_costs_hc_discounted_patient <- aggregate(patient_event_history_update$total_costs_hc_discounted,list(SIMID=patient_event_history_update$SIMID),sum)
  
  
  ### summary of total costs
  
  mean_total_costs_societal      <- round(mean(total_costs_societal_patient$x),0)
  mean_total_costs_societal_disc <- round(mean(total_costs_societal_discounted_patient$x),0)
  
  mean_total_costs_hc <- round(mean(total_costs_hc_patient$x),0)
  mean_total_costs_hc_disc <- round(mean(total_costs_hc_discounted_patient$x),0)
  
  
  ####################################
  ### Save patient history outputs ###
  ####################################
  
  ### Think about how we want to do this and automathize it
  #write.csv(patient_event_history_update,"Model - simulation results/simulation_17dic_psa_base_case_predicted_2.csv")
  
  
  #############################
  ### Return model outcomes ###
  #############################
  
  return(list(mean_annual_fev1_decline       = mean_annual_fev1_decline,
              mean_life_expectancy           = mean_life_expectancy,
              mean_mod_exa_rate              = mean_mod_exa_rate,
              mean_sev_exa_rate              = mean_sev_exa_rate,
              mean_CWE_TOT_change            = mean_CWE_TOT_change,
              mean_SGACT_change              = mean_SGACT_change,
              mean_cough_rate                = mean_cough_rate,
              mean_time_cough                = mean_time_cough, # july 2018
              mean_breathless_rate           = mean_breathless_rate,
              mean_time_breath               = mean_time_breath, #july 2018
              mean_pneu_rate                 = mean_pneu_rate,
              mean_pneu_hosp_rate            = mean_pneu_hosp_rate,
              mean_SGTOT_change              = mean_SGTOT_change,
              mean_qalys                     = mean_qalys,
              mean_qalys_disc                = mean_qalys_disc,
              mean_total_costs_societal      = mean_total_costs_societal,
              mean_total_costs_societal_disc = mean_total_costs_societal_disc,
              mean_total_costs_hc            = mean_total_costs_hc,
              mean_total_costs_hc_disc       = mean_total_costs_hc_disc))
  
  ### add cost components to return list
  
  
} #end COPD_model_simulation function



#####################
### Run the model ###
#####################

# 
# ### Call the model function with the appropriate inputs in the correct order
# # 1. Select whether the model will be based on "observed" (1) or "predicted" (0) regression equations
# # 2. Select the number of patients in the simulation
# 
# 
# COPD_simulation_deterministic_results <- COPD_model_simulation(0, # 0 = predicted equations / 1 = observed equations
#                                                                500, # Patient size
#                                                                0, # 0 = deterministic run / 1 = probabilistic run
#                                                                1, # Treatment effect: factor on time to exacerbation (report 1.3)
#                                                                1, # Treatment effect: factor on probability of severe exacerbation (report NA)
#                                                                1, # Treatment effect: factor on fev1 change (report 0.8)
#                                                                1, # Treatment effect: factor on cwe score change (report 1.2)
#                                                                0, # Treatment effect: absolute change in SGRQ activity score (report -4)
#                                                                1, # Treatment effect: factor on probability of cough/sputum (report 0.8)
#                                                                1, # Treatment effect: factor on probability of shortness of breath (report 0.8)
#                                                                0, # Treatment effect: absolute change in SGRQ total score (report -4)
#                                                                "", # Subgroup
#                                                                177) # random seed input: for the report I used seed =177
# 
# 
# COPD_simulation_deterministic_results
# 
# # For the PE paper
# # Run all scenarios with 500 patients
# # 1. base case
# # 2. scenario tx effect tte factor 1.15, medication cost 1.50 and same seed (=177) for patients
# # 3. scenario tx effect tte factor 1.15, medication cost 1.50 and different seed (=24) for patients
# # 4. scenario tx effect tte factor 1.15, medication cost 1.50, same seed (=177) for patients and no seeds for events (including life expectancy)
# # 5. base case not adjusting for remaining life expectancy
# 
# 
# 
# #write.csv(COPD_simulation_deterministic_results, "PE paper/simulation_scenario4_PE_paper.csv")



#######################
### PSA SIMULATION  ###
#######################

### This is basically calling the simulation function multiple times and getting the average results

COPD_model_PSA <- function(psa_size_input,
                           run_obs_input,
                           patient_size_input,
                           exac_treatment_effect_tte_input,
                           exac_treatment_effect_sevexa_input,
                           fev1_treatment_effect_input,
                           cwe_treatment_effect_input,
                           sgact_treatment_effect_input,
                           coughsputum_treatment_effect_input,
                           breathless_treatment_effect_input,
                           sgtot_treatment_effect_input,
                           subgroup_input){
  
  psa_history <- matrix(nrow=psa_size_input,ncol=19) #ncol is th enumber of variables we want to save. We fix this
  colnames(psa_history) <- c("Annual FEV1 decline",
                             "Life expectancy",
                             "Moderate exacerbation rate",
                             "Severe exacerbation rate",
                             "Change in CWE score",
                             "Change in SGRQ act. score",
                             "Change in SGRQ total score",
                             "Cough/sputum rate",
                             "Cough/sputum time", # July 2018
                             "Shortness of breath rate",
                             "Shortness of breath time", # July 2018
                             "Pneumonia rate",
                             "Hospitalisation after pneu. rate",
                             "QALYs",
                             "QALYs (discounted)",
                             "Costs (societal)",
                             "Costs (societal discounted)",
                             "Costs (Health-care)",
                             "Costs (Health-care discounted)")
  
  for(i in 1:psa_size_input){
    
    print(i)
    
    current_psa <- COPD_model_simulation(run_obs_input,
                                         patient_size_input,
                                         1, #note 1 here hardcoded
                                         exac_treatment_effect_tte_input,
                                         exac_treatment_effect_sevexa_input,
                                         fev1_treatment_effect_input,
                                         cwe_treatment_effect_input,
                                         sgact_treatment_effect_input,
                                         coughsputum_treatment_effect_input,
                                         breathless_treatment_effect_input,
                                         sgtot_treatment_effect_input,
                                         subgroup_input,
                                         i)
    
    psa_history[i,] <- c(current_psa$mean_annual_fev1_decline,
                         current_psa$mean_life_expectancy,
                         current_psa$mean_mod_exa_rate,
                         current_psa$mean_sev_exa_rate,
                         current_psa$mean_CWE_TOT_change,
                         current_psa$mean_SGACT_change,
                         current_psa$mean_SGTOT_change,
                         current_psa$mean_cough_rate,
                         current_psa$mean_time_cough, # July 2018
                         current_psa$mean_breathless_rate,
                         current_psa$mean_time_breath, # July 2018
                         current_psa$mean_pneu_rate,
                         current_psa$mean_pneu_hosp_rate,
                         current_psa$mean_qalys,
                         current_psa$mean_qalys_disc,
                         current_psa$mean_total_costs_societal,
                         current_psa$mean_total_costs_societal_disc,
                         current_psa$mean_total_costs_hc,
                         current_psa$mean_total_costs_hc_disc)
    
  }
  
  return(list(psa_history=psa_history))
  
}



### Change the settings below in order to run different PSA scenarios

# PSA for base case tio
# PSA for comparator tte 1.15 cost 1.50

# PSA base case tio no random seed for PSA
# PSA for comparator tte 1.15 cost 1.50 no random seed for PSA

run_obs                      <- 0
psa_size                     <- 300
patient_size                 <- 50
treatment_effect_tte         <- 1
treatment_effect_sevexa      <- 1
treatment_effect_fev1        <- 1
cwe_treatment_effect         <- 1
sgact_treatment_effect       <- 0
coughsputum_treatment_effect <- 1
breathless_treatment_effect  <- 1
sgtot_treatment_effect       <- 0
subgroup                     <- " "

init_time <- proc.time()

COPD_model_PSA_output <- COPD_model_PSA(psa_size,
                                        run_obs,
                                        patient_size,
                                        treatment_effect_tte,
                                        treatment_effect_sevexa,
                                        treatment_effect_fev1,
                                        cwe_treatment_effect,
                                        sgact_treatment_effect,
                                        coughsputum_treatment_effect,
                                        breathless_treatment_effect,
                                        sgtot_treatment_effect,
                                        subgroup)



tot_time <- proc.time() - init_time #in seconds

(tot_time/60) #minutes

(tot_time/60)/60 #hours

PSA_summary_table <- data.frame(COPD_model_PSA_output$psa_history)

colnames(PSA_summary_table) <- c("Annual FEV1 decline",
                                 "Life expectancy",
                                 "Moderate exacerbation rate",
                                 "Severe exacerbation rate",
                                 "Change in CWE score",
                                 "Change in SGRQ act. score",
                                 "Change in SGRQ total score",
                                 "Cough/sputum rate",
                                 "Cough/sputum time",
                                 "Shortness of breath rate",
                                 "Shortness of breath time",
                                 "Pneumonia rate",
                                 "Hospitalisation after pneu. rate",
                                 "QALYs",
                                 "QALYs (discounted)",
                                 "Costs (societal)",
                                 "Costs (societal discounted)",
                                 "Costs (Health-care)",
                                 "Costs (Health-care discounted)")

PSA_summary_table <- rbind(PSA_summary_table,colMeans(PSA_summary_table))

PSA_summary_table <- rbind(PSA_summary_table,lapply(PSA_summary_table[1:psa_size,], quantile, probs=c(0.025,0.975), name=FALSE))

rownames(PSA_summary_table) <- c(1:psa_size,"Mean","2.5%", "97.5%")

PSA_summary_table

write.csv(PSA_summary_table, paste("PE paper/",Sys.Date(),"_",psa_size,"x",patient_size,"_base_case_no_seed_PE_paper",'.csv', sep=''))

# 
# # # ############
# # # ## CEAC   ##
# # # ############
# # # 
# # # # This is one way of computing the CEAC: with the INMB (only valid for 2 treatments)
# # # 
# # # tio_psa_data  <- read.csv(paste0(wd,c("/Model - simulation results/PSA/21nov_20x10_dynagito_tio.csv")),sep=",")
# # # eff_psa_data  <- read.csv(paste0(wd,c("/Model - simulation results/PSA/21nov_20x10_dynagito_tio_olo.csv")),sep=",")
# # # 
# # # psa_size <- 20
# # # diffutil <- eff_psa_data[1:psa_size,14] - tio_psa_data[1:psa_size,14] #mind the bottom rows!!!!
# # # diffcost <- eff_psa_data[1:psa_size,18] - tio_psa_data[1:psa_size,18]
# # # 
# # # 
# # # k <- 0
# # # maxthreshold <- 80000
# # # 
# # # thresholds  <- c()
# # # CEAC_points <- c()
# # # 
# # # while (k < maxthreshold){
# # # 
# # #   inmb <- k*diffutil - diffcost
# # # 
# # #   thresholds <-  append(thresholds,k)
# # #   CEAC_points <- append(CEAC_points, length(inmb[inmb>0])/psa_size)
# # # 
# # #   k <- k+5
# # # 
# # # }
# # # 
# # # 
# # # # then we can plot the CEAC as follows
# # # #jpeg('C:/Users/Isaac/Dropbox/COPD/Model - simulation results/PSA/CEAC_100x100_base_case_vs_tte_scenario.jpg',width = 480, height = 480, units = "px", pointsize = 12,quality = 1000)
# # # windows()
# # # plot(thresholds, CEAC_points, pch=1, xlab=expression(lambda), ylab="Probability that a treatment is optimal",
# # #      col="black",ylim=c(0,1.15),type="l", las=1, lty=1,cex.main=1.4,cex.lab=1.4)
# # # title(main="Acceptability curves",cex.main=1.4)
# # # points(thresholds, 1-CEAC_points, pch=20, col="black",ylim=c(0,1.15),type="l",lty=2)
# # # #dev.off()
# # # 
# # # ##############
# # # # CE plane   #
# # # ##############
# # # 
# # # ### For the moment it is not incremental!!!
# # # 
# # # limy <- max(abs(min(diffcost)), abs(max(diffcost)))
# # # limx <- max(abs(min(diffutil)), abs(max(diffutil)))
# # # 
# # # 
# # # jpeg('Model - simulation results/PSA/CE_plane_100x100_base_case_vs_tte_scenario.jpg',width = 480, height = 480, units = "px", pointsize = 12,quality = 1000)
# # # # and now we plot
# # # windows()
# # # plot(diffutil,diffcost,
# # #      type="p",pch=18,xlab=expression(paste(Delta,Q)),ylab=expression(paste(Delta,"C(£)")),
# # #      #type="p",pch=18,xlab="QALYs",ylab="Costs (£)",
# # #      bty="l",las=1, xlim=c(-limx,limx),ylim=c(-limy,limy))
# # # 
# # # title(main="CE plane")
# # # abline(0,0)
# # # abline(0,0,v=0)
# # # points(mean(diffutil),mean(diffcost),col="green",pch=15)
# # # 
# # # dev.off()
# # # 
# # # 
# # # #points(mean(diffutil),mean(diffcost),col="green",pch=15)
# # # #text(0.1,-200, c("ICER=9870"),cex=1.4)
# # # #legend
# # # 
# # # 
# # # 
# # # # round(c(mean(COPD_model_PSA_output$psa_history[,1]),quantile(COPD_model_PSA_output$psa_history[,1],c(0.025,0.975))),4)
# # # # round(c(mean(COPD_model_PSA_output$psa_history[,2]),quantile(COPD_model_PSA_output$psa_history[,2],c(0.025,0.975))),4)
# # # #
# # # #
