########################################################################################
########## Patient-level discrete event simulation model for COPD ######################
########################################################################################

# Original work Copyright (C) 2019 Isaac Corro Ramos

# Before the simulation code starts, make sure that all the packages below are installed in your computer.
# Then load the packages.
library(lattice)
library(MASS)
library(survival)
library(plyr)

# When you are reading external files and exporting results you may set a working directory.
# This can be for example the folder where you saved some previous results.
# Do not forget to change this path into your personal working directory.
wd <- setwd("C:/Users/Isaac/Dropbox/COPD")

# This model makes use of previously estimated regression equations to predict time to events and COPD-related intermediate outcomes.
# These regression equations were estimated somewhere else and the results (e.g. the regression coefficients) can be 
# found in the paper by Hoogendoorn et al. (see README file).
# The simulation model reads these coefficients from previously saved .csv files. With these regression coefficients
# and patient characteristics (that will be read from another file), the functions below can be used to predict time to events 
# and COPD-related intermediate outcomes. 
# You can choose to read these functions from a different R file. In that case, you would need a statement like this:
# source("COPD model simulation - auxiliar functions.R"). 
# However, in this case, we decided to show all functions in the same file.

# Notation: as a general rule, names for function input parameters end with "_input".

# Three events can occur in the simulation: death, exacerbation and pneumonia.

### TIME TO DEATH (AT BASELINE)
# Time to death (at baseline) is a function of the previously estimated regression coefficients and patient characteristics.
# It is assumed to follow a Weibull distribution.

predicted_mortality_weibull <- function(regression_coefficents_mortality_weibull_input, 
                                        patient_characteristics_mortality_weibull_input){
  patient_characteristics_mortality_weibull_input <- as.numeric(patient_characteristics_mortality_weibull_input)
  ### the last element in the coefficients inputs is Log(scale) in survreg.
  ### Thus, exp(tail(mortality_weibull_regression_coef$Value,n=1)) is scale in survreg.
  ### BUT scale in survreg is 1/shape in rweibull.
  shape_mortality_weibull <- 1/exp(tail(regression_coefficents_mortality_weibull_input,n=1))
  scale_mortality_weibull <- exp(sum(head(regression_coefficents_mortality_weibull_input,n=-1)*c(1,patient_characteristics_mortality_weibull_input)) )
  return(list(shape_mortality_weibull=shape_mortality_weibull,
              scale_mortality_weibull=scale_mortality_weibull))
}

### TIME TO EXACERBATION
# Time to exacerbation is a function of the previously estimated regression coefficients and patient characteristics.
# It is assumed to follow a Weibull distribution.

predicted_exacerbation_weibull <- function(regression_coefficents_exacerbation_weibull_input, 
                                           patient_characteristics_exacerbation_weibull_input){
  patient_characteristics_exacerbation_weibull_input <- as.numeric(patient_characteristics_exacerbation_weibull_input)
  ### the last element in the coefficients inputs is Log(scale) in survreg.
  ### Thus, exp(tail(exacerbation_weibull_regression_coef$Value,n=1)) is scale in survreg.
  ### BUT scale in survreg is 1/shape in rweibull.
  shape_exacerbation_weibull <- 1/exp(tail(regression_coefficents_exacerbation_weibull_input,n=1))
  scale_exacerbation_weibull <- exp(sum(head(regression_coefficents_exacerbation_weibull_input,n=-1)*c(1,patient_characteristics_exacerbation_weibull_input)) )
  return(list(shape_exacerbation_weibull=shape_exacerbation_weibull,
              scale_exacerbation_weibull=scale_exacerbation_weibull))
}

### EXACERBATION SEVERITY
# The probability of experiencing a severe exacerbation is a function of the previously estimated regression coefficients 
# and patient characteristics. The function return also the log odds.  

predicted_exacerbation_severity <- function(regression_coefficents_exacerbation_severity_input, 
                                            patient_characteristics_exacerbation_severity_input){
  patient_characteristics_exacerbation_severity_input <- as.numeric(patient_characteristics_exacerbation_severity_input)
  log.oods.exacerbation_severity <- sum(regression_coefficents_exacerbation_severity_input*c(1,patient_characteristics_exacerbation_severity_input)) # this is log(ODDS)
  p.exacerbation_severity <- exp(log.oods.exacerbation_severity)/(1+exp(log.oods.exacerbation_severity))
  return(list(log.oods.exacerbation_severity=log.oods.exacerbation_severity,p.exacerbation_severity=p.exacerbation_severity))
}

### TIME TO PNEUMONIA
# Time to pneumonia is a function of the previously estimated regression coefficients and patient characteristics.
# It is assumed to follow a Weibull distribution.

predicted_pneumonia_weibull <- function(regression_coefficents_pneumonia_weibull_input, 
                                        patient_characteristics_pneumonia_weibull_input){
  patient_characteristics_pneumonia_weibull_input <- as.numeric(patient_characteristics_pneumonia_weibull_input)
  ### the last element in the coefficients inputs is Log(scale) in survreg.
  ### Thus, exp(tail(pneumonia_weibull_regression_coef$Value,n=1)) is scale in survreg.
  ### BUT scale in survreg is 1/shape in rweibull.
  shape_pneumonia_weibull <- 1/exp(tail(regression_coefficents_pneumonia_weibull_input,n=1))
  scale_pneumonia_weibull <- exp(sum(head(regression_coefficents_pneumonia_weibull_input,n=-1)*c(1,patient_characteristics_pneumonia_weibull_input)) )
  return(list(shape_pneumonia_weibull=shape_pneumonia_weibull,
              scale_pneumonia_weibull=scale_pneumonia_weibull))
}

### PNEUMONIA SEVERITY
# The probability of experiencing a pneumonia leading to hospitalisation is a function of the previously estimated regression 
# coefficients and patient characteristics. The function return also the log odds.  

predicted_pneumonia_hosp <- function(regression_coefficents_pneumonia_hosp_input, 
                                     patient_characteristics_pneumonia_hosp_input){
  patient_characteristics_pneumonia_hosp_input <- as.numeric(patient_characteristics_pneumonia_hosp_input)
  log.oods.pneumonia.hosp <- sum(regression_coefficents_pneumonia_hosp_input*c(1,patient_characteristics_pneumonia_hosp_input)) # this is log(ODDS)
  p.pneumonia.hosp        <- exp(log.oods.pneumonia.hosp)/(1+exp(log.oods.pneumonia.hosp))
  return(list(log.oods.pneumonia.hosp=log.oods.pneumonia.hosp,p.pneumonia.hosp=p.pneumonia.hosp)) 
}

# Six intermediate outcomes are included in the simulation: lung function, exercise capacity, symptoms (shortness of breath 
# and cough/sputum), physical activity and disease-specific quality of life.

### LUNG FUNCTION
# Lung function defined as FEV1 is predicted as a function of the previously estimated regression coefficients, 
# patient characteristics and possibly a treatmemt effect.

predicted_fev1 <- function(regression_coefficents_fev1_input, 
                           patient_characteristics_fev1_input,
                           fev1_treatment_effect_input){
  patient_characteristics_fev1 <- as.numeric(patient_characteristics_fev1_input)
  x <- c(1,patient_characteristics_fev1[1]*fev1_treatment_effect_input, #this is anlyear
         patient_characteristics_fev1[17], #this is FEVA_BL
         patient_characteristics_fev1[18], #this is modexac in patient characteristics
         patient_characteristics_fev1[19], #this is sevexac in patient characteristics
         patient_characteristics_fev1[1]*fev1_treatment_effect_input*patient_characteristics_fev1[2:17])
  fev_1 <- sum(regression_coefficents_fev1_input*x)
  return(list(fev_1=fev_1))
}

# FEVPPA is calculated from FEV1 and it is different for FEMALEs and males
#Equation 1B males (ECSC 1993 used by BI): if FEMALE=0 FEV1pred =  0.0430*htstd - 0.0290*AGE_TIME - 2.490.
#Equation 1B FEMALEs (ECSC 1993 used by BI):  if FEMALE=1 FEV1pred =  0.0395*htstd - 0.0250*AGE_TIME - 2.600.
#Equation 1C: FEVPPA = FEVA / FEV1pred *100

FEVPPA_calc <- function(FEMALE_input,height_input,age_input,feva_input){
  if(FEMALE_input==0){FEV1_pred <- 0.0430*height_input-0.0290*age_input-2.490}
  if(FEMALE_input==1){FEV1_pred <- 0.0395*height_input-0.0250*age_input-2.600}
  else{FEV1_pred <- FEMALE_input*(0.0395*height_input-0.0250*age_input-2.600)+(1-FEMALE_input)*(0.0430*height_input-0.0290*age_input-2.49)}
  FEVPPA <- max(0,100*(feva_input/FEV1_pred)) # Added max(0,x) to avoid negative numbers
  return(list(FEVPPA=FEVPPA,FEV1_pred=FEV1_pred))
}

### EXERCISE CAPACITY 
# Estimated as continuous exercise capacity (defined as treadmill test in seconds).

predicted_cwe_tot <- function(regression_coefficents_cwe_tot_input, 
                              patient_characteristics_cwe_tot_input){
  patient_characteristics_cwe_tot_input <- as.numeric(patient_characteristics_cwe_tot_input)
  totexa <- max(tail(patient_characteristics_cwe_tot_input,2))
  cwe_tot <- sum(regression_coefficents_cwe_tot_input*c(1,head(patient_characteristics_cwe_tot_input,-2),totexa))
  return(list(cwe_tot=cwe_tot))
}

### SHORTNESS OF BREATH
# The probability of experiencing shortness of breath is a function of the previously estimated regression coefficients 
# and patient characteristics. The function return also the log odds.  

predicted_breathless <- function(regression_coefficents_breathless_input, 
                                 patient_characteristics_breathless_input){
  patient_characteristics_breathless_input <- as.numeric(patient_characteristics_breathless_input)
  log.oods.breathless <- sum(regression_coefficents_breathless_input*c(1,patient_characteristics_breathless_input)) 
  p.breathless <- exp(log.oods.breathless)/(1+exp(log.oods.breathless))
  return(list(log.oods.breathless=log.oods.breathless,p.breathless=p.breathless))
}

### COUGH/SPUTUM
# The probability of experiencing cough/sputum is a function of the previously estimated regression coefficients 
# and patient characteristics. The function return also the log odds.  

predicted_coughsputum <- function(regression_coefficents_coughsputum_input, 
                                  patient_characteristics_coughsputum_input){
  patient_characteristics_coughsputum_input <- as.numeric(patient_characteristics_coughsputum_input)
  log.oods.coughsputum <- sum(regression_coefficents_coughsputum_input*c(1,patient_characteristics_coughsputum_input)) # this is log(ODDS)
  p.coughsputum <- exp(log.oods.coughsputum)/(1+exp(log.oods.coughsputum))
  return(list(log.oods.coughsputum=log.oods.coughsputum,p.coughsputum=p.coughsputum))
}

### PHYSICAL ACTIVITY
# Predicted as SGRQ Activity score. It is a function of the previously estimated regression coefficients 
# and patient characteristics.

predicted_SGACT <- function(regression_coefficents_SGACT_input, 
                            patient_characteristics_SGACT_input){
  patient_characteristics_SGACT_input <- as.numeric(patient_characteristics_SGACT_input)
  SGACT <- sum(regression_coefficents_SGACT_input*c(1,patient_characteristics_SGACT_input))
  return(list(SGACT=SGACT))
}

### QUALITY OF LIFE
# Disease specific quality of life defined as SGRQ total score. It is a function of the previously estimated regression coefficients 
# and patient characteristics.

predicted_SGTOT <- function(regression_coefficents_SGTOT_input, 
                            patient_characteristics_SGTOT_input){
  patient_characteristics_SGTOT_input <- as.numeric(patient_characteristics_SGTOT_input)
  SGTOT <- sum(regression_coefficents_SGTOT_input*c(1,patient_characteristics_SGTOT_input))
  return(list(SGTOT=SGTOT))
}

# The main simulation starts below. The code is used to 1) simulate patients’ clinical history, 2) calculate costs and 
# 3) calculate QALYs. Patients’ clinical histories are simulated first and, based on these, costs and QALYs are subsequently 
# calculated. Note that this could be implemented as three independent functions but in this tutorial we decided to show everything
# as one larger function called COPD_model_simulation. If a probabilistic sensitivity analysis is conducted, this function
# is basically called multiple times.

# The input parameters of the COPD_model_simulation function are the following:
# 1. patient_size_input = number of patients included in the simulation.
# 2. run_PSA_input = runs the model in probabilistic mode. Otherwise, deterministic.

# Eight treatment effect parameters:
# 3. exac_treatment_effect_tte_input = variable to indicate increase (or decrease) in time to exacerbation. Default should be 1.
# 4. exac_treatment_effect_sevexaprob_input = variable to indicate increase (or decrease) in probability of experiencing a severe exacerbation. Default should be 1.
# 5. fev1_treatment_effect_input = variable to indicate a reduction or increase in lung function decline. Default should be 0.
# 6. cwe_treatment_effect_input = variable to indicate increase (or decrease) in exercise capacity. Default should be 0.
# 7. sgact_treatment_effect_input = variable to indicate increase (or decrease) in physical activity. Default should be 0. 
# 8. coughsputum_treatment_effect_input = variable to indicate increase (or decrease) in probability of experiencing cough/sputum. Default should be 1.
# 9. breathless_treatment_effect_input variable to indicate increase (or decrease) in probability of experiencing shortness of breath. Default should be 1.
# 10. sgtot_treatment_effect_input = variable to indicate increase (or decrease) in quality of life. Default should be 0.  

# Other parameters 
# 11. seed_input = A random seed is used to ensure consistency in the model results as explained in the MDM paper (see README file).
COPD_model_simulation <- function(patient_size_input,
                                  run_PSA_input,
                                  exac_treatment_effect_tte_input,
                                  exac_treatment_effect_sevexaprob_input,
                                  fev1_treatment_effect_input,
                                  cwe_treatment_effect_input,
                                  sgact_treatment_effect_input,
                                  coughsputum_treatment_effect_input,
                                  breathless_treatment_effect_input,
                                  sgtot_treatment_effect_input,
                                  seed_input){
  #############
  ### SETUP ###
  #############
  
  ### Read regression coefficients from csv files.
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
  
  ### Assign predictors (explanatory variables) for the regression equations used in the model (further details in the ViH or MDM papers - see README file)
  fev1_predictors        <- c("ANLYEAR","FEMALE","AGE","BMI_CLASS_2","BMI_CLASS_3","SMOKER","SMPKY","OTHER_CVD", "REVERSIBILITY","DIABETES","DEPRESSION","HEART_FAILURE","ASTHMA","EMPHDIA","EOS_yn","ICS","FEVA_BL","MODEXAC_yn","SEVEXAC_yn")
  cwe_tot_predictors     <- c("FEMALE","BMI_CLASS_2","BMI_CLASS_3","SMOKER","SMPKY","OTHER_CVD","REVERSIBILITY","DIABETES","DEPRESSION","HEART_FAILURE","ASTHMA","EMPHDIA","ICS", "EOS_yn","AGE_TIME","FEVPPA","lag_SGACT","lag_CWE_TOT","MODEXAC_yn","SEVEXAC_yn")
  breathless_predictors  <- c("FEMALE","BMI_CLASS_2","BMI_CLASS_3","SMOKER","SMPKY_SCALED","OTHER_CVD","REVERSIBILITY_SCALED","DIABETES","DEPRESSION","HEART_FAILURE","ASTHMA","EMPHDIA","ICS","EOS_yn","ANLYEAR_SCALED","AGE_SCALED","FEVPPA_SCALED","MODEXAC_yn","SEVEXAC_yn","SGACT_SCALED","CWE_TOT_SCALED","lag_BREATHLESS_yn") 
  coughsputum_predictors <- c("FEMALE","BMI_CLASS_2","BMI_CLASS_3","SMOKER","SMPKY_SCALED","OTHER_CVD","REVERSIBILITY_SCALED","DIABETES","DEPRESSION","HEART_FAILURE","ASTHMA","EMPHDIA","ICS", "EOS_yn","ANLYEAR_SCALED","AGE_SCALED","FEVPPA_SCALED","MODEXAC_yn","SEVEXAC_yn","SGACT_SCALED","CWE_TOT_SCALED","lag_COUGHSPUTUM_yn")
  SGACT_predictors       <- c("FEMALE","BMI_CLASS_2","BMI_CLASS_3","SMOKER","SMPKY","OTHER_CVD","REVERSIBILITY","DIABETES","DEPRESSION","HEART_FAILURE","ASTHMA","EMPHDIA","ICS","EOS_yn","ANLYEAR", "AGE","FEVPPA","lag_SGACT","CWE_TOT","lag_BREATHLESS_yn","lag_COUGHSPUTUM_yn","lag_SGTOT")
  SGTOT_predictors       <- c("FEMALE","BMI_CLASS_2","BMI_CLASS_3","SMOKER","SMPKY","OTHER_CVD","REVERSIBILITY","DIABETES","DEPRESSION","HEART_FAILURE","ASTHMA","EMPHDIA","ICS","EOS_yn","ANLYEAR", "AGE","FEVPPA","lag_SGTOT","MODEXAC_yn","SEVEXAC_yn","SGACT","CWE_TOT","BREATHLESS_yn","COUGHSPUTUM_yn","PNEU_yn") 
  mortality_predictors   <- c("FEMALE","AGE","BMI_CLASS_2","BMI_CLASS_3","SMOKER","SMPKY","OTHER_CVD","REVERSIBILITY","DIABETES","DEPRESSION","HEART_FAILURE","ASTHMA","EMPHDIA","ICS","EOS_yn","FEVPPA","PREV_SEVEXAC_yn","SGACT","CWE_TOT","BREATHLESS_yn","COUGHSPUTUM_yn","SGTOT")
  exacerbation_predictors <- c("FEMALE","BMI_CLASS_2","BMI_CLASS_3","SMOKER","SMPKY","OTHER_CVD","REVERSIBILITY","DIABETES","DEPRESSION","HEART_FAILURE","ASTHMA","EMPHDIA","ICS","EOS_yn","lag_FEVPPA","lag_SGACT","PREV_TOTEXAC_yn","PREV_SEVEXAC_yn","lag_SGTOT","AGE_TIME")
  exacerbation_severity_predictors <- c("FEMALE","BMI_CLASS_2","BMI_CLASS_3","SMOKER","SMPKY_SCALED","OTHER_CVD","REVERSIBILITY_SCALED","DIABETES","DEPRESSION","HEART_FAILURE","ASTHMA","EMPHDIA","ICS","EOS_yn","lag_FEVPPA_SCALED","lag_SGACT_SCALED","PREV_TOTEXAC_yn", "PREV_SEVEXAC_yn", "lag_SGTOT_SCALED","AGE_TIME_SCALED")
  pneumonia_predictors <- c("FEMALE","BMI_CLASS_2","BMI_CLASS_3","SMOKER","SMPKY","OTHER_CVD","REVERSIBILITY","DIABETES","DEPRESSION","HEART_FAILURE","ASTHMA","EMPHDIA","ICS","EOS_yn","AGE_TIME")
  pneumonia_hosp_predictors <- c("FEMALE", "BMI_CLASS_2","BMI_CLASS_3","SMOKER","SMPKY_SCALED","OTHER_CVD","REVERSIBILITY_SCALED","DIABETES","DEPRESSION","HEART_FAILURE","ASTHMA","EMPHDIA","ICS","EOS_yn","AGE_TIME_SCALED")
  
  ### Read the data with the patient characteristics and select the complete cases (no missing values in patient characteristics are allowed)
  baseline_characteristics_run <- read.csv(paste0(wd,c("/Model - datasets/baseline_characteristics_predicted_data.csv")),sep=",")
  complete_cases <- baseline_characteristics_run[colnames(baseline_characteristics_run)][complete.cases(baseline_characteristics_run[colnames(baseline_characteristics_run)]),]
  
  ### Indicate the patient characteristics that we will save during the simulation. 
  history_characteristics <- c("SIMID","PTID","ANLYEAR","AGE_TIME","FEVA","FEVPPA","SEVEXAC_yn","MODEXAC_yn","CWE_TOT","SGACT","SGTOT","COUGHSPUTUM_yn","BREATHLESS_yn","PNEU_yn","pneu_hosp_yn","dead",
                               "FEMALE","AGE","BMI_CLASS_2","BMI_CLASS_3","SMOKER","SMPKY","OTHER_CVD",
                               "REVERSIBILITY", "DIABETES","DEPRESSION","HEART_FAILURE","ASTHMA","EMPHDIA","EOS_yn","ICS","FEVA_BL","HTSTD",
                               "lag_SGACT","lag_CWE_TOT","lag_BREATHLESS_yn","lag_COUGHSPUTUM_yn","lag_SGTOT",
                               "SMPKY_SCALED","REVERSIBILITY_SCALED","ANLYEAR_SCALED","AGE_SCALED","FEVPPA_SCALED","SGACT_SCALED","CWE_TOT_SCALED")
  
  ##################################################
  ########## MAIN PART I: simulate events ##########
  ##################################################
  
  # Regression coefficients for PSA -- once per PSA iteration: if the model is run in probabilistic mode, the regression coefficients are randomly drawn from multivariate normal distributions
  # Treatment effect parameters are randomly drawn from uniform distributions, allowing a 10% variation from the deterministic value. 
  if(run_PSA_input == 1){
    lung_function_regression_coef         <- mvrnorm(1,lung_function_regression_coef,lung_function_cov_matrix[,-1])
    cwe_tot_regression_coef               <- mvrnorm(1,cwe_tot_regression_coef,cwe_tot_cov_matrix[,-1])
    SGACT_regression_coef                 <- mvrnorm(1,SGACT_regression_coef,SGACT_cov_matrix[,-1])
    SGTOT_regression_coef                 <- mvrnorm(1,SGTOT_regression_coef,SGTOT_cov_matrix[,-1])
    exacerbation_weibull_regression_coef  <- mvrnorm(1,exacerbation_weibull_regression_coef,exacerbation_weibull_cov_matrix[,-1])
    pneumonia_weibull_regression_coef     <- mvrnorm(1,pneumonia_weibull_regression_coef,pneumonia_weibull_cov_matrix[,-1])
    exacerbation_severity_regression_coef <- mvrnorm(1,exacerbation_severity_regression_coef,exacerbation_severity_cov_matrix[,-1])
    pneumonia_hosp_regression_coef        <- mvrnorm(1,pneumonia_hosp_regression_coef,pneumonia_hosp_cov_matrix[,-1])
    breathless_regression_coef            <- mvrnorm(1,breathless_regression_coef,breathless_cov_matrix[,-1])
    coughsputum_regression_coef           <- mvrnorm(1,coughsputum_regression_coef,coughsputum_cov_matrix[,-1])
    mortality_weibull_regression_coef     <- mvrnorm(1,mortality_weibull_regression_coef,mortality_weibull_cov_matrix[,-1]) 
    
    ### Treatment effect parameters must come at the end to ensure the same random seed per treatment arm  
    exac_treatment_effect_tte_input        <- runif(1,min(c(exac_treatment_effect_tte_input-(exac_treatment_effect_tte_input-1)*0.1,exac_treatment_effect_tte_input+(exac_treatment_effect_tte_input-1)*0.1)),max(c(exac_treatment_effect_tte_input-(exac_treatment_effect_tte_input-1)*0.1,exac_treatment_effect_tte_input+(exac_treatment_effect_tte_input-1)*0.1)))
    exac_treatment_effect_sevexaprob_input <- runif(1,min(c(exac_treatment_effect_sevexaprob_input-(exac_treatment_effect_sevexaprob_input-1)*0.1,exac_treatment_effect_sevexaprob_input+(exac_treatment_effect_sevexaprob_input-1)*0.1)),max(c(exac_treatment_effect_sevexaprob_input-(exac_treatment_effect_sevexaprob_input-1)*0.1,exac_treatment_effect_sevexaprob_input+(exac_treatment_effect_sevexaprob_input-1)*0.1)))
    fev1_treatment_effect_input            <- runif(1,min(c(fev1_treatment_effect_input-(fev1_treatment_effect_input-1)*0.1,fev1_treatment_effect_input+(fev1_treatment_effect_input-1)*0.1)),max(c(fev1_treatment_effect_input-(fev1_treatment_effect_input-1)*0.1,fev1_treatment_effect_input+(fev1_treatment_effect_input-1)*0.1)))
    } # end if regression coef PSA
  
  # Next we sample with replacement from the pool of patients. IMPORTANT: Set a random seed to be able to replicate the results.
  # This first seed is used to draw the same patients when the function is called multiple times: e.g. patients have to be the same per treatment arm.
  set.seed(seed_input) 
  simulation_baseline_patients <-  complete_cases[sample(nrow(complete_cases), patient_size_input, replace = TRUE), ]
  
  # Add a simulation ID variable to summarize results after the simulation is finished. This is needed because patients might be repeated in the simulation, and then aggregating results by patient ID would not be correct.
  SIMID <- rep("NA",patient_size_input)
  simulation_baseline_patients <- cbind(simulation_baseline_patients,SIMID)
  
  # Create the simulation patient history table (for now is just empty)
  simulation_patients_history <- simulation_baseline_patients[FALSE,c(history_characteristics)]
  
  # Choose the 1st patient who will enter the loop
  patient_index <- 1
  
  # Begin the loop on the simulation size (i.e. the number of patients you want to simulate)
  for(patient_index in 1:patient_size_input){
    
    # Tip: print the patient index to know how advanced is the simulation.
    if(run_PSA_input == 0){print(patient_index)}
    
    # Pick the current patient from those selected fomr baseline
    current_patient <- simulation_baseline_patients[patient_index,]
    current_patient$SIMID <- patient_index
    
    # Predict only continuous variables at baseline. Otherwise, we may get "strange" values after the first event is simulated.    
    current_patient$FEVA          <- max(0,predicted_fev1(lung_function_regression_coef,current_patient[fev1_predictors],fev1_treatment_effect_input)$fev_1) 
    baseline_FEVPPA_calc          <- FEVPPA_calc(current_patient$FEMALE,current_patient$HTSTD,current_patient$AGE_TIME,current_patient$FEVA)
    current_patient$FEVPPA        <- baseline_FEVPPA_calc$FEVPPA
    current_patient$FEV1pred      <- baseline_FEVPPA_calc$FEV1_pred
    current_patient$FEVPPA_SCALED <- (current_patient$FEVPPA - attr(scale(baseline_characteristics_run$FEVPPA),"scaled:center"))/attr(scale(baseline_characteristics_run$FEVPPA),"scaled:scale")
    # Baseline predicted CWE (min observed is 42. here for now we truncate at 0)
    current_patient$CWE_TOT        <- cwe_treatment_effect_input*max(0,predicted_cwe_tot(cwe_tot_regression_coef,current_patient[cwe_tot_predictors])$cwe_tot) #max(0,predicted_cwe_tot(cwe_tot_regression_coef$Value,current_patient[cwe_tot_predictors])$cwe_tot)
    current_patient$CWE_TOT_SCALED <- (current_patient$CWE_TOT - attr(scale(baseline_characteristics_run$CWE_TOT),"scaled:center"))/attr(scale(baseline_characteristics_run$CWE_TOT),"scaled:scale")
    current_patient$SGACT        <- min(100,max(0,sgact_treatment_effect_input+predicted_SGACT(SGACT_regression_coef,current_patient[SGACT_predictors])$SGACT))
    current_patient$SGACT_SCALED <- (current_patient$SGACT - attr(scale(baseline_characteristics_run$SGACT),"scaled:center"))/attr(scale(baseline_characteristics_run$SGACT),"scaled:scale")
    current_patient$SGTOT        <- min(100,max(0,sgtot_treatment_effect_input+predicted_SGTOT(SGTOT_regression_coef,current_patient[SGTOT_predictors])$SGTOT))
    current_patient$SGTOT_SCALED <- (current_patient$SGTOT - attr(scale(baseline_characteristics_run$SGTOT),"scaled:center"))/attr(scale(baseline_characteristics_run$SGTOT),"scaled:scale")
        
    # Save the characteristics to be used in the simulation history (not those that are stable, only those changing) 
    simulation_patients_history <- rbind(simulation_patients_history,current_patient[history_characteristics])
        
    ######################################
    # Sample from mortality distribution #
    ######################################
    
    # These random seeds are used to draw the life expectancy. They ensure that patients in different arms have consistent time to death.
    # For example, if no treatment effect has been applied, time to death should be the same.
    # For the PSA, use a different seed for each patient in the loop: e.g. loop size is 500x100 so 100 is number of patients.
    # This can be changed accordingly or just select a large enough value for the seed.
    if(run_PSA_input == 0){set.seed(current_patient$SIMID)}else{set.seed((seed_input*100)+current_patient$SIMID)} 
    
    # Mortality at baseline: Weibull distribution
    baseline_remaining_life_exp_parameters <- predicted_mortality_weibull(mortality_weibull_regression_coef,current_patient[mortality_predictors])
    baseline_remaining_life_exp            <- rweibull(1,baseline_remaining_life_exp_parameters$shape_mortality_weibull,baseline_remaining_life_exp_parameters$scale_mortality_weibull)/365
    baseline_remaining_life_exp_mean       <- baseline_remaining_life_exp_parameters$scale_mortality_weibull*gamma(1+1/baseline_remaining_life_exp_parameters$shape_mortality_weibull)/365
    
    # Save the sampled life expectancy. This will be used as reference for adjusting the remaining life expectancy after events (see MDM paper for details -- README file).
    current_remaining_life_exp  <- baseline_remaining_life_exp
    lag_current_mortality       <- baseline_remaining_life_exp_mean 
    
    #######################################################
    # Start the "timed" simulation (while loop = clock)   #
    #######################################################
    
    # Initialize the index for events 
    current_event   <- 1
    
    # Set random seeds for each time to exacerbation and penumonia. They ensure that patients in different arms have consistent time to events.
    # For example, if no treatment effect has been applied, time to event should be the same. 
    # We used 100 but can be anything large enough (no patient will experience 100 events in the simulation)
    factor_for_seed <- 100 
    
    # While loop starts for alive patients
    while(current_patient$dead==0){
      
      # Set new random seeds to draw the time to events while the patient is alive 
      if(run_PSA_input == 0){set.seed(factor_for_seed*current_patient$SIMID + current_event)}
      else{set.seed(factor_for_seed*seed_input*current_patient$SIMID + current_event)}
      
      ###########################################
      # Sample exacerbation and pneumonia time  #
      ###########################################
      
      # Exacerbation is Weibull
      current_exacerbation_parameters <- predicted_exacerbation_weibull(exacerbation_weibull_regression_coef,current_patient[exacerbation_predictors])
      current_exacerbation_time       <- rweibull(1,current_exacerbation_parameters$shape_exacerbation_weibull,current_exacerbation_parameters$scale_exacerbation_weibull)
      current_exacerbation_time       <- exac_treatment_effect_tte_input*current_exacerbation_time
      
      # Pneumonia is Weibull
      current_pneumonia_parameters <- predicted_pneumonia_weibull(pneumonia_weibull_regression_coef,current_patient[pneumonia_predictors])
      current_pneumonia_time       <- rweibull(1,current_pneumonia_parameters$shape_pneumonia_weibull,current_pneumonia_parameters$scale_pneumonia_weibull)
      
      ####################################
      # Update patient characteristics   #
      ####################################
      
      # We first copy all the previous characteristics
      current_patient_update <- current_patient
      
      # Update first the time lagged characteristics (these are needed for prediction of future events/intermediate characteristics)
      current_patient_update$lag_FEVPPA        <- current_patient$FEVPPA
      current_patient_update$lag_FEVPPA_SCALED <- (current_patient$FEVPPA - attr(scale(baseline_characteristics_run$FEVPPA),"scaled:center"))/attr(scale(baseline_characteristics_run$FEVPPA),"scaled:scale")
      current_patient_update$lag_SGACT         <- current_patient$SGACT
      current_patient_update$lag_SGACT_SCALED  <- (current_patient$SGACT - attr(scale(baseline_characteristics_run$SGACT),"scaled:center"))/attr(scale(baseline_characteristics_run$SGACT),"scaled:scale")
      current_patient_update$lag_CWE_TOT       <- current_patient$CWE_TOT
      current_patient_update$lag_BREATHLESS_yn  <- current_patient$BREATHLESS_yn
      current_patient_update$lag_COUGHSPUTUM_yn <- current_patient$COUGHSPUTUM_yn
      current_patient_update$lag_SGTOT         <- current_patient$SGTOT
      current_patient_update$lag_SGTOT_SCALED  <- (current_patient$SGTOT - attr(scale(baseline_characteristics_run$SGTOT),"scaled:center"))/attr(scale(baseline_characteristics_run$SGTOT),"scaled:scale")
      current_patient_update$PREV_SEVEXAC_yn <- current_patient$SEVEXAC_yn
      current_patient_update$PREV_TOTEXAC_yn <- ifelse(current_patient$SEVEXAC_yn==1 | current_patient$MODEXAC_yn==1,1,0)
            
      # Next update depending on the event occurred: if an exacerbation happened, then update at exacerbation time the following characteristics
      if(min(current_remaining_life_exp,current_exacerbation_time,current_pneumonia_time)==current_exacerbation_time){
        # Patient is still alive
        current_patient_update$dead <- 0
        # Update time-dependent characteristics
        current_patient_update$ANLYEAR         <- current_patient$ANLYEAR  + current_exacerbation_time
        current_patient_update$ANLYEAR_SCALED  <- (current_patient_update$ANLYEAR - 1.529411)/1.376549 # hard-coded based on dataset
        current_patient_update$lag_ANLYEAR     <- current_exacerbation_time
        current_patient_update$AGE_TIME        <- current_patient$AGE_TIME + current_exacerbation_time
        current_patient_update$AGE_TIME_SCALED <- (current_patient_update$AGE_TIME - attr(scale(baseline_characteristics_run$AGE_TIME),"scaled:center"))/attr(scale(baseline_characteristics_run$AGE_TIME),"scaled:scale")
        # Assumption: if age + current life expectancy is > 100 years then we force death
        if(current_patient_update$AGE_TIME>100){current_patient_update$dead <- 1}
        
        # Update exacerbation status. Decide first whether the exacerbation was moderate or severe.
        # Calculate first the probability of the exacerbation being severe
        current_severity_prob   <- exac_treatment_effect_sevexaprob_input*predicted_exacerbation_severity(exacerbation_severity_regression_coef,current_patient_update[exacerbation_severity_predictors])$p.exacerbation_severity
        # Then sample from a Bernoulli distirbution (severe = yes/no)
        current_severity_sample <- rbinom(1,1,current_severity_prob)
        # If the exacerbation was severe then update the following
        if(current_severity_sample == 1){
          current_patient_update$MODEXAC_yn <- 0
          current_patient_update$SEVEXAC_yn <- 1
          
          # Assumption: additional death risk because of severe exacerbation 
          current_sevexa_death_prob   <- 0.063 # hard-coded based on data. See ViH paper (README file)
          # Simulate whether patient dies from severe exacerbation
          current_sevexa_death_sample <- rbinom(1,1,current_sevexa_death_prob)
          if(current_sevexa_death_sample == 1){current_patient_update$dead <- 1}else{current_patient_update$dead <- 0}
        }        
        else{ # If the exacerbation was moderate, just update status
          current_patient_update$MODEXAC_yn <- 1
          current_patient_update$SEVEXAC_yn <- 0
        }
        # Update pneumonia status too 
        current_patient_update$PNEU_yn      <- 0
        current_patient_update$pneu_hosp_yn <- 0
      }
            
      # If a pneumonia happened, then update at pneumonia time the following characteristics:
      if(min(current_remaining_life_exp,current_exacerbation_time,current_pneumonia_time)==current_pneumonia_time){
        
        # Update time-dependent characteristics
        current_patient_update$ANLYEAR         <- current_patient$ANLYEAR  + current_pneumonia_time
        current_patient_update$ANLYEAR_SCALED  <- (current_patient_update$ANLYEAR - 1.529411)/1.376549 # hard-coded based on dataset
        current_patient_update$lag_ANLYEAR     <- current_pneumonia_time
        current_patient_update$AGE_TIME        <- current_patient$AGE_TIME + current_pneumonia_time
        current_patient_update$AGE_TIME_SCALED <- (current_patient_update$AGE_TIME - attr(scale(baseline_characteristics_run$AGE_TIME),"scaled:center"))/attr(scale(baseline_characteristics_run$AGE_TIME),"scaled:scale")
        
        # Force death at 100 years
        if(current_patient_update$AGE_TIME>100){current_patient_update$dead <- 1}
        
        # Update pneumonia status 
        current_patient_update$PNEU_yn      <- 1
        current_patient_update$pneu_hosp_yn <- rbinom(1,1,predicted_pneumonia_hosp(pneumonia_hosp_regression_coef,current_patient_update[pneumonia_hosp_predictors])$p.pneumonia.hosp)
        
        # Patient might die because of pneumonia after hospitalisation
        if(current_patient_update$pneu_hosp_yn==1){
          current_pneu_death_prob   <- 1607/19786 ### hard-coded based on dataset
          current_pneu_death_sample <- rbinom(1,1,current_pneu_death_prob)
          if(current_pneu_death_sample==1){current_patient_update$dead <- 1}else{current_patient_update$dead <- 0}
        }
      }
            
      # If death  happened then update age and finish the simulation for this patient
      if(min(current_remaining_life_exp,current_exacerbation_time)==current_remaining_life_exp){
        # Patient is dead
        current_patient_update$dead <- 1
        # Update time-relatedcharacteristics
        current_patient_update$ANLYEAR          <- current_patient$ANLYEAR  + current_remaining_life_exp
        current_patient_update$ANLYEAR_SCALED   <- (current_patient_update$ANLYEAR - 1.529411)/1.376549 # hardcoded based on dataset
        current_patient_update$lag_ANLYEAR      <- current_remaining_life_exp
        current_patient_update$AGE_TIME         <- current_patient$AGE_TIME + current_remaining_life_exp
        current_patient_update$AGE_TIME_SCALED  <- (current_patient_update$AGE_TIME - attr(scale(baseline_characteristics_run$AGE_TIME),"scaled:center"))/attr(scale(baseline_characteristics_run$AGE_TIME),"scaled:scale")
        # Update exacerbation and pneumonia status
        current_patient_update$MODEXAC_yn   <- 0
        current_patient_update$SEVEXAC_yn   <- 0
        current_patient_update$PNEU_yn      <- 0
        current_patient_update$pneu_hosp_yn <- 0
      }
            
      ################################################################
      # Update continuous variables depending on the event occurred  #
      ################################################################
      
      # Update FEV1, CWE, SGACT, symptoms and SGTOT. The order is important here: see ViH for further details (README file)
      # Update FEV1-related variables
      current_patient_update$FEVA          <- max(0,predicted_fev1(lung_function_regression_coef,current_patient_update[fev1_predictors],fev1_treatment_effect_input)$fev_1) 
      current_FEVPPA_calc                  <- FEVPPA_calc(current_patient_update$FEMALE,current_patient_update$HTSTD,current_patient_update$AGE_TIME,current_patient_update$FEVA)
      current_patient_update$FEVPPA        <- current_FEVPPA_calc$FEVPPA
      current_patient_update$FEV1pred      <- current_FEVPPA_calc$FEV1_pred
      current_patient_update$FEVPPA_SCALED <- (current_patient_update$FEVPPA - attr(scale(baseline_characteristics_run$FEVPPA),"scaled:center"))/attr(scale(baseline_characteristics_run$FEVPPA),"scaled:scale")
      # Force death if FEV1 < 0.2
      if(current_patient_update$FEVA < 0.2){current_patient_update$dead <- 1}
      
      # Update CWE (the minimum value observed in the dataset was 42. In the code we truncate at 0)
      current_patient_update$CWE_TOT        <- cwe_treatment_effect_input*max(0,predicted_cwe_tot(cwe_tot_regression_coef,current_patient_update[cwe_tot_predictors])$cwe_tot) #max(0,predicted_cwe_tot(cwe_tot_regression_coef$Value,current_patient_update[cwe_tot_predictors])$cwe_tot)
      current_patient_update$CWE_TOT_SCALED <- (current_patient_update$CWE_TOT - attr(scale(baseline_characteristics_run$CWE_TOT),"scaled:center"))/attr(scale(baseline_characteristics_run$CWE_TOT),"scaled:scale")
      
      # Update SGACT: range of values [0, 100]
      current_patient_update$SGACT        <- min(100,max(0,sgact_treatment_effect_input+predicted_SGACT(SGACT_regression_coef,current_patient_update[SGACT_predictors])$SGACT))
      current_patient_update$SGACT_SCALED <- (current_patient_update$SGACT - attr(scale(baseline_characteristics_run$SGACT),"scaled:center"))/attr(scale(baseline_characteristics_run$SGACT),"scaled:scale")
      
      # Update breathlesness and coughsputum status (sample from Bernoulli distributions)
      current_patient_update$BREATHLESS_yn  <- rbinom(1,1,breathless_treatment_effect_input*predicted_breathless(breathless_regression_coef,current_patient_update[breathless_predictors])$p.breathless)
      current_patient_update$COUGHSPUTUM_yn <- rbinom(1,1,coughsputum_treatment_effect_input*predicted_coughsputum(coughsputum_regression_coef,current_patient_update[coughsputum_predictors])$p.coughsputum)
      
      # Update SGTOT: range of values [0, 100]
      current_patient_update$SGTOT        <- min(100,max(0,sgtot_treatment_effect_input+predicted_SGTOT(SGTOT_regression_coef,current_patient_update[SGTOT_predictors])$SGTOT))
      current_patient_update$SGTOT_SCALED <- (current_patient_update$SGTOT - attr(scale(baseline_characteristics_run$SGTOT),"scaled:center"))/attr(scale(baseline_characteristics_run$SGTOT),"scaled:scale")
      
      # When all characteristics have been updated, we add these to the patient history maxtrix
      simulation_patients_history <- rbind(simulation_patients_history,current_patient_update[history_characteristics])
      
      # And update the current patient which will "go up" in the while loop again
      current_patient <- current_patient_update
            
      # Adjust remaining life expectancy according to updated patient characteristics. See MDM paper for details.
      current_mortality_parameters <- predicted_mortality_weibull(mortality_weibull_regression_coef,current_patient[mortality_predictors])
      current_mortality            <- current_mortality_parameters$scale_mortality_weibull*gamma(1+1/current_mortality_parameters$shape_mortality_weibull)/365
      ### Adjust remaining life expectancy according to improvement or worsened in condition with respect to baseline or previous time period
      current_remaining_life_exp <- max(0,(current_remaining_life_exp - current_patient$lag_ANLYEAR)*(current_mortality/lag_current_mortality))
      
      # Update mortality curve and event index
      lag_current_mortality    <- current_mortality 
      current_event            <- current_event + 1
    } #end while loop 
    # Move to another patient
    patient_index <- patient_index + 1
  } #end for loop in number of patients
  
  #################################################################################
  ########## MAIN PART II: Update intermediate outcomes after every year ##########
  #################################################################################
  
  # We were interested in updating the intermediate outcomes after every year because it is commonly done in COPD, e.g. annual decline in FEV1 is one of the most commonly used outcomes.
  # However, this may not be relevant for all models and, therefore, not needed. Removing this part from the model would make it simpler and faster since this part basically runs a loop on the previously simulated outcomes.
  # The simulated patient characteristics of interest are the following:
  patient_characteristics_saved <- c("SIMID","PTID","ANLYEAR","AGE_TIME","FEVA","FEVPPA","SEVEXAC_yn","MODEXAC_yn","SGACT","SGTOT","COUGHSPUTUM_yn","BREATHLESS_yn","PNEU_yn","pneu_hosp_yn","dead")
  
  # These characteristics are saved in a currently empty matrix
  patient_event_history_update <- simulation_patients_history[FALSE,c(patient_characteristics_saved)]
  
  # The above defined empt matrix will be filled-in in the loop below.
  for(i in 1:max(simulation_patients_history$SIMID)){
    
    # Loop counter for the PSA progress (you may delete this if not interested)
    if(run_PSA_input == 0){print(i+patient_size_input)}
    
    # Slect the first simulated patient and calculate how many years passed between simulated events
    current_patient_event_history <- simulation_patients_history[which(simulation_patients_history$SIMID == i),]
    where    <- tail(1:nrow(current_patient_event_history),-1) 
    how_many <- floor(diff(current_patient_event_history$ANLYEAR))
    
    # Then add to the history matrix one row for each year that passed between events
    current_patient_event_history_update <- data.frame(matrix(,ncol=ncol(current_patient_event_history),nrow=max(where)+sum(how_many)))
    colnames(current_patient_event_history_update) <- c(history_characteristics)
    current_patient_event_history_update[c(1,cumsum(how_many)+where),] <- current_patient_event_history[1:nrow(current_patient_event_history),]
    current_patient_event_history_update$SIMID <- current_patient_event_history$SIMID[1]
    current_patient_event_history_update$PTID  <- current_patient_event_history$PTID[1]
    
    # Then, for each year between events, re-calculate all the intemediate outcomes
    for(j in 2:(nrow(current_patient_event_history_update)-1)){
      
      if(is.na(current_patient_event_history_update[j,]$ANLYEAR)==TRUE){
        current_patient_event_history_update[j,]$FEMALE <- current_patient_event_history_update[j-1,]$FEMALE
        current_patient_event_history_update[j,]$AGE <- current_patient_event_history_update[j-1,]$AGE
        current_patient_event_history_update[j,]$AGE_SCALED <- current_patient_event_history_update[j-1,]$AGE_SCALED
        current_patient_event_history_update[j,]$BMI_CLASS_2 <- current_patient_event_history_update[j-1,]$BMI_CLASS_2  
        current_patient_event_history_update[j,]$BMI_CLASS_3 <- current_patient_event_history_update[j-1,]$BMI_CLASS_3
        current_patient_event_history_update[j,]$SMOKER <- current_patient_event_history_update[j-1,]$SMOKER
        current_patient_event_history_update[j,]$SMPKY <- current_patient_event_history_update[j-1,]$SMPKY
        current_patient_event_history_update[j,]$SMPKY_SCALED <- current_patient_event_history_update[j-1,]$SMPKY_SCALED
        current_patient_event_history_update[j,]$OTHER_CVD <- current_patient_event_history_update[j-1,]$OTHER_CVD  
        current_patient_event_history_update[j,]$REVERSIBILITY <- current_patient_event_history_update[j-1,]$REVERSIBILITY
        current_patient_event_history_update[j,]$REVERSIBILITY_SCALED <- current_patient_event_history_update[j-1,]$REVERSIBILITY_SCALED 
        current_patient_event_history_update[j,]$DIABETES <- current_patient_event_history_update[j-1,]$DIABETES
        current_patient_event_history_update[j,]$DEPRESSION <- current_patient_event_history_update[j-1,]$DEPRESSION
        current_patient_event_history_update[j,]$HEART_FAILURE <- current_patient_event_history_update[j-1,]$HEART_FAILURE
        current_patient_event_history_update[j,]$ASTHMA <- current_patient_event_history_update[j-1,]$ASTHMA
        current_patient_event_history_update[j,]$EMPHDIA <- current_patient_event_history_update[j-1,]$EMPHDIA
        current_patient_event_history_update[j,]$EOS_yn <- current_patient_event_history_update[j-1,]$EOS_yn
        current_patient_event_history_update[j,]$ICS <- current_patient_event_history_update[j-1,]$ICS
        current_patient_event_history_update[j,]$FEVA_BL <- current_patient_event_history_update[j-1,]$FEVA_BL
        current_patient_event_history_update[j,]$HTSTD <- current_patient_event_history_update[j-1,]$HTSTD
        current_patient_event_history_update[j,]$ANLYEAR <- current_patient_event_history_update[j-1,]$ANLYEAR + 1 
        current_patient_event_history_update[j,]$ANLYEAR_SCALED <- (current_patient_event_history_update[j,]$ANLYEAR - 1.529411)/1.376549
        current_patient_event_history_update[j,]$AGE_TIME <- current_patient_event_history_update[j-1,]$AGE_TIME + 1 
        current_patient_event_history_update[j,]$SEVEXAC_yn <- 0
        current_patient_event_history_update[j,]$MODEXAC_yn <- 0
        current_patient_event_history_update[j,]$PNEU_yn <- 0
        current_patient_event_history_update[j,]$pneu_hosp_yn <- 0
        current_patient_event_history_update[j,]$dead <- 0
        current_patient_event_history_update[j,]$FEVA <- max(0,predicted_fev1(lung_function_regression_coef,current_patient_event_history_update[j,][fev1_predictors],fev1_treatment_effect_input)$fev_1) 
        current_patient_event_history_update_FEVPPA_calc <- FEVPPA_calc(current_patient_event_history_update[j,]$FEMALE,current_patient_event_history_update[j,]$HTSTD,current_patient_event_history_update[j,]$AGE_TIME,current_patient_event_history_update[j,]$FEVA)
        current_patient_event_history_update[j,]$FEVPPA <- current_patient_event_history_update_FEVPPA_calc$FEVPPA
        current_patient_event_history_update[j,]$FEVPPA_SCALED <- (current_patient_event_history_update[j,]$FEVPPA - attr(scale(baseline_characteristics_run$FEVPPA),"scaled:center"))/attr(scale(baseline_characteristics_run$FEVPPA),"scaled:scale")
        current_patient_event_history_update[j,]$lag_SGACT <- current_patient_event_history_update[j-1,]$SGACT 
        current_patient_event_history_update[j,]$lag_CWE_TOT <- current_patient_event_history_update[j-1,]$CWE_TOT
        current_patient_event_history_update[j,]$CWE_TOT <- max(0,predicted_cwe_tot(cwe_tot_regression_coef,current_patient_event_history_update[j,][cwe_tot_predictors])$cwe_tot) #max(0,predicted_cwe_tot(cwe_tot_regression_coef$Value,current_patient_event_history_update[j,][cwe_tot_predictors])$cwe_tot)    
        current_patient_event_history_update[j,]$CWE_TOT_SCALED <- (current_patient_event_history_update[j,]$CWE_TOT - attr(scale(baseline_characteristics_run$CWE_TOT),"scaled:center"))/attr(scale(baseline_characteristics_run$CWE_TOT),"scaled:scale")
        current_patient_event_history_update[j,]$lag_BREATHLESS_yn  <- current_patient_event_history_update[j-1,]$BREATHLESS_yn
        current_patient_event_history_update[j,]$lag_COUGHSPUTUM_yn <- current_patient_event_history_update[j-1,]$COUGHSPUTUM_yn
        current_patient_event_history_update[j,]$lag_SGTOT <- current_patient_event_history_update[j-1,]$SGTOT
        current_patient_event_history_update[j,]$SGACT <- min(100,max(0,predicted_SGACT(SGACT_regression_coef,current_patient_event_history_update[j,][SGACT_predictors])$SGACT))
        current_patient_event_history_update[j,]$SGACT_SCALED <- (current_patient_event_history_update[j,]$SGACT - attr(scale(baseline_characteristics_run$SGACT),"scaled:center"))/attr(scale(baseline_characteristics_run$SGACT),"scaled:scale")
        current_patient_event_history_update[j,]$BREATHLESS_yn <- rbinom(1,1,rbinom(1,1,predicted_breathless(breathless_regression_coef,current_patient_event_history_update[j,][breathless_predictors])$p.breathless))
        current_patient_event_history_update[j,]$COUGHSPUTUM_yn <- rbinom(1,1,rbinom(1,1,predicted_coughsputum(coughsputum_regression_coef,current_patient_event_history_update[j,][coughsputum_predictors])$p.coughsputum))
        current_patient_event_history_update[j,]$SGTOT <- min(100,max(0,predicted_SGTOT(SGTOT_regression_coef,current_patient_event_history_update[j,][SGTOT_predictors])$SGTOT))
      }
    } #end for loop per patient
    # Update the simulated clinical history with all the annual intermediate outcomes 
    patient_event_history_update <- rbind(patient_event_history_update,current_patient_event_history_update[,c(patient_characteristics_saved)])
  }
  
  #################################################################
  ########## MAIN PART III: Calculate aggregated results ##########
  #################################################################
  
  # We first made additional columns for the "diff" variables, which will calculate the difference between two consecutive outcomes (at time t minus at time t-1)
  patient_event_history_update$diff_ANLYEAR <- "NA"
  patient_event_history_update$diff_FEVA    <- "NA"
  patient_event_history_update$diff_SGTOT   <- "NA"
  patient_event_history_update$diff_SGACT   <- "NA"
  patient_event_history_update$diff_CWE_TOT <- "NA"
  
  # Calculate the "diff" variables to calculate the change in outcomes per year
  diff_ANLYEAR <- ddply(patient_event_history_update, "SIMID", summarize, diff_ANLYEAR = c(0,diff(ANLYEAR)))
  diff_FEVA    <- ddply(patient_event_history_update, "SIMID", summarize, diff_FEVA    = c(0,diff(FEVA)))
  diff_SGTOT   <- ddply(patient_event_history_update, "SIMID", summarize, diff_SGTOT   = c(0,diff(SGTOT)))
  diff_SGACT   <- ddply(patient_event_history_update, "SIMID", summarize, diff_SGACT   = c(0,diff(SGACT)))
  diff_CWE_TOT <- ddply(patient_event_history_update, "SIMID", summarize, diff_CWE_TOT = c(0,diff(CWE_TOT)))
  
  patient_event_history_update$diff_ANLYEAR <- diff_ANLYEAR$diff_ANLYEAR
  patient_event_history_update$diff_FEVA    <- diff_FEVA$diff_FEVA
  patient_event_history_update$diff_SGTOT   <- diff_SGTOT$diff_SGTOT
  patient_event_history_update$diff_SGACT   <- diff_SGACT$diff_SGACT
  patient_event_history_update$diff_CWE_TOT <- diff_CWE_TOT$diff_CWE_TOT
  
  # Calculate model outcomes
  
  # Adjust the number of exacerbations (no exacerbation is allowed when pneumonia happened): this is to correct for potential double-counting at baseline
  patient_event_history_update[which(patient_event_history_update$PNEU_yn==1),"MODEXAC_yn"] <- 0
  patient_event_history_update[which(patient_event_history_update$PNEU_yn==1),"SEVEXAC_yn"] <- 0
  
  # Lung function: Mean FEV1 decline per year
  mean_annual_fev1_decline <- round(mean(patient_event_history_update$diff_FEVA)/mean(patient_event_history_update$diff_ANLYEAR),4)
  # Life expectancy
  mean_life_expectancy <- round(mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"]),4)
  # Moderate exacerbations
  modexac_year <-function(low_limit,upper_limit){sum(patient_event_history_update[which(patient_event_history_update$ANLYEAR>low_limit & patient_event_history_update$ANLYEAR<=upper_limit),"MODEXAC_yn"])}
  mean_mod_exa_rate <- round(modexac_year(0,1000)/(patient_size_input*mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"])),4)
  # Severe exacerbations
  sevexac_year <-function(low_limit,upper_limit){sum(patient_event_history_update[which(patient_event_history_update$ANLYEAR>low_limit & patient_event_history_update$ANLYEAR<=upper_limit),"SEVEXAC_yn"])}
  mean_sev_exa_rate <- round(sevexac_year(0,1000)/(patient_size_input*mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"])),4)
  # Exercise capacity
  mean_CWE_TOT_change <- round(mean(patient_event_history_update$diff_CWE_TOT)/mean(patient_event_history_update$diff_ANLYEAR),4)
  # SGRQ activity score 
  mean_SGACT_change <- round(mean(patient_event_history_update$diff_SGACT)/mean(patient_event_history_update$diff_ANLYEAR),4)
  # Cough/sputum
  # Total number of cough/sputum per patient during lifetime
  cum_coughsputum <- aggregate(patient_event_history_update$COUGHSPUTUM_yn, list(Patient = patient_event_history_update$SIMID), sum)
  # Rate per year
  mean_cough_rate <- round(mean(cum_coughsputum$x)/mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"]),4)
  # Create new time variable for symptoms
  diff_ANLYEAR_sym <- ddply(patient_event_history_update, "SIMID", summarize, diff_ANLYEAR_sym = c(diff(ANLYEAR),0))
  # Create new data frame for symptoms
  simulation_clinical_history_sym <- cbind(patient_event_history_update[,c("ANLYEAR", "COUGHSPUTUM_yn", "BREATHLESS_yn","dead")], diff_ANLYEAR_sym)
  simulation_dead                 <- simulation_clinical_history_sym[which(simulation_clinical_history_sym$dead == 1),c("SIMID", "ANLYEAR")]
  simulation_cough                <- simulation_clinical_history_sym[,c("SIMID", "diff_ANLYEAR_sym")]
  simulation_cough_time           <- simulation_clinical_history_sym$COUGHSPUTUM_yn * simulation_clinical_history_sym$diff_ANLYEAR_sym
  simulation_cough                <- cbind(simulation_cough, simulation_cough_time)
  mean_time_cough                 <- mean(aggregate(simulation_cough$simulation_cough_time, list(Patient = simulation_cough$SIMID), sum)$x/simulation_dead$ANLYEAR)
    
  # Shortness of breath
  # Total number of cough/sputum per patient during lifetime
  cum_breathless <- aggregate(patient_event_history_update$BREATHLESS_yn, list(Patient = patient_event_history_update$SIMID), sum)
  # Rate per year
  mean_breathless_rate <- round(mean(cum_breathless$x)/mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"]),4)
  simulation_breath      <- simulation_clinical_history_sym[,c("SIMID", "diff_ANLYEAR_sym")]
  simulation_breath_time <- simulation_clinical_history_sym$BREATHLESS_yn * simulation_clinical_history_sym$diff_ANLYEAR_sym
  simulation_breath      <- cbind(simulation_breath, simulation_breath_time)
  mean_time_breath       <- mean(aggregate(simulation_breath$simulation_breath_time, list(Patient = simulation_breath$SIMID), sum)$x/simulation_dead$ANLYEAR)
   
  # Adverse events (pneumonia)
  # Total number of pneumonias per patient during lifetime
  cum_pneu <- aggregate(patient_event_history_update$PNEU_yn, list(Patient = patient_event_history_update$SIMID), sum)
  # Rate per year
  mean_pneu_rate <- round(mean(cum_pneu$x)/mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"]),4)
  # Total number of pneumonias leading to hospitalisation per patient during lifetime
  cum_pneu_hosp <- aggregate(patient_event_history_update$pneu_hosp_yn, list(Patient = patient_event_history_update$SIMID), sum)
  # Rate per year
  mean_pneu_hosp_rate <- round(mean(cum_pneu_hosp$x)/mean(patient_event_history_update[which(patient_event_history_update$dead==1),"ANLYEAR"]),4)
  
  # SGRQ total score 
  mean_SGTOT_change <- round(mean(patient_event_history_update$diff_SGTOT)/mean(patient_event_history_update$diff_ANLYEAR),4)
  
  # QALYs are gender dependent: gender needed to calculate utilities. Read this from the baseline characteristics file and add column with gender to the simulated results 
  baseline_characteristics <- read.csv("Model - datasets/baseline_characteristics_predicted_data.csv",sep=",")
  patient_event_history_update <- merge(patient_event_history_update,baseline_characteristics[c("PTID","FEMALE")],by.x = "PTID",by.y = "PTID")
  # Order the dataset by "SIMID". This is very important because after merging the order is lost.
  patient_event_history_update <- patient_event_history_update[order(patient_event_history_update$SIMID,patient_event_history_update$ANLYEAR),]
  
  # Simulate the utilities and add them to the results 
  utilities <- round(0.9617-0.0013*(patient_event_history_update$SGTOT)-0.0001*((patient_event_history_update$SGTOT)^2) + ifelse(patient_event_history_update$FEMALE==0,0.0231,0),4)
  patient_event_history_update$utilities <- utilities
  # To calculate QALYS we multiply utilities by diff_ANLYEAR, then discount them
  QALYs <- patient_event_history_update$diff_ANLYEAR*patient_event_history_update$utilities
  discount_rate_QALYs <- ifelse(patient_event_history_update$ANLYEAR<=30,0.035,0.035) # UK = 0.035 always #Dynagito 0.04 0.02 #hardcoded
  QALYs_discounted    <- QALYs/(1+discount_rate_QALYs)^patient_event_history_update$ANLYEAR
  # Add QALYs to the results table
  patient_event_history_update$QALYs <- QALYs
  patient_event_history_update$QALYs_discounted <- QALYs_discounted
  # Finally, aggregate QALYS per patient and compute the average
  QALYs_patient <- aggregate(patient_event_history_update$QALYs,list(SIMID=patient_event_history_update$SIMID),sum)
  mean_qalys <- round(mean(QALYs_patient$x),4)
  QALYs_patient_discounted <- aggregate(patient_event_history_update$QALYs_discounted,list(SIMID=patient_event_history_update$SIMID),sum)
  mean_qalys_disc <- round(mean(QALYs_patient_discounted$x),4)
    
  # Costs are split into categories
  
  # Treatment costs: the model allows for a  distinction between health care and societal treatment costs, even though in the example below they are assumed to be the same.
  treatment_price_year_hc <- 1.12*1*365.25 # Hard-coded here
  treatment_price_year_societal <- treatment_price_year_hc
  # Add treatment costs to simulation results: discounted and undiscounted, health care and societal perspectives 
  patient_event_history_update$treatment_costs_hc <- patient_event_history_update$diff_ANLYEAR*treatment_price_year_hc
  discount_rate_costs <- ifelse(patient_event_history_update$ANLYEAR<=30,0.035,0.035) # Discount rates hard-coded
  patient_event_history_update$treatment_costs_hc_discounted <- patient_event_history_update$treatment_costs_hc/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  patient_event_history_update$treatment_costs_societal <- patient_event_history_update$diff_ANLYEAR*treatment_price_year_societal
  patient_event_history_update$treatment_costs_societal_discounted <- patient_event_history_update$treatment_costs_societal/(1+discount_rate_costs)^patient_event_history_update$ANLYEAR
  
  # Exacerbation costs 
  
  # Retirement age
  retirement_age <- 65 # hard-coded
  
  # Societal persp. & health care use: cost items
  exacerbation_costs_hc_use_row_names <- c("Primary care visits", "Secondary care visits", "Hospital days", "Ambulance rides", "ER visits", "Course antibiotics","Course oral steroids","Work days lost", "Distance primary care clinic", "Distance specialist clinic")
  # Read costs from Excel file depending on the selected options
  exacerbation_costs_hc_use <- read.csv("Model - datasets/Costs/exacerbation costs hc use.csv",sep=";",row.names = exacerbation_costs_hc_use_row_names)
  
  
  # Societal persp. & average: when health care use is not available the model allows inputting just average costs
  exacerbation_costs_average_row_names <- c("Societal", "Societal (retired)", "Health care", "Health care (retired)")
  exacerbation_costs_average <- read.csv("Model - datasets/Costs/exacerbation costs average.csv",sep=";",row.names = exacerbation_costs_average_row_names)
  
  #Assign the type of exacerbation costs read from the Excel
  mod_exa_costs_societal_average <- exacerbation_costs_average[1,1]   
  sev_exa_costs_societal_average <- exacerbation_costs_average[1,2]   
  mod_exa_costs_societal_average_retired <- exacerbation_costs_average[2,1]   
  sev_exa_costs_societal_average_retired <- exacerbation_costs_average[2,2]   
  
  # Health care persp. & average
  mod_exa_costs_hc_average <- exacerbation_costs_average[3,1]   
  sev_exa_costs_hc_average <- exacerbation_costs_average[3,2]   
  mod_exa_costs_hc_average_retired <- exacerbation_costs_average[4,1]   
  sev_exa_costs_hc_average_retired <- exacerbation_costs_average[4,2]   
    
  # The following functions calculate the cost per moderate and severe exacerbation depending on the cost types and perspective considered 
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
    
    ifelse(cost_type=="Health care use",sum(mod_exa_primary_care_cost,mod_exa_secondary_care_cost,mod_exa_hospital_days_cost,mod_exa_ambulance_ride_cost,mod_exa_ER_visit_cost,
                                            mod_exa_antibiotics_cost,mod_exa_steroids_cost,mod_exa_work_days_lost_cost,mod_exa_distance_pcc_cost,mod_exa_distance_sc_cost),mod_exa_average_cost)
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
    
    ifelse(cost_type=="Health care use", sum(sev_exa_primary_care_cost,sev_exa_secondary_care_cost,sev_exa_hospital_days_cost,sev_exa_ambulance_ride_cost,sev_exa_ER_visit_cost,
                                             sev_exa_antibiotics_cost,sev_exa_steroids_cost,sev_exa_work_days_lost_cost,sev_exa_distance_pcc_cost,sev_exa_distance_sc_cost), sev_exa_average_cost)
  }
    
  # Now calculate the exacerbation costs associated to exacerbations in the simulation. This has to be age-dependent when the perspective is societal (productivity lossses).
  # And apply it to the model results
  sev_exacerbation_costs_calc_sim <- function(index){
    if(patient_event_history_update[index,"ANLYEAR"]>0){sev_exa_cost(patient_event_history_update[index,"AGE_TIME"])*patient_event_history_update[index,"SEVEXAC_yn"]}
    else{0}
  }
  
  mod_exacerbation_costs_calc_sim <- function(index){
    if(patient_event_history_update[index,"ANLYEAR"]>0){mod_exa_cost(patient_event_history_update[index,"AGE_TIME"])*patient_event_history_update[index,"MODEXAC_yn"]}
    else{0}
  }
  
  # Choose first perspective and type of cost data: these should be input parameters!!!
  cost_type <- "Health care use" # Choose between "Health care use" or "Average" 
  perspective <- "Societal" # Choose between "Societal" or "Health care"
  
  # These are the simulated exacerbation costs
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
  gpvisits_stable_predictors <- c("FEMALE","BMI_CLASS_2","BMI_CLASS_3","SMOKER","SMPKY","OTHER_CVD"
                                  ,"DIABETES","DEPRESSION","HEART_FAILURE","ASTHMA","EMPHDIA","ICS","EOS_yn")
  
  gpvisits_predictors <- c(gpvisits_stable_predictors,"AGE_TIME","FEVPPA","SGACT","COUGHSPUTUM_yn","BREATHLESS_yn","SGTOT","MODEXAC_yn","SEVEXAC_yn")
  
  
  ### Merge the two datsets: I'm creating another file here because I don't want to save all the patient
  ### characteristics in the simulated results file
  patient_event_history_update_maintenance_costs <- merge(patient_event_history_update[,1:19],baseline_characteristics[c("PTID",gpvisits_stable_predictors)],by.x = "PTID",by.y = "PTID")
  
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
  
  #maintenance_costs
  
  ### Adjust then the number of gp and spec visits and multiply by the time interval where occurs
  num_gpvisits_observed <- 0.422
  num_gpvisits_country  <- maintenance_costs[1,1]
  patient_event_history_update$adjusted_num_gpvisits <- (patient_event_history_update$num_gpvisits*num_gpvisits_country/num_gpvisits_observed)*patient_event_history_update$diff_ANLYEAR
  
  
  num_specvisits_observed <- 0.545
  num_specvisits_country  <- maintenance_costs[2,1]
  patient_event_history_update$adjusted_num_specvisits <- (patient_event_history_update$num_specvisits*num_specvisits_country/num_specvisits_observed)*patient_event_history_update$diff_ANLYEAR
  
  ### Add the costs per visit
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
    
    if(patient_event_history_update[index,"PNEU_yn"]==1){
      if(patient_event_history_update[index,"pneu_hosp_yn"]==1){
        sev_exa_cost(patient_event_history_update[index,"AGE_TIME"])
      }else{
        mod_exa_cost(patient_event_history_update[index,"AGE_TIME"])
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
  
  
  ### Return model outcomes 
  
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
} #end COPD_model_simulation function

# To run the model simply call the model function with the appropriate inputs in the correct order. For example,
# the line below will run the model for 500 patients, deterministically, without treatment effects and with a 
# random seed equal to 177.

COPD_simulation_deterministic_results <- COPD_model_simulation(500, # Patient size
                                                               0, # 0 = deterministic run / 1 = probabilistic run
                                                               1, # Treatment effect: factor on time to exacerbation 
                                                               1, # Treatment effect: factor on probability of severe exacerbation
                                                               1, # Treatment effect: factor on fev1 change
                                                               1, # Treatment effect: factor on cwe score change
                                                               0, # Treatment effect: absolute change in SGRQ activity score
                                                               1, # Treatment effect: factor on probability of cough/sputum 
                                                               1, # Treatment effect: factor on probability of shortness of breath
                                                               0, # Treatment effect: absolute change in SGRQ total score
                                                               177) # random seed input


COPD_simulation_deterministic_results


# To run a probabilistic sensitivity analysis (PSA) we can use the function below, which is basically calling 
# the simulation function multiple times and getting the average results.

COPD_model_PSA <- function(psa_size_input,
                           patient_size_input,
                           exac_treatment_effect_tte_input,
                           exac_treatment_effect_sevexa_input,
                           fev1_treatment_effect_input,
                           cwe_treatment_effect_input,
                           sgact_treatment_effect_input,
                           coughsputum_treatment_effect_input,
                           breathless_treatment_effect_input,
                           sgtot_treatment_effect_input){
  
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
    
    current_psa <- COPD_model_simulation(patient_size_input,
                                         1, #note 1 here hardcoded because the main function will run in probabilistic mode
                                         exac_treatment_effect_tte_input,
                                         exac_treatment_effect_sevexa_input,
                                         fev1_treatment_effect_input,
                                         cwe_treatment_effect_input,
                                         sgact_treatment_effect_input,
                                         coughsputum_treatment_effect_input,
                                         breathless_treatment_effect_input,
                                         sgtot_treatment_effect_input,
                                         i) # note i is the random seed, which will differ per PSA iteration
    
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
} # end PSA function
