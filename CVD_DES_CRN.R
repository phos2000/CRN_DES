## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(warning = FALSE, message = FALSE)
load.lib<-c("reshape2","tidyverse","ggrepel","dampack","kableExtra","flexsurv","readxl","haven","janitor","here","parallel", "Rcpp")
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)
library(simmer)
select <- dplyr::select

options("scipen"=1000, "digits"=4)
sourceCpp("cost_utility_cvd.cpp")

## ----inputs-------------------------------------------------------------------------------------------------------------------
params = list(
  ## basic/general from Stroke paper, 2023
  strategy = 0, # or 1
  vN = 1000,
  vHorizon = 100, # all people should terminate before horizon
  discount_rate = 0.03,
  wtp = 100000,
  CPI2017 = 245.12,
  CPI2019 = 255.657,
  CPI2020 = 258.811,
  
  ## events from paper 2017
  rrStatinsASCVD = 0.79,
  rrStatinsASCVD_low = 0.77,
  rrStatinsASCVD_upp = 0.81,
  rrStatinsASCVD_dist = "lognormal",
  
  hrAfterASCVD = 1.9,
  hrAfterASCVD_low = 1.60,
  hrAfterASCVD_upp = 2.40, 
  hrAfterASCVD_dist = "lognormal", # log-normal
  
  mortalityFirstYearASCVD_female_young = 0.09,
  mortalityFirstYearASCVD_female_young_dist = "beta",
  
  mortalityFirstYearASCVD_female_old = 0.3,
  mortalityFirstYearASCVD_female_old_dist = "beta",
  
  mortalityFirstYearASCVD_male_young = 0.14,
  mortalityFirstYearASCVD_male_young_dist = "beta",
  
  mortalityFirstYearASCVD_male_old = 0.25,
  mortalityFirstYearASCVD_male_old_dist = "beta",
  
  pMildStatinAdverse = 0.047,
  pMildStatinAdverse_low = 0.047*0.85, 
  pMildStatinAdverse_upp = 0.047*1.15,
  pMildStatinAdverse_dist = "beta",
  
  pMajorStatinAdverse = 0.006/100,
  pMajorStatinAdverse_low = 0.006/100*0.85,
  pMajorStatinAdverse_upp = 0.006/100*1.15,
  pMajorStatinAdverse_dist = "beta",
  
  pMajorStatinAdverseBeingFatal = 0.09/100,
  pMajorStatinAdverseBeingFatal_low = 0.09/100*0.85,
  pMajorStatinAdverseBeingFatal_upp = 0.09/100*1.15,
  pMajorStatinAdverseBeingFatal_dist = "beta",

  ## utility
  uHealthy = 1,
  uAfterASCVD = 0.773,
  uAfterASCVD_low = 0.773*0.85,
  uAfterASCVD_upp = 0.773*1.15,
  uAfterASCVD_dist = "beta",
  
  uHealthyStatin = 0.996,
  uHealthyStatin_low = 0.991,
  uHealthyStatin_upp = 1.000,
  uHealthyStatin_dist = "beta",
  
  uPenaltyMildStatinAdverse = -0.0055,
  uPenaltyMajorStatinAdverse = -0.0383,
  
  uDeath = 0,
  
  ## cost(2017)
  cAnnualFU_afterASCVD = 3917,
  cAnnualFU_afterASCVD_low = 2611,
  cAnnualFU_afterASCVD_upp = 6528,
  cAnnualFU_afterASCVD_dist = "gamma",
  
  cNonFatalASCVD = 49348,
  cNonFatalASCVD_low = 49348*0.85,
  cNonFatalASCVD_upp = 49348*1.15,
  cNonFatalASCVD_dist = "gamma",
  
  cFatalASCVD = 16760,
  cFatalASCVD_low = 16760*0.85,
  cFatalASCVD_upp = 16760*1.15,
  cFatalASCVD_dist = "gamma",
  
  cAnnualStatin = 84,
  
  cMildStatinAdverse = 215,
  cMildStatinAdverse_low = 196,
  cMildStatinAdverse_upp = 326,
  cMildStatinAdverse_dist = "gamma",
  
  cMajorStatinAdverse = 8486,
  cMajorStatinAdverse_low = 6920, 
  cMajorStatinAdverse_upp = 10314,
  cMajorStatinAdverse_dist = "gamma",
  
  annual_discount_rate = 0.03,
  
  ### source：Agency for Healthcare Research and Quality. Mean expenditure per person by age groups, United States, 1996 to 2020. Medical Expenditure Panel Survey.(2020)
  cHealthcareNonCVD_18_44 = 3901,
  cHealthcareNonCVD_45_64 = 8693,
  cHealthcareNonCVD_65 = 12441

)
## inflation
params$rInflation2017 = params$CPI2020 / params$CPI2017
params$rInflation2019 = params$CPI2020 / params$CPI2019

## leave some spaces for PSA


## ----func---------------------------------------------------------------------------------------------------------------------
# simply change from prob to rate than to event
ProbToRate = function(prob, t){
  -log(1-prob)/t
}

RateToProb = function(rate, t){
  1 - exp(-rate*t)
}
# # what this t should be -- 10 years?
# # time unit: year
# rate = ProbToRate(pcr()/100)
# days = rexp(1,rate)

params$cont_discount_rate <- -log(1- params$annual_discount_rate) # Yearly Time Scale
discounted <- function(undiscounted, start_year, end_year, rate = params$cont_discount_rate)
{
  undiscounted/ (exp(rate*(end_year-start_year)))
}

discount <- function(value, A, ar=params$annual_discount_rate) value / (1+ar)^A


## ----generator----------------------------------------------------------------------------------------------------------------
rCRN = 1:8
names(rCRN) = c("nhanesNo", "rDeathWithoutCVD", "rGetCVD", "rDeathOfASCVD", "rDeathAfterASCVD", "rStatinsMildAdverse", "rStatinsMajorAdverse", "rDeathOfStatinsAdverse")
# 
# params_CRN = params
# 
# params_CRN$nRN = params_CRN$vN * length(rCRN)
# params_CRN$randomNums = runif(params_CRN$nRN)


## ----source-------------------------------------------------------------------------------------------------------------------
# source("nonCVD_life_table.R")


## ----life-table-gompertz------------------------------------------------------------------------------------------------------
# # male
# 
# mltNonCVD_long = data.frame(Age = rep(NA, ceiling(sum(mlt2020_nonCVD$dx))), Death = rep(1, ceiling(sum(mlt2020_nonCVD$dx))))
# 
# start = 1
# for (age in 0:(nrow(mlt2020_nonCVD)-1)) {
#   n_row = round(mlt2020_nonCVD[age+1,"dx"])
#   mltNonCVD_long[seq(start,(n_row + start - 1)),"Age"] = age
#   start = start + n_row
# }
# 
# mltNonCVD_long = mltNonCVD_long %>% drop_na()
# 
# parameters_m = data.frame(
#   Age = mlt2020_nonCVD$Age,
#   shape=rep(NA,nrow(mlt2020_nonCVD)),
#   rate=rep(NA,nrow(mlt2020_nonCVD)))
# 
# 
# for (age in 0:nrow(mlt2020_nonCVD)-1) {
#   surv.data <- with(mltNonCVD_long[mltNonCVD_long$Age >= age,], Surv(Age, Death, origin=age))
#   surv.model <- flexsurvreg(surv.data ~ 1, dist="gompertz")
# 
#   parameters_m$shape[age+1] <- surv.model$coefficients[1]
#   parameters_m$rate[age+1] <- exp(surv.model$coefficients[2])
# }
# 
# # if start since 45-yr
# # pgompertz(cycle, parameters_m$shape[46], parameters_m$rate[46])
# 
# # female
# 
# fltNonCVD_long = data.frame(Age = rep(NA, ceiling(sum(flt2020_nonCVD$dx))), Death = rep(1, ceiling(sum(flt2020_nonCVD$dx))))
# 
# start = 1
# for (age in 0:(nrow(flt2020_nonCVD)-1)) {
#   n_row = round(flt2020_nonCVD[age+1,"dx"])
#   fltNonCVD_long[seq(start,(n_row + start - 1)),"Age"] = age
#   start = start + n_row
# }
# 
# fltNonCVD_long = fltNonCVD_long %>% drop_na()
# 
# parameters_f = data.frame(
#   Age = flt2020_nonCVD$Age,
#   shape=rep(NA,nrow(flt2020_nonCVD)),
#   rate=rep(NA,nrow(flt2020_nonCVD)))
# 
# 
# for (age in 0:nrow(flt2020_nonCVD)-1) {
#   surv.data <- with(fltNonCVD_long[fltNonCVD_long$Age >= age,], Surv(Age, Death, origin=age))
#   surv.model <- flexsurvreg(surv.data ~ 1, dist="gompertz")
# 
#   parameters_f$shape[age+1] <- surv.model$coefficients[1]
#   parameters_f$rate[age+1] <- exp(surv.model$coefficients[2])
# }

ageAtDeath <- function(currentAge, gender, inputs, patientID)
{
  # patientID = get_attribute(env, "patientID")
  currentAge_f <- floor(currentAge) #negative time to secular death issue
  
  # Male 1 FEMALE 2
  if(gender == 2){
    shape <- parameters_f$shape[currentAge_f+1]
    rate  <- parameters_f$rate[currentAge_f+1]
  } else {
    shape <- parameters_m$shape[currentAge_f+1]
    rate  <- parameters_m$rate[currentAge_f+1]
  }
  
  min(currentAge + ifelse(is.null(inputs$randomNums), rgompertz(1, shape, rate), qgompertz(inputs$randomNums[(patientID-1)*length(rCRN) + rCRN["rDeathWithoutCVD"]], shape, rate)), 100)
}



## ----gompertz-hr--------------------------------------------------------------------------------------------------------------
# # male
# mlt2020_hr = data.frame(Age = 0:(nrow(mlt2020_nonCVD)-1), qx = sapply(sapply(mlt2020_nonCVD$qx, ProbToRate, t = 1) * inputs$hrAfterASCVD, RateToProb, t = 1), lx = rep(0,nrow(mlt2020_nonCVD)), dx = rep(NA,nrow(mlt2020_nonCVD)))
# 
# mlt2020_hr[1,"lx"] = 100000
# 
# for (i in 1:nrow(mlt2020_hr)) {
#   if (mlt2020_hr[i, "lx"] > 0) {
#     mlt2020_hr[i, "dx"] = round(mlt2020_hr[i, "qx"] * mlt2020_hr[i, "lx"])
#     mlt2020_hr[i+1, "lx"] = mlt2020_hr[i, "lx"] - mlt2020_hr[i, "dx"]
#   }
# }
# mlt2020_hr  = mlt2020_hr %>% filter(!is.na(dx))
# 
# mltHr_long = data.frame(Age = rep(NA, 100000), Death = rep(1, 100000))
# 
# start = 1
# for (age in 0:(nrow(mlt2020_hr)-1)) {
#   n_row = mlt2020_hr[age+1,"dx"]
#   mltHr_long[seq(start,(n_row + start - 1)),"Age"] = age
#   start = start + n_row
# }
# 
# parameters_m_hr = data.frame(
#   Age = mlt2020_hr$Age,
#   shape=rep(NA,nrow(mlt2020_hr)),
#   rate=rep(NA,nrow(mlt2020_hr)))
# 
# for (age in 0:nrow(mlt2020_hr)-1) {
#   surv.data <- with(mltHr_long[mltHr_long$Age >= age,], Surv(Age, Death, origin=age))
#   surv.model <- flexsurvreg(surv.data ~ 1, dist="gompertz")
# 
#   parameters_m_hr$shape[age+1] <- surv.model$coefficients[1]
#   parameters_m_hr$rate[age+1] <- exp(surv.model$coefficients[2])
# }
# 
# # female
# 
# flt2020_hr = data.frame(Age = 0:(nrow(flt2020_nonCVD)-1), qx = sapply(sapply(flt2020_nonCVD$qx, ProbToRate, t = 1) * inputs$hrAfterASCVD, RateToProb, t = 1), lx = rep(0,nrow(flt2020_nonCVD)), dx = rep(NA,nrow(flt2020_nonCVD)))
# 
# flt2020_hr[1,"lx"] = 100000
# 
# for (i in 1:nrow(flt2020_hr)) {
#   if (flt2020_hr[i, "lx"] > 0) {
#     flt2020_hr[i, "dx"] = round(flt2020_hr[i, "qx"] * flt2020_hr[i, "lx"])
#     flt2020_hr[i+1, "lx"] = flt2020_hr[i, "lx"] - flt2020_hr[i, "dx"]
#   }
# }
# flt2020_hr  = flt2020_hr %>% filter(!is.na(dx))
# 
# fltHr_long = data.frame(Age = rep(NA, 100000), Death = rep(1, 100000))
# 
# start = 1
# for (age in 0:(nrow(flt2020_hr)-1)) {
#   n_row = flt2020_hr[age+1,"dx"]
#   fltHr_long[seq(start,(n_row + start - 1)),"Age"] = age
#   start = start + n_row
# }
# 
# parameters_f_hr = data.frame(
#   Age = flt2020_hr$Age,
#   shape=rep(NA,nrow(flt2020_hr)),
#   rate=rep(NA,nrow(flt2020_hr)))
# 
# for (age in 0:nrow(flt2020_hr)-1) {
#   surv.data <- with(fltHr_long[fltHr_long$Age >= age,], Surv(Age, Death, origin=age))
#   surv.model <- flexsurvreg(surv.data ~ 1, dist="gompertz")
# 
#   parameters_f_hr$shape[age+1] <- surv.model$coefficients[1]
#   parameters_f_hr$rate[age+1] <- exp(surv.model$coefficients[2])
# }
# 
ageAtDeath_afterASCVD = function(currentAge, gender, inputs, patientID)
{
  # patientID = get_attribute(env,"patientID")
  currentAge_f <- floor(currentAge) #negative time to secular death issue

  # Male 1 FEMALE 2
  if(gender == 2)
  {
    shape <- parameters_f_hr$shape[currentAge_f+1]
    rate  <- parameters_f_hr$rate[currentAge_f+1]
    age = min(currentAge + ifelse(is.null(inputs$randomNums), rgompertz(1, shape, rate), qgompertz(inputs$randomNums[(patientID-1)*length(rCRN) + rCRN["rDeathAfterASCVD"]], shape, rate)), max(parameters_f_hr$Age))
  }
  else
  {
    shape <- parameters_m_hr$shape[currentAge_f+1]
    rate  <- parameters_m_hr$rate[currentAge_f+1]
    age = min(currentAge + ifelse(is.null(inputs$randomNums), rgompertz(1, shape, rate), qgompertz(inputs$randomNums[(patientID-1)*length(rCRN) + rCRN["rDeathAfterASCVD"]], shape, rate)), max(parameters_m_hr$Age))
  }

  age
}
# 
# save(parameters_f, parameters_m, parameters_f_hr, parameters_m_hr, file = "life table/MarkovCVD_gompertz.RData")

load("life table/MarkovCVD_gompertz.RData")


## ----Pooled-Cohort------------------------------------------------------------------------------------------------------------
source('pcr.R')
# # the probability are percents!!


## ----nhanes-clearance---------------------------------------------------------------------------------------------------------
# ## SEQN: Respondent sequence number
# # Examination: Body Measures
# ## BMXBMI - Body Mass Index (kg/m**2)
# ## 12.3 to 86.2	Range of Values
# ## .	Missing
# nhanes_data_1 = data.frame(read_xpt("./NHANES/BMX_J.XPT", col_select = c("SEQN","BMXBMI")))
# # Lab: Glycohemoglobin
# ## LBXGH: Glycohemoglobin (%)
# nhanes_data_2 = data.frame(read_xpt("./NHANES/GHB_J.XPT", col_select = c("SEQN","LBXGH")))
# # Examination: Blood Pressure
# ## BPXSY1 - Systolic: Blood pres (1st rdg) mm Hg
# ## BPXDI1 - Diastolic: Blood pres (1st rdg) mm Hg
# ## 0 to 136	Range of Values
# ## .	Missing
# nhanes_data_3 = data.frame(read_xpt("./NHANES/BPX_J.XPT", col_select = c("SEQN","BPXDI1","BPXSY1")))
# # Demographic Variables and Sample Weights
# ## RIAGENDR - Gender
# ## 1	Male
# ## 2	Female	
# ## .	Missing
# ## RIDAGEYR - Age in years at screening
# ## 0 to 79	Range of Values
# ## 80	80 years of age and over
# ## .	Missing
# 
# ## RIDRETH1 - Race/Hispanic origin
# ## 1	Mexican American	1367	1367	
# ## 2	Other Hispanic	820	2187	
# ## 3	Non-Hispanic White	3150	5337	
# ## 4	Non-Hispanic Black	2115	7452	
# ## 5	Other Race - Including Multi-Racial	1802	9254	
# ## .	Missing
# 
# ## WTMEC2YR - Full sample 2 year MEC exam weight
# ## 2566.1838545 to 419762.83649	Range of Values
# ## 0	Not MEC Examined
# ## .	Missing
# nhanes_data_4 = data.frame(read_xpt("./NHANES/DEMO_J.XPT", col_select = c("SEQN","RIDAGEYR","RIAGENDR","RIDRETH1","WTMEC2YR"))) #WTMEC2YR is the sample weights variable
# # Questionnaire: Diabetes
# ## DIQ010 - Doctor told you have diabetes
# ## 1	Yes
# ## 2	No
# ## 3	Borderline
# ## 7	Refused
# ## 9	Don't know
# ## .	Missing
# nhanes_data_5 = data.frame(read_xpt("./NHANES/DIQ_J.XPT", col_select = c("SEQN","DIQ010")))
# # Lab: Cholesterol - High - Density Lipoprotein
# ## LBDHDD: Direct HDL-Cholesterol (mg/dL)
# nhanes_data_6 = data.frame(read_xpt("./NHANES/HDL_J.XPT", col_select = c("SEQN","LBDHDD")))
# # Questionnaire: medical conditions
# ## MCQ160d: Ever told you had angina/angina pectoris
# ## MCQ160e: Ever told you had heart attack
# ## MCQ160f: Ever told you had a stroke
# ## 1	Yes
# ## 2	No
# ## 7	Refused	
# ## 9	Don't know
# ## .	Missing
# nhanes_data_7 = data.frame(read_xpt("./NHANES/MCQ_J.XPT", col_select = c("SEQN","MCQ160D","MCQ160F","MCQ160E")))
# # Questionnaire: Smoking - Cigarette Use
# ## SMQ040: Do you now smoke cigarettes?
# ## 1	Every day
# ## 2	Some days
# ## 3	Not at all	
# ## 7	Refused
# ## 9	Don't know
# ## .	Missing
# nhanes_data_8 = data.frame(read_xpt("./NHANES/SMQ_J.XPT", col_select = c("SEQN","SMQ040")))
# # Lab: Cholesterol - Total
# ## LBXTC: Total Cholesterol (mg/dL)
# nhanes_data_9 = data.frame(read_xpt("./NHANES/TCHOL_J.XPT", col_select = c("SEQN","LBXTC")))
# nhanes_data_10 = data.frame(read_xpt("./NHANES/BPQ_J.XPT", col_select = c("SEQN", "BPQ050A")))
# # Lab: Cholesterol - Low-Density Lipoproteins (LDL)
# ## LBDLDL - LDL-Cholesterol, Friedewald (mg/dL)
# nhanes_data_11 = data.frame(read_xpt("./NHANES/TRIGLY_J.XPT", col_select = c("SEQN", "LBDLDL")))
# 
# nhanes_data_raw <- list(nhanes_data_1,nhanes_data_2,nhanes_data_3,nhanes_data_4,
#                         nhanes_data_5,nhanes_data_6,nhanes_data_7,nhanes_data_8,
#                         nhanes_data_9,nhanes_data_10,nhanes_data_11)
# nhanes_data_raw <- data.frame(reduce(nhanes_data_raw, full_join, by='SEQN'))
# nhanes_pop = nhanes_data_raw %>%
#   # whether needs to check missing data
#   filter(WTMEC2YR != 0 & RIDAGEYR >= 40 & RIDAGEYR <= 79 & LBXTC >= 30 & LBXTC <= 500 & LBDHDD >= 5 & LBDHDD <= 200 & BPXSY1 >= 60 & BPXSY1 <= 200) %>%
#   drop_na() %>%
#   mutate(RIDRETH1 = ifelse(RIDRETH1 == 3, 1, ifelse(RIDRETH1 == 4, 2, ifelse(RIDRETH1 %in% c(1,2), 3, 4))))
# 
# nhanes_pop$cum_weight = cumsum(nhanes_pop$WTMEC2YR / sum(nhanes_pop$WTMEC2YR))
# 
race = 1:4
names(race) = c("White", "African American", "Hispanic", "Other")
race_convert = function(num) {
  return(names(race[num]))
}

bp_treatment = 1:2
names(bp_treatment) = c("Yes","No")

smoker = 1:3
names(smoker) = c("Every day","Some days","Not at all")

diabetic = 1:2
names(diabetic) = c("Yes","No")
# 
# save(nhanes_pop, nhanes_data_raw, file = "./NHANES/nhanes_full.RData")

load("./NHANES/nhanes_full.RData")


## ----get-pop------------------------------------------------------------------------------------------------------------------
# for CRN, add an identifier
initialize_patient <- function(traj, inputs)
{
        traj %>%
                seize("time_in_model")       %>%
                set_attribute("patientID", function()
                  # patientID is 0-index
                  as.numeric(str_remove(get_name(env), "patient"))) %>%
                # NHANES version
                set_attribute("nhanesNo",         function() ifelse(
                  is.null(inputs$randomNums),
                  sample(1:nrow(nhanes_pop), 1, prob=1/nhanes_pop$WTMEC2YR),
                  which(nhanes_pop$cum_weight>=inputs$randomNums[(get_attribute(env,"patientID"))*length(rCRN) + rCRN["nhanesNo"]])[1])) %>%
                set_attribute("Gender",     function() nhanes_pop$RIAGENDR[get_attribute(env, "nhanesNo")]) %>%
                set_attribute("Age",        function() nhanes_pop$RIDAGEYR[get_attribute(env, "nhanesNo")]) %>% 
                set_attribute("AgeInitial", function() get_attribute(env, "Age")) %>%
                set_attribute("Race",       function() nhanes_pop$RIDRETH1[get_attribute(env, "nhanesNo")]) %>%
                set_attribute("TotChol",    function() nhanes_pop$LBXTC[get_attribute(env, "nhanesNo")]) %>%
                set_attribute("HdlChol",    function() nhanes_pop$LBDHDD[get_attribute(env, "nhanesNo")]) %>%
                set_attribute("LdlChol",   function()
nhanes_pop$LBDLDL[get_attribute(env, "nhanesNo")]) %>%
                set_attribute("SystolicBp", function() nhanes_pop$BPXSY1[get_attribute(env, "nhanesNo")]) %>%
                set_attribute("BpTreatment",function() nhanes_pop$BPQ050A[get_attribute(env, "nhanesNo")]) %>%
                set_attribute("Smoker",      function() nhanes_pop$SMQ040[get_attribute(env, "nhanesNo")]) %>%
                set_attribute("Diabetic",    function() nhanes_pop$DIQ010[get_attribute(env, "nhanesNo")]) %>%
                set_attribute("PrsZ",       function() 0) %>%
                set_attribute("PCErisk",    function()
  pcr(gender = ifelse(get_attribute(env, "Gender") == 1, 'M', 'F'),
             get_attribute(env,"AgeInitial"),
             race_convert(get_attribute(env, "Race")),
             get_attribute(env, "TotChol"),
             get_attribute(env, "HdlChol"),
             get_attribute(env, "SystolicBp"),
             (get_attribute(env, "BpTreatment") == 1),
             (get_attribute(env, "Smoker") %in% c(1,2)),
             (get_attribute(env, "Diabetic") == 1),
             get_attribute(env, "PrsZ")
             )[1]/100) %>%
                set_attribute("Strategy", function() inputs$strategy) %>%
                set_attribute("aUseStatins", function()
                  ifelse(get_attribute(env, "Strategy") == 1 &
                           (get_attribute(env, "LdlChol") >= 190 |
                            get_attribute(env, "PCErisk") >= 0.075 |
                            get_attribute(env, "Diabetic") == 1), 1, 0))
}


## ----events-------------------------------------------------------------------------------------------------------------------
counters = c("get_CVD", "death_without_CVD", "death_of_ASCVD", "death_after_ASCVD", "statins_mild_adverse_event", "statins_major_adverse_event", "death_of_statins_adverse_event", "time_in_model")

### Note:
### As only simulate all the things once, lots of places I use the initial age instead of the updated age. 

terminate <- function(traj, inputs)
{
  traj %>%
  branch(function() 1,
         continue=FALSE,
         trajectory()               %>%
           release("time_in_model") %>%
           timeout(0)
        )
}

# # Statins events
# years_till_use_statins <- function(inputs){
#   ifelse(get_attribute(env, "Strategy") == 1 &
#                            (get_attribute(env, "LdlChol") >= 190 |
#                             get_attribute(env, "PCErisk") >= 0.075 |
#                             get_attribute(env, "Diabetic") == 1), 0, Inf)
# }
# 
# event_use_statins <- function(traj, inputs){
#   traj %>%
#     mark("use_statins")
# }

years_till_statins_mild_adverse <- function(inputs){
  if (get_attribute(env, "aUseStatins") == 1 & is.na(get_attribute(env, "aStatinsMildAdverse"))) {
      ifelse(is.null(inputs$randomNums),
      ifelse(runif(1)<inputs$pMildStatinAdverse,0, Inf),
      ifelse(inputs$randomNums[(get_attribute(env,"patientID"))*length(rCRN) + rCRN["rStatinsMildAdverse"]]<inputs$pMildStatinAdverse, 0, Inf)
      )
  } else Inf
}

event_statins_mild_adverse <- function(traj, inputs){
  traj %>% 
    mark("statins_mild_adverse_event")
}

years_till_statins_major_adverse <- function(inputs){
  if (get_attribute(env, "aUseStatins") == 1 & is.na(get_attribute(env, "aStatinsMajorAdverse"))) {
      ifelse(is.null(inputs$randomNums),
      ifelse(runif(1)<inputs$pMajorStatinAdverse,0, Inf),
      ifelse(inputs$randomNums[(get_attribute(env,"patientID"))*length(rCRN) + rCRN["rStatinsMajorAdverse"]]<inputs$pMajorStatinAdverse, 0, Inf)
      )
  } else Inf
}

event_statins_major_adverse <- function(traj, inputs){
  traj %>% 
    mark("statins_major_adverse_event")
}

years_till_death_of_statins_adverse <- function(inputs){
  time_to_major_adverse = get_attribute(env, "aStatinsMajorAdverse")
  if (!is.infinite(time_to_major_adverse)){
    ifelse(is.null(inputs$randomNums),
    ifelse(runif(1)<inputs$pMajorStatinAdverseBeingFatal, time_to_major_adverse +0, Inf),
    ifelse(inputs$randomNums[(get_attribute(env,"patientID"))*length(rCRN) + rCRN["rDeathOfStatinsAdverse"]] < inputs$pMajorStatinAdverseBeingFatal, time_to_major_adverse+0, Inf))
  } else Inf
}

event_death_of_statins_adverse <- function(traj, inputs){
   traj %>% 
    set_attribute("AgeDeathStatinsAdverse", function() get_attribute(env,"aDeathOfStatinsAdverse") + get_attribute(env,"AgeInitial")) %>%
    branch(
    function() 1,
    continue=c(TRUE),
    trajectory("Death associated with major adverse events") %>%
      mark("death_of_statins_adverse_event") %>%
      terminate(inputs)
  )
}

## CVD events

years_till_death_without_CVD <- function(inputs)
{
  age       <- get_attribute(env, 'Age')
  gender    <- get_attribute(env, 'Gender')
  patientID <- get_attribute(env, 'patientID')
  
  deathAge = ageAtDeath(age, gender, inputs, patientID)
  
  ## the if-else here make sure that if get CVD firstly, will never die of 1.0 BG mortality
  
  return(ifelse(get_attribute(env, "aGetCVD") < deathAge - age, Inf, deathAge - age))
}

event_death_without_CVD <- function(traj, inputs)
{
  cat("event_Death() ", now(env), "\n")
  
  traj %>% 
    set_attribute("AgeDeathNonCVD", function() get_attribute(env,"aDeathWithoutCVD") + get_attribute(env,"AgeInitial")) %>%
    branch(
    function() 1,
    continue=c(TRUE),
    trajectory("Death without CVD") %>%
      mark("death_without_CVD") %>%
      terminate(inputs)
  )
}

years_till_get_CVD <- function(inputs)
{
  gender = ifelse(get_attribute(env, "Gender") == 1, "M", "F")
  race_pcr = race_convert(get_attribute(env, "Race")) 
  
  prob = get_attribute(env, "PCErisk")
  rr = ifelse(inputs$strategy == 1, inputs$rrStatinsASCVD, 1)
  
  # 10-year probability to 1-year rate
  years = ifelse(is.null(inputs$randomNums), rexp(1,ProbToRate(prob, 10)*rr), qexp(inputs$randomNums[(get_attribute(env,"patientID"))*length(rCRN) + rCRN["rGetCVD"]], ProbToRate(prob, 10)*rr))
  
  return(ifelse(is.na(get_attribute(env, "aGetCVD")), years, Inf))
}

event_get_CVD = function(traj, inputs) {
  cat("event_get_CVD() ", now(env), "\n")

  traj %>%
    set_attribute("AgeCVD", function() get_attribute(env,"AgeInitial") + get_attribute(env, "aGetCVD")) %>%
    branch(
    function() 1,
    continue=c(TRUE),
    trajectory("Get CVD") %>%
      mark("get_CVD"))
}

years_till_death_of_ASCVD = function(inputs) {
  
  ageCVD    <- get_attribute(env, 'aGetCVD') + get_attribute(env, "AgeInitial")
  gender    <- get_attribute(env, 'Gender')
  
  ## getCVD maybe sth really big.
  
  age_max = ifelse(gender == 1, max(parameters_m_hr$Age), max(parameters_f_hr$Age))
  
  if (ageCVD >= age_max - 1) {
    # not able to use ageAtDeath(ageCVD+1, gender)
    years = get_attribute(env, 'aGetCVD')
  } else {
    ## ageCVD under the age limitation -- able to use gomepertz
    # whether ASCVD died in first year happened
    happen = 0
    if (gender == 1) {
      ## men
      if (ageCVD < 65) {
        happen = ifelse(is.null(inputs$randomNums),
                        sample(0:1, 1, prob = c(1-inputs$mortalityFirstYearASCVD_male_young, inputs$mortalityFirstYearASCVD_male_young)),
                        ifelse(inputs$randomNums[(get_attribute(env,"patientID"))*length(rCRN) + rCRN["rDeathOfASCVD"]] < inputs$mortalityFirstYearASCVD_male_young, 1, 0))
        } else {
          happen = ifelse(is.null(inputs$randomNums),
                        sample(0:1, 1, prob = c(1-inputs$mortalityFirstYearASCVD_male_old, inputs$mortalityFirstYearASCVD_male_old)),
                        ifelse(inputs$randomNums[(get_attribute(env,"patientID"))*length(rCRN) + rCRN["rDeathOfASCVD"]] < inputs$mortalityFirstYearASCVD_male_old, 1, 0))
          }
      } else {
        ## women
        if (ageCVD < 65) {
          happen = ifelse(is.null(inputs$randomNums),
                        sample(0:1, 1, prob = c(1-inputs$mortalityFirstYearASCVD_female_young, inputs$mortalityFirstYearASCVD_female_young)),
                        ifelse(inputs$randomNums[(get_attribute(env,"patientID"))*length(rCRN) + rCRN["rDeathOfASCVD"]] < inputs$mortalityFirstYearASCVD_female_young, 1, 0))
          } else {
            happen = ifelse(is.null(inputs$randomNums),
                        sample(0:1, 1, prob = c(1-inputs$mortalityFirstYearASCVD_female_old, inputs$mortalityFirstYearASCVD_female_old)),
                        ifelse(inputs$randomNums[(get_attribute(env,"patientID"))*length(rCRN) + rCRN["rDeathOfASCVD"]] < inputs$mortalityFirstYearASCVD_female_old, 1, 0))
          }
      }
    years = ifelse(happen == 1, get_attribute(env, "aGetCVD") + 0.5, Inf)
  }
  
  return(years)
}

event_death_of_ASCVD = function(traj, inputs) {
  traj %>%
    set_attribute("AgeDeathCVD", function() get_attribute(env, "AgeInitial") + now(env)) %>%
    branch(
    function() 1,
    continue=c(TRUE),
    trajectory("Death of first-year ASCVD") %>%
      mark("death_of_ASCVD") %>%
      terminate(inputs)
    )
  }

years_till_death_after_ASCVD <- function(inputs) {
  # whether die of ASCVD(in <= 1 yr)
  if (is.infinite(get_attribute(env, 'aDeathOfASCVD'))) {
  
    ageCVD    <- get_attribute(env, 'aGetCVD') + get_attribute(env, "AgeInitial")
    gender    <- get_attribute(env, 'Gender')
    patientID <- get_attribute(env, 'patientID')

    deathAge = ageAtDeath_afterASCVD(ageCVD + 1, gender, inputs, patientID)
  
    # GetCVD + ageInitial maybe a impossible number
    years = deathAge - get_attribute(env, "AgeInitial")
  
  } else {
    years = Inf
  }
  
  return(years)
}

event_death_after_ASCVD = function(traj, inputs) {
  cat("event_death_with_CVD", now(env), "\n")
  
  traj %>% 
    set_attribute("AgeDeathCVD", function() get_attribute(env, "AgeInitial")+now(env)) %>%
    branch(
    function() 1,
    continue=c(TRUE),
    trajectory("Death after ASCVD event") %>%
      mark("death_after_ASCVD") %>%
      terminate(inputs)
    )
}

event_registry <- list(
        list(name          = "Get CVD",
             attr          = "aGetCVD",
             time_to_event = years_till_get_CVD,
             func          = event_get_CVD,
             reactive      = FALSE),
        list(name          = "Death without CVD",
             attr          = "aDeathWithoutCVD",
             time_to_event = years_till_death_without_CVD,
             func          = event_death_without_CVD,
             reactive      = FALSE),
        list(name          = "Death of ASCVD",
             attr          = "aDeathOfASCVD",
             time_to_event = years_till_death_of_ASCVD,
             func          = event_death_of_ASCVD,
             reactive      = FALSE),
        list(name          = "Death after ASCVD",
             attr          = "aDeathAfterASCVD",
             time_to_event = years_till_death_after_ASCVD,
             func          = event_death_after_ASCVD,
             reactive      = FALSE),
        list(name          = "Halt at time horizon",
             attr          = "aTerminate",
             time_to_event = function(inputs) inputs$vHorizon,
             func          = terminate,
             reactive      = FALSE),
        # list(name          = "Use Statins",
        #      attr          = "aUseStatins",
        #      time_to_event = years_till_use_statins,
        #      func          = event_use_statins,
        #      reactive      = FALSE),
        list(name          = "Mild Adverse for Statins",
             attr          = "aStatinsMildAdverse",
             time_to_event = years_till_statins_mild_adverse,
             func          = event_statins_mild_adverse,
             reactive      = FALSE),
        list(name          = "Major Adverse for Statins",
             attr          = "aStatinsMajorAdverse",
             time_to_event = years_till_statins_major_adverse,
             func          = event_statins_major_adverse,
             reactive      = FALSE),
        list(name          = "Death of Major Adverse for Statins",
             attr          = "aDeathOfStatinsAdverse",
             time_to_event = years_till_death_of_statins_adverse,
             func          = event_death_of_statins_adverse,
             reactive      = FALSE)
)
## aAgeCVD is not immediately(only get after this happen), but aGetCVD happens at 0.

## in this way there'll be 2 aGetCVD
## because if getCVD is inside the horizon, there'll be repeats. for the first version, they just died/terminated. 


## ----run----------------------------------------------------------------------------------------------------------------------
source("DES.R")

env = simmer("CVD")

exec.simulation <- function(inputs)
{
        #set.seed(11451428)
        env  <<- simmer("CVD")
        traj <- simulation(env, inputs)
        env %>% create_counters(counters)
        
        env %>%
                add_generator("patient", traj, at(rep(0, inputs$vN)), mon=2) %>%
                run(inputs$vHorizon+1e-6) %>% # Simulate just past horizon
                wrap()
        
        get_mon_arrivals(env, per_resource = T)
}

results <- NULL
attributes <- NULL

for (strategy in 0:1){
  params$strategy = strategy

  run <- exec.simulation(params)
  run$strategy <- strategy

  at <- arrange(get_mon_attributes(env),name,key,time) #obtain attributes data
  at$strategy <- strategy

  if(is.null(results)) { results <- run } else  {results <- rbind(results, run)}
  if(is.null(attributes)) { attributes <- at } else  {attributes <- rbind(attributes, at)}
  rm(run)
  rm(at)
}



## ----cost-utility-func--------------------------------------------------------------------------------------------------------
# set the cHealthcareNonCVD to be linear
x <- c(31.5, 55, 80)
y <- c(3901, 8693, 12441)
params$model_bgcost <- lm(y ~ x)

x <- c(50, 60, 70, 80, 90)
y <- c((10.1+4.2)/2/100, (21.4+8.9)/2/100, (34.6+20.0)/2/100, (59.2+40.2)/2/100, (74.4+65.2)/2/100)
params$model_cvdpct <- lm(y ~ x)

cHealthy = function(inputs, age) {
  x <- c(31.5, 55, 80)
  y <- c(3901, 8693, 12441)

  model_bgcost <- lm(y ~ x)
  intercept_bgcost = coef(model_bgcost)[1]
  slope_bgcost = coef(model_bgcost)[2]

  x <- c(50, 60, 70, 80, 90)
  y <- c((10.1+4.2)/2/100, (21.4+8.9)/2/100, (34.6+20.0)/2/100, (59.2+40.2)/2/100, (74.4+65.2)/2/100)

  model_cvdpct <- lm(y ~ x)
  intercept_cvdpct = coef(model_cvdpct)[1]
  slope_cvdpct = coef(model_cvdpct)[2]
  
  cvdpct = max(intercept_cvdpct + slope_cvdpct*age,0)

  (intercept_bgcost + slope_bgcost*age - cvdpct*(inputs$cAnnualFU_afterASCVD + inputs$cNonFatalASCVD/10)) / (1 - cvdpct)
}

## ----discounted---------------------------------------------------------------------------------------------------------------
basic_CVD_costs_discounted = function(inputs, death_after_ASCVD, death_of_ASCVD, death_without_CVD, get_CVD, time_in_model, initial_age) {
  cost = 0
  ##cHealthcareNonCVD
  cost = cost + discounted(sum(sapply(initial_age:floor(min(ifelse(is.na(get_CVD),Inf,get_CVD),time_in_model)+initial_age), function(x) cHealthy(inputs, x))), start_year = 0, end_year = min(ifelse(is.na(get_CVD),Inf,get_CVD),time_in_model))
  
  # FatalASCVD
  if(!is.na(get_CVD) & !is.na(death_of_ASCVD)) {
    cost = cost + discount(inputs$rInflation2017 * inputs$cFatalASCVD, A = death_of_ASCVD)
  }
  
  # Non-Fatal ASCVD
  if(!is.na(get_CVD) & is.na(death_of_ASCVD)) {
    cost = cost + discount(inputs$rInflation2017 * inputs$cNonFatalASCVD, A= get_CVD) + discounted(inputs$rInflation2017 * inputs$cAnnualFU_afterASCVD * (time_in_model - get_CVD - 0.5), start_year = get_CVD + 0.5, end_year = time_in_model)
  }
  return(cost)
}

adjustment_statins_costs_discounted = function(inputs, statins_use, statins_mild_adverse_event, statins_major_adverse_event, get_CVD, time_in_model){
  if (statins_use == 1){
    discounted(inputs$cAnnualStatin * min(ifelse(is.na(get_CVD),Inf,get_CVD),time_in_model) * inputs$rInflation2017, start_year = 0, end_year = min(ifelse(is.na(get_CVD),Inf,get_CVD),time_in_model)) +
    ifelse(!is.na(statins_major_adverse_event), discount(inputs$cMajorStatinAdverse, A=0), ifelse(!is.na(statins_mild_adverse_event), discount(inputs$cMildStatinAdverse,A=0), 0))
  } else 0
}

adjusted_statins_utilities_discounted = function(inputs, statins_use, statins_mild_adverse_event, statins_major_adverse_event, get_CVD, time_in_model){
  if (statins_use == 1) {
    if (!is.na(statins_major_adverse_event)){
      uStatins = inputs$uHealthyStatin - inputs$uPenaltyMajorStatinAdverse
    } else if (!is.na(statins_mild_adverse_event)) {
      uStatins = inputs$uHealthyStatin - inputs$uPenaltyMildStatinAdverse
    } else {
      uStatins = inputs$uHealthyStatin
    }
    
    if (is.na(get_CVD)) { discounted(uStatins * time_in_model, start_year = 0, end_year = time_in_model)}
    else { discounted(uStatins * get_CVD, start_year = 0, end_year = get_CVD) + discounted(inputs$uAfterASCVD * (time_in_model - get_CVD), start_year = get_CVD, end_year = time_in_model)}
  
    } else {
    if (is.na(get_CVD)) {
      return(discounted(inputs$uHealthy * time_in_model, start_year = 0, end_year = time_in_model))
    } else {
      return(discounted(inputs$uHealthy * get_CVD, start_year = 0, end_year = get_CVD) + discounted(inputs$uAfterASCVD * (time_in_model - get_CVD), start_year = get_CVD, end_year = time_in_model))
      }
  }
}

## ----draws--------------------------------------------------------------------------------------------------------------------
source("random_draws_for_CVD_PSAs.R")

RIPS_CVD_run = function(inputs){
  results <- NULL
  attributes <- NULL
  #outputs <- list()

  for (strategy in 0:1){
    inputs$strategy = strategy
  
    run <- exec.simulation(inputs)
    run$strategy <- strategy
  
    at <- arrange(get_mon_attributes(env),name,key,time) #obtain attributes data
    at$strategy <- strategy
  
    if(is.null(results)) { results <- run } else  {results <- rbind(results, run)}
    if(is.null(attributes)) { attributes <- at } else  {attributes <- rbind(attributes, at)}
    rm(run)
    rm(at)
  }
  
  repeatSets = attributes %>% filter((key == "aGetCVD" | key == "aStatinsMildAdverse" | key == "aStatinsMajorAdverse") & is.infinite(value))

  casted = setdiff(attributes, repeatSets) %>% select(-replication, -time) %>% spread(key, value)
  
  traces = results %>% select(-start_time, -activity_time, -replication) %>% spread(resource, end_time) %>% 
  left_join(casted %>% select(name, strategy, statins_use = aUseStatins, initial_age = AgeInitial), by = c("name", "strategy"))

  if (is.null(traces$statins_major_adverse_event)) {traces$statins_major_adverse_event = NA}
  
  # patientWcea = traces %>% 
  #   rowwise() %>%
  #   mutate(
  #     cost_disc = basic_CVD_costs_discounted(inputs, death_after_ASCVD, death_of_ASCVD, death_without_CVD, get_CVD, time_in_model, initial_age) + adjustment_statins_costs_discounted(inputs, statins_use, statins_mild_adverse_event, statins_major_adverse_event, get_CVD, time_in_model),
  #     util_disc = adjusted_statins_utilities_discounted(inputs, statins_use, statins_mild_adverse_event, statins_major_adverse_event, get_CVD, time_in_model))

  cost_util = calc_cost_util(death_of_ASCVD_vec = traces$death_of_ASCVD,
                             get_CVD_vec = traces$get_CVD, 
                             time_in_model_vec = traces$time_in_model, 
                             initial_age_vec = traces$initial_age,
                             statins_use_vec = traces$statins_use, 
                             statins_mild_adverse_event_vec = traces$statins_mild_adverse_event, 
                             statins_major_adverse_event_vec = traces$statins_major_adverse_event, 
                             n_pop = nrow(traces),
                             cr = inputs$cont_discount_rate,
                             ar = inputs$annual_discount_rate,
                             bgcost_coef_vec = coef(inputs$model_bgcost), 
                             cvdpct_coef_vec = coef(inputs$model_cvdpct),
                             cAnnualFU_afterASCVD = inputs$cAnnualFU_afterASCVD,
                             cNonFatalASCVD = inputs$cNonFatalASCVD,
                             cFatalASCVD = inputs$cFatalASCVD,
                             cAnnualStatin = inputs$cAnnualStatin,
                             cMildStatinAdverse = inputs$cMildStatinAdverse,
                             cMajorStatinAdverse = inputs$cMajorStatinAdverse,
                             rInflation2017 = inputs$rInflation2017,
                             uHealthy = inputs$uHealthy,
                             uAfterASCVD = inputs$uAfterASCVD,
                             uHealthyStatin = inputs$uHealthyStatin,
                             uPenaltyMildStatinAdverse = inputs$uPenaltyMildStatinAdverse,
                             uPenaltyMajorStatinAdverse = inputs$uPenaltyMajorStatinAdverse)
  
  patientWcea = cbind(traces, cost_util)
  
  strategyWcea = patientWcea %>% 
    group_by(strategy) %>%
      summarize(
        time_in_model = mean(time_in_model),
        cost = mean(cost_disc), 
        QALY = mean(util_disc)) %>%
    mutate(strategy = ifelse(strategy == 1, "Statins(2013 ACC/AHA)", "Status Quo"))

  df_cea = calculate_icers(cost = strategyWcea$cost,
                         effect = strategyWcea$QALY,
                         strategies = strategyWcea$strategy) %>%
    left_join(strategyWcea %>% select(Strategy = strategy, LE = time_in_model))

  return(df_cea)
} 