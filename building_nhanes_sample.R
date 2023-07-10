## Load packages
  library(tidyverse)
  library(janitor)
  library(dplyr)
  library(readxl)
  library(haven)
  library(mice)

# Read SAS -----------------------------------------------------------------------

  # ----nhanes-clearance---------------------------------------------------------------------------------------------------------
  ## SEQN: Respondent sequence number
  # Examination: Body Measures
  ## BMXBMI - Body Mass Index (kg/m**2)
  ## 12.3 to 86.2	Range of Values
  ## .	Missing
  nhanes_data_1 = data.frame(read_xpt("./NHANES/BMX_J.XPT", col_select = c("SEQN","BMXBMI")))
  # Lab: Glycohemoglobin
  ## LBXGH: Glycohemoglobin (%)
  nhanes_data_2 = data.frame(read_xpt("./NHANES/GHB_J.XPT", col_select = c("SEQN","LBXGH")))
  # Examination: Blood Pressure
  ## BPXSY1 - Systolic: Blood pres (1st rdg) mm Hg
  ## BPXDI1 - Diastolic: Blood pres (1st rdg) mm Hg
  ## 0 to 136	Range of Values
  ## .	Missing
  nhanes_data_3 = data.frame(read_xpt("./NHANES/BPX_J.XPT", col_select = c("SEQN","BPXDI1","BPXSY1")))
  # Demographic Variables and Sample Weights
  ## RIAGENDR - Gender
  ## 1	Male
  ## 2	Female
  ## .	Missing
  ## RIDAGEYR - Age in years at screening
  ## 0 to 79	Range of Values
  ## 80	80 years of age and over
  ## .	Missing

  ## RIDRETH1 - Race/Hispanic origin
  ## 1	Mexican American	1367	1367
  ## 2	Other Hispanic	820	2187
  ## 3	Non-Hispanic White	3150	5337
  ## 4	Non-Hispanic Black	2115	7452
  ## 5	Other Race - Including Multi-Racial	1802	9254
  ## .	Missing

  ## WTMEC2YR - Full sample 2 year MEC exam weight
  ## 2566.1838545 to 419762.83649	Range of Values
  ## 0	Not MEC Examined
  ## .	Missing
  nhanes_data_4 = data.frame(read_xpt("./NHANES/DEMO_J.XPT", col_select = c("SEQN","RIDAGEYR","RIAGENDR","RIDRETH1","WTMEC2YR"))) #WTMEC2YR is the sample weights variable
  # Questionnaire: Diabetes
  ## DIQ010 - Doctor told you have diabetes
  ## 1	Yes
  ## 2	No
  ## 3	Borderline
  ## 7	Refused
  ## 9	Don't know
  ## .	Missing
  nhanes_data_5 = data.frame(read_xpt("./NHANES/DIQ_J.XPT", col_select = c("SEQN","DIQ010")))
  # Lab: Cholesterol - High - Density Lipoprotein
  ## LBDHDD: Direct HDL-Cholesterol (mg/dL)
  nhanes_data_6 = data.frame(read_xpt("./NHANES/HDL_J.XPT", col_select = c("SEQN","LBDHDD")))
  # Questionnaire: medical conditions
  ## MCQ160d: Ever told you had angina/angina pectoris
  ## MCQ160e: Ever told you had heart attack
  ## MCQ160f: Ever told you had a stroke
  ## 1	Yes
  ## 2	No
  ## 7	Refused
  ## 9	Don't know
  ## .	Missing
  nhanes_data_7 = data.frame(read_xpt("./NHANES/MCQ_J.XPT", col_select = c("SEQN","MCQ160D","MCQ160F","MCQ160E")))
  # Questionnaire: Smoking - Cigarette Use
  ## SMQ040: Do you now smoke cigarettes?
  ## 1	Every day
  ## 2	Some days
  ## 3	Not at all
  ## 7	Refused
  ## 9	Don't know
  ## .	Missing
  nhanes_data_8 = data.frame(read_xpt("./NHANES/SMQ_J.XPT", col_select = c("SEQN","SMQ040")))
  # Lab: Cholesterol - Total
  ## LBXTC: Total Cholesterol (mg/dL)
  nhanes_data_9 = data.frame(read_xpt("./NHANES/TCHOL_J.XPT", col_select = c("SEQN","LBXTC")))
  nhanes_data_10 = data.frame(read_xpt("./NHANES/BPQ_J.XPT", col_select = c("SEQN", "BPQ050A")))
  # Lab: Cholesterol - Low-Density Lipoproteins (LDL)
  ## LBDLDL - LDL-Cholesterol, Friedewald (mg/dL)
  nhanes_data_11 = data.frame(read_xpt("./NHANES/TRIGLY_J.XPT", col_select = c("SEQN", "LBDLDL")))

  nhanes_data_raw <- list(nhanes_data_1,nhanes_data_2,nhanes_data_3,nhanes_data_4,
                          nhanes_data_5,nhanes_data_6,nhanes_data_7,nhanes_data_8,
                          nhanes_data_9,nhanes_data_10,nhanes_data_11)
  nhanes_data_raw <- data.frame(reduce(nhanes_data_raw, full_join, by='SEQN'))

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


# Filter and Adjust -------------------------------------------------------
  
## Filter the data to only respondents over 40 years old and < 80 yrs old
## NHANES age: 0 to 79, 80 and over
  nhanes_pop = nhanes_data_raw %>%
    # whether needs to check missing data
    filter(WTMEC2YR != 0 & 
             RIDAGEYR >= 40 & 
             RIDAGEYR <= 79 & 
             LBXTC >= 30 & 
             LBXTC <= 500 & 
             LBDHDD >= 5 & 
             LBDHDD <= 200 & 
             BPXSY1 >= 60 & 
             BPXSY1 <= 200) %>%
    mutate(RIDRETH1 = ifelse(RIDRETH1 == 3, 1, ifelse(RIDRETH1 == 4, 2, ifelse(RIDRETH1 %in% c(1,2), 3, 4))))
  
  # check sample size now
  nrow(nhanes_pop)

## Calculating missingness in NHANES data
  missingness = (colMeans(is.na(nhanes_pop)))*100
  
## Imputing data
  
  init = mice(nhanes_pop, maxit=0) 
  predM = init$predictorMatrix
  predM[,"SEQN"]=0
  nhanes_imputed_test <- mice(nhanes_pop, predictorMatrix = predM, m = 5, seed = 11)
  nhanes_pop <- complete(nhanes_imputed_test)
  
## Sample with replacement
  nhanes_pop$cum_weight = cumsum(nhanes_pop$WTMEC2YR / sum(nhanes_pop$WTMEC2YR))
  
  # save(nhanes_pop, nhanes_data_raw, file = "./NHANES/nhanes_full.RData")
  
#   nhanes_sample = sample_n(nhanes_pop[!is.na(nhanes_pop$WTMEC2YR),], size = 1000000, replace = TRUE, weight = na.omit(nhanes_pop$WTMEC2YR))
#   
#   
#   nhanes_sample = nhanes_sample %>%
#     rowwise() %>%
#     mutate(category = sample(1:4, size = 1, 
#                              replace = TRUE, prob = p_cat_among_stenosis)) %>%
#     mutate(category = ifelse(stenosis_occurrence == 0, 0, category)) %>%
#     mutate(pce = pcr(gender = ifelse(male == 1, 'M', 'F'),
#                      age = RIDAGEYR,
#                      race = ifelse(race == 1, 'African American', 'Other'), 
#                      tot_chol = LBXTC,
#                      hdl_chol = LBDHDD,
#                      systolic_bp = BPXSY1,
#                      bp_treatment = (BPQ050A == 1),
#                      smoker = (SMQ040 == 1),
#                      diabetic = (DIQ010 == 1),
#                      prs_z = 0)[[1]])
#   
#   
# # Save baseline pop
#   save(nhanes_sample, file = "basepop.RData")
  
  
  