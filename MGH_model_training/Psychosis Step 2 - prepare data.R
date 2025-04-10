
## Define Global Functions

showCountandRate <- function(m,n,n_digits=1){ return(paste(m,' (',round(100 * m / n, digits = n_digits),'%)',sep='')) }

getSummaryStatsTbl <- function(roc.res){
  res = data.frame()
  for(spec in c(0.8,0.9,0.95,0.99)){
    res = rbind(res, round(coords(roc.res, x = spec, input = 'specificity', ret = c('spec','sen','ppv','npv', 'tp','fp','tn','fn'), transpose = F), digits = 2))
  }
  res$AUC = round(roc.res$auc, digits = 2)
  for(col in c('tp','fp','tn','fn')) res[,col] <- round(res[,col], digits = 0)
  return(res)
}


prepareModelDataUsingORCutoff <- function(df, odds_ratios, cutoff_lower = 0.01, cutoff_upper = 1.5, min_cases = 5){
  all_ids = unique(df$PATIENT_ID)
  selected_concepts = odds_ratios$concept_code[(odds_ratios$cases > min_cases) & (odds_ratios$OR < cutoff_lower | odds_ratios$OR >= cutoff_upper)]
  df$concept_type_code = paste(df$concept_type, df$concept_code,sep='_')
  df = df[df$concept_type_code %in% selected_concepts, ]
  df = rbind(df, data.frame(PATIENT_ID = all_ids, date = NA, age = NA, concept_type = 'stub', concept_code = 'stub', isCase = NA, time_till_index = NA, concept_type_code = 'stub_stub'))  
  # make sure we include all subjects (use a 'stub' concept code for that)
  tbl = as.data.frame.matrix(table(df$PATIENT_ID, df$concept_type_code))
  return(tbl)
}


addMissingLabel <- function(a, missing_label){
  a = as.character(a)
  a[is.na(a)] <- missing_label
  a = factor(a)
  return(a)
}

reduceLevels <- function(a, n_levels = 50){
  s = summary(a, n_levels)
  a = as.character(a)
  a[!(a %in% names(s))] <- 'Other'
  a = factor(a)
  return(a)
}

readAllDatasetSegments <- function(path_name){
  all_files = list.files(path_name, pattern = '*.csv', full.names = T)
  res = data.frame()
  count = 1
  for(fn in all_files){
    print(paste('Reading file',count,'/',length(all_files)))
    df = read.delim(fn, sep = '\t'); print(paste('number of rows =', nrow(df)))
    res = rbind(res, df)
    count = count + 1
  } 
  return(res)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# MAIN CODE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setwd("/data/js95/PNGU/RPDRml_v2021.1/")

## Prepare case definition dataset
# (this dataset includes 1.3M patients born between 1965-2005)
demographics <- readRDS('/data/js95/PNGU/RPDRml_v2021.1/sample2/RPDRml_demographics_sample2.RDs')


# load cases --------------------------------

SELECTED_OUTCOME = 'broad' # there are 3 possible case-definitions: broad, narroaw, and hybrid. We use 'broad' for now
cases = readRDS('./sample2/cases_broad_sample2.RDs')

# mark the index date as the date in which the patient first met the FEP case definition (by default use Table-1, if not available then Table-2)
cases$index_date = cases$table1_start
idxs = which(!is.na(cases$table2_start) & (cases$index_date > cases$table2_start))
if(length(idxs) > 0) cases$index_date[idxs] <- cases$table2_start[idxs]

# extract the age of the patient at the index event
cases$index_age = round(as.numeric(difftime(cases$index_date, cases$birth_date, units = 'days')) / 365, digits = 1)
saveRDS(cases, './sample2/cases.RDs') # cases = readRDS('cases.RDs')

# load controls --------------------------------

# There are 2 options for a control group: either patients meeting the controls case-def ('CONTROLS') or all patients who are not cases ('NON_CASES)
CONTROLS_TYPE = 'CONTROLS' #  NON_CASES

# 1) all non-cases
if (CONTROLS_TYPE == 'NON_CASES'){ 
  controls <- demographics_all[!(demographics_all$PATIENT_ID %in% cases$PATIENT_ID), ]
} else 
  # 2) all controls who meet controls-criteria
{ 
  controls = readRDS('./sample2/controls_sample2.RDs') 
}

# for controls also add a surrogate for an 'index event' which is the last documented information they have
# add 'time till index' for controls as well - this will be the time till last follow-up
# (this can be done only after the controls_last_timestamp is generated, later in code)

# get last visit from 'encounters' dataset
encounters = readRDS('./sample2/encounters2.RDs') # derived from diagnosis dataset
controls_last_timestamp = encounters[order(encounters$subject_num, encounters$sstart_date), ]
controls_last_timestamp = controls_last_timestamp[rev(!duplicated(rev(controls_last_timestamp$subject_num))), ]

controls$index_date = controls_last_timestamp$sstart_date[match(controls$subject_num, controls_last_timestamp$subject_num)]
controls$index_age = round(as.numeric(difftime(controls$index_date, controls$sbirth_date, units = 'days')) / 365, digits = 1)
saveRDS(controls, './sample2/controls.RDs')


# down-sample the controls for computational purposes
controls_ids_sample = sample(controls$subject_num, nrow(cases)*4, replace = F)

saveRDS(controls[controls$subject_num %in% controls_ids_sample, ], './sample2/controls_subsample.RDs')


## Load Files

setwd("/data/js95/PNGU/RPDRml_v2021.1/sample2/")
# load case definitions
cases = readRDS('cases.RDs')
cases_ids = cases$subject_num
controls = readRDS('controls.RDs')
controls_ids = controls$subject_num

controls.sample = readRDS('controls_subsample.RDs')
controls.sample_ids = controls.sample$PATIENT_ID

# subjects included in study = cases + controls
included_ids = c(cases_ids, controls_ids)
included_ids.sample = c(cases_ids, controls.sample_ids)


# load all data available for modeling
demographics_all <- readRDS('./RPDRml_demographics_sample2.RDs')
diagnosis = readRDS('./RPDRml__DX_sample2.RDs')
meds = readRDS('./RPDRml__MED_sample2.RDs') 
encounters = readRDS('./encounters2.RDs') # derived from diagnosis dataset
labs = readRDS('./RPDRml__LAB_sample2.RDs')  
procs = readRDS('./RPDRml__PROC_sample2.RDs')  

# # leave only information for included patients
# PARTNERS - was already done before
demographics <- demographics_all[demographics_all$subject_num %in% included_ids, ]
# diagnosis <- diagnosis[diagnosis$PATIENT_ID %in% included_ids, ]
# meds <- meds[meds$PATIENT_ID %in% included_ids, ]
# encounters <- encounters[encounters$PATIENT_ID %in% included_ids, ]
# labs <- labs[labs$PATIENT_ID %in% included_ids, ]
# procs <- procs[procs$PATIENT_ID %in% included_ids, ]
# insurance <- insurance[insurance$PATIENT_ID %in% included_ids, ]


# create a table including all index events
index_events = rbind(data.frame(subject_num = cases$subject_num, index_date = cases$index_date, isCase = 1),
                     data.frame(subject_num = controls$subject_num, index_date = controls$index_date, isCase = 0))
index_events$birth_date = demographics$sbirth_date[match(index_events$subject_num, demographics$subject_num)]
index_events$age_at_index = round(as.numeric(difftime(index_events$index_date, index_events$birth_date, units = 'days')) / 365, digits = 1)
saveRDS(index_events, 'index_events.RDs')

index_events = readRDS('index_events.RDs')



## Run preliminary descriptive statistics over each dataset
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Encounters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dim(encounters) # 43,484,336 | subject_num, encounter_num, sstart_date, site
paste(colnames(encounters), collapse = ', ')
summary(encounters$sstart_date)
# all_dates = encounters$sstart_date
# all_dates = all_dates[order(all_dates)]

# plot the dates of all available encounters, cases vs. controls, between given date-range (1975-2023)
# require(ggplot2)
# start_date = strptime('1975-01-01', '%Y-%m-%d')
# end_date = strptime('2023-01-01', '%Y-%m-%d')
# ggplot(data = encounters[encounters$START_DATE > start_date & encounters$START_DATE < end_date, ], aes(x = START_DATE, fill = 'Encounters')) + 
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', aes(y=..density..)) +
#   scale_fill_manual(values=c("#69b3a2", "#404080")) +
#   theme_minimal() +
#   labs(fill="")

# see how many encounters will we exclude with given cut-off dates
# start_date = strptime('1985-01-01', '%Y-%m-%d')
# end_date = strptime('2022-07-31', '%Y-%m-%d')
# sum(encounters$START_DATE < start_date)
# sum(encounters$START_DATE > end_date)
# encounters = encounters[encounters$START_DATE >= start_date & encounters$START_DATE <= end_date, ] #  17,874,169

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Diagnosis
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dim(diagnosis) # 159,960,185 | subject_num, encounter_num, concept_code, sstart_date, diagnosis_type, site, epic_fhir_ws, table1_narrow, table1_broad, table1_hybrid, table2, table3, table4
paste(colnames(diagnosis), collapse = ', ')

dim(meds) # 84,476,112 | subject_num, encounter_num, concept_code, sstart_date, site, epic_fhir_ws
paste(colnames(meds), collapse = ', ')
summary(meds$sstart_date) # 1997-12-20 till 2021-12-21

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Labs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dim(labs) # 316,432,428 | subject_num, encounter_num, concept_code, sstart_date, site, valueflag, valtype, nval, epic_fhir_ws
paste(colnames(labs), collapse = ', ')
summary(labs$sstart_date)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Procedures
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dim(procs) # 66655894 | PATIENT_ID, ENCOUNTER_NUM, START_DATE, END_DATE, CPT_CODE, CPT_NAME
paste(colnames(procs), collapse = ', ')
summary(procs$START_DATE) # 1997-12-20 till 2021-12-19

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##  Merge datasets for predictions - create model data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#
#
# Function: getModelRawData
# Input variables:
#   patient_id: id of patient
#   concept_date: date when concept was documented
#   concept_code: the code used to describe the concept (any text format)
#   demographics: the demographics dataset in order to extract birth date
#   index_events: the index dates for both cases (index event) and controls (last documented code)
#   cases_ids: the cases dataset to identify cases and also calculate time till index
#   concept_type: the prefix to add for the concept for final model table (e.g. 'dx', 'labs', etc.)
#
# Output: a dataset with the following format - PATIENT_ID | date | age | concept_type | concept_code | isCase | time_till_index
#
getModelRawData <- function(patient_id, concept_date, concept_code, demographics, index_events, cases_ids, concept_type){
  birth_dates = demographics$sbirth_date[match(patient_id, demographics$subject_num)]
  index_dates = index_events$index_date[match(patient_id, index_events$subject_num)]
  time_till_index = round(as.numeric(difftime(index_dates, concept_date, units = 'days')), digits = 0)
  age_at_concept = round(as.numeric(difftime(concept_date, birth_dates, units = 'days')) / 365, digits = 1)
  res = data.frame(subject_num = patient_id, date = concept_date, age = age_at_concept, concept_type = concept_type, concept_code = concept_code, isCase = (patient_id %in% cases_ids), time_till_index = time_till_index)  
  res = res[order(res$subject_num, res$date), ]
  res = res[!duplicated(paste(res$subject_num, res$concept_code)), ]
  return(res)
}

index_events = readRDS('index_events.RDs')

#
# Convert each dataset into the same format that can be taken into the final model
#

# filter to make sure our datasets only include patients that are included in the study (within index_events table)
diagnosis = diagnosis[diagnosis$subject_num %in% index_events$subject_num, ]
labs = labs[labs$subject_num %in% index_events$subject_num, ]
meds = meds[meds$subject_num %in% index_events$subject_num, ]
procs = procs[procs$subject_num %in% index_events$subject_num, ]
encounters = encounters[encounters$subject_num %in% index_events$subject_num, ]

# diagnosis
diagnosis.model <- getModelRawData(diagnosis$subject_num, diagnosis$sstart_date, diagnosis$concept_code, demographics, index_events, cases_ids, concept_type = 'dx')
saveRDS(diagnosis.model, './model_datasets/diagnosis.model.RDs')
# labs
labs.model <- getModelRawData(labs$subject_num, labs$sstart_date, labs$concept_code, demographics, index_events, cases_ids, concept_type = 'labs')
saveRDS(labs.model, './model_datasets/labs.model.RDs')
# meds
meds.model <- getModelRawData(meds$subject_num, meds$sstart_date, meds$concept_code, demographics, index_events, cases_ids, concept_type = 'meds')
saveRDS(meds.model, './model_datasets/meds.model.RDs')
# procedures
procs.model <- getModelRawData(procs$subject_num, procs$sstart_date, procs$concept_code, demographics, index_events, cases_ids, concept_type = 'procs')
saveRDS(procs.model, './model_datasets/procs.model.RDs')
# encounters
encounters.model <- getModelRawData(encounters$subject_num, encounters$sstart_date,encounters$site, demographics, index_events, cases_ids, concept_type = 'encounters')
saveRDS(encounters.model, './model_datasets/encounters.model.RDs')

# merge all datasets into one big dataset for later modeling                              
model.data <- rbind(diagnosis.model, labs.model, meds.model, procs.model, encounters.model) 
model.data = model.data[order(model.data$subject_num, model.data$date), ]

saveRDS(model.data, './model_datasets/model.data.RDs')



## Clear all datasets to free memory
# free memory after creating model.data data frame
rm(labs); rm(labs.model)
rm(procs); rm(procs.model)
rm(meds); rm(meds.model)
rm(diagnosis); rm(diagnosis.model)
gc()

## Get the cohorts for each landmark time-point (currently by age)

model.data = readRDS('./model_datasets/model.data.RDs')

#
# getPredictionCohort - truncates all information after prediction age and exclude patients 
# who met the case definition prior to the time of prediction. 
#
# input: 
#   model.data - the dataset prepared earlier, ready for modeling, including all data
#   PRED_AGE - the age at which prediction is made at (landmark), everything after is truncated. 
#              if PRED_AGE=NA then will predict using all data up to 30 days prior to index_date
#   train_size - ratio of training/testing (default 80%/20%) - this will only be a flag in dataset
# output: the model.data dataset just truncated after PRED_AGE and without excluded patients (also mark training/testing)
#
#
getPredictionCohort <- function(model.data, PRED_AGE, index_events, train_size = 0.8)
{
  # use only data available prior to the time when prediction made. 
  # If no age is given then use all data available 30-days and more prior to case definition
  if (!is.na(PRED_AGE)){
    pred.data = model.data[model.data$age < PRED_AGE, ]   
  } else {
    pred.data = model.data[is.na(model.data$time_till_index) | (!is.na(model.data$time_till_index) & model.data$time_till_index >= 30), ] 
  }
  
  all_possible_ids = unique(pred.data$subject_num)
  
  # Exclude patients who met the case definition prior to the prediction time (for cases) or lack follow-up data (controls)
  pred.data$age_at_index = index_events$age_at_index[match(pred.data$subject_num, index_events$subject_num)]
  if (!is.na(PRED_AGE)){
    excluded_ids1 = unique(pred.data$subject_num[is.na(pred.data$age_at_index) | (pred.data$age_at_index < PRED_AGE)])
  } else {
    excluded_ids1 = c()
  }
  
  
  # remove all information collected after the index event
  excluded_ids2 = unique(pred.data$subject_num[!is.na(pred.data$time_till_index) & (pred.data$time_till_index < 0)])
  
  # Exclude patients with insufficient follow-up (basically applies only for controls. If not defined than without follow-up till age 35)
  excluded_ids3 = c() # we remove this condition for now as we use survival analysis
  
  # Generate model-date using all data available up to the time of prediction
  pred.data <- pred.data[!(pred.data$subject_num %in% unique(c(excluded_ids1, excluded_ids2, excluded_ids3))), ]  
  
  training_ids = sample(unique(pred.data$subject_num), floor(train_size * sum(!duplicated(pred.data$subject_num))), replace = F)
  pred.data$isTraining = (pred.data$subject_num %in% training_ids)
  # Identify cases that meet the case definition within the selected time-window for prediction
  # (this was also removed as we work with survival analysis)
  
  return(pred.data)
}

# for each landmark time (ages 15,20,25,30) create the relevant model.data dataset
for (PRED_AGE in c(15, 20,25,30)){
  print(paste('PRED_AGE',PRED_AGE, sep = '='))
  pred.data = getPredictionCohort(model.data, PRED_AGE, index_events)
  saveRDS(pred.data, paste('./model_datasets/pred.data.age_',PRED_AGE,'.RDs',sep=''))
}




## Convert Raw data (long) into Model-Ready (flat/wide) Datasets


#
# genORTable - create a table of odds-ratio to outcome for each concept code
# input: 
#   concept_code: the independent risk factors (taken from model.data dataset)
#   case_def: binary flag (True/False) indicating if the code is taken from case (True) or control (False)
#   min_count: the minimum number of times each concept needs to appear in order to be included in final table
# output: a table with the odds ratio for each concept, including the following columns: 
#   concept_code, cases (cases count), controls (controls count), total count (cases+controls), OR (odds-ratio)
#   * the dataset is ordered from highest to lowest OR
#
genORTable <- function(concept_code, case_def, min_count = 20){
  
  total_cases = sum(case_def); total_controls = sum(!case_def)
  concept_code.cases = tapply(case_def, concept_code, sum)
  concept_code.controls = tapply(!case_def, concept_code, sum)
  
  res = data.frame(concept_code = names(concept_code.cases), cases = concept_code.cases, controls = concept_code.controls)
  res$total_count = res$cases + res$controls
  res <- res[res$total_count > min_count, ]
  
  res$OR <- round((res$cases / (total_cases - res$cases)) / (res$controls / (total_controls - res$controls)), digits = 1)
  res <- res[order(res$OR, res$cases, decreasing = T), ]
  return(res)
}

#
# prepareModelData - Convert Raw data (long) into Model-Ready (flat/wide) Datasets
# input:
#   df: the model.data data-frame (after preparing for prediction time/landmark)
#   odds_ratios: the list of odds-ratios of each feature as generated by the genORTable function
#   num_features_to_include: the number of features to include in model (selected by odds-ratio)
# output: a data.frame of all patients (rows) X all features (columns) - each with the number of occurrences of the code
#
# df = pred.data.unique; num_features_to_include = 1000

# df = pred.data; num_features_to_include = 3000; min_OR = 1.1
prepareModelData <- function(df, odds_ratios, demographics, num_features_to_include = NA, min_OR = 1.5){
  all_ids = unique(df$subject_num)
  if (!is.na(num_features_to_include)){
    odds_ratios = odds_ratios[order(abs(odds_ratios$OR), decreasing = T), ]
    tmp = odds_ratios[1:min(num_features_to_include, nrow(odds_ratios)), ]
    if(!is.na(min_OR)){
      selected_concepts = tmp$concept_code[abs(tmp$OR) > min_OR]    
    }
  } else {
    selected_concepts = odds_ratios$concept_code[odds_ratios$OR > min_OR]  
  }
  
  df$concept_type_code = paste(df$concept_type, df$concept_code,sep='_')
  df = df[df$concept_type_code %in% selected_concepts, c('subject_num', 'concept_type_code')]
  
  # to make sure that all patients appear in final table, create empty (stub) rows for patients with missing data
  missing_ids = all_ids[!(all_ids %in% df$subject_num)]
  df = rbind(df, data.frame(subject_num = missing_ids, concept_type_code = 'stub_stub'))
  
  tbl = as.data.frame.matrix(table(df$subject_num, df$concept_type_code))
  # remove the 'stub' column
  tbl = tbl[, -grep('stub_stub', colnames(tbl))]
  
  # add demographics to the model data 
  tbl$SEX = demographics$gender[match(rownames(tbl), as.character(demographics$subject_num))]
  tbl$MARITAL_STATUS = demographics$marital_status[match(rownames(tbl), as.character(demographics$subject_num))]
  tbl$race = factor(as.character(demographics$race[match(rownames(tbl), as.character(demographics$subject_num))]))
  tbl$ethnicity = factor(as.character(demographics$ethnicity[match(rownames(tbl), as.character(demographics$subject_num))]))
  tbl$veteran = factor(as.character(demographics$veteran[match(rownames(tbl), as.character(demographics$subject_num))]))
  tbl$public_payer = factor(as.character(demographics$public_payer[match(rownames(tbl), as.character(demographics$subject_num))]))
  
  return(tbl)
}


# for each landmark time (age at prediction) convert the raw (long) format to flat/wide format
for (PRED_AGE in c(15, 20,25,30)) {
  print(paste('PRED_AGE',PRED_AGE, sep = '='))
  pred.data = readRDS(paste('./model_datasets/pred.data.age_',PRED_AGE,'.RDs',sep=''))
  
  
  # Calculate the odds-ratio of each concept-code (count once for each subject), use only training-set
  training = pred.data[pred.data$isTraining == TRUE, ]
  
  training_ids = unique(training$subject_num)
  
  # training.unique <- training[!duplicated(paste(training$subject_num, training$concept_type, training$concept_code)), ] # --> this was done when creating the model.dataset
  odds_ratios = genORTable(concept_code = paste(training$concept_type, training$concept_code,sep='_'), 
                           case_def = training$isCase, 
                           min_count = 20)
  saveRDS(odds_ratios, paste('./model_datasets/odds_ratio.age_',PRED_AGE,'.RDs',sep=''))
  # odds_ratios = readRDS(paste('./model_datasets/odds_ratio.age_',PRED_AGE,'.RDs',sep=''))
  
  # limit the number of features in the model to include only features with relatively high OR (calculated using training set)
  # also convert the flat long concepts dataset into a single-row-per-patient format
  pred.data.flat = prepareModelData(pred.data,odds_ratios, demographics, num_features_to_include = 3000, min_OR = 1.1) #1.5
  
  # add case definition  
  pred.data.flat$isCase <- factor(ifelse(rownames(pred.data.flat) %in% as.character(cases$subject_num), 1, 0))
  
  # add trainin/testing flag
  pred.data.flat$isTraining <- ifelse(rownames(pred.data.flat) %in% as.character(training_ids), 1, 0)
  
  pred.data.flat$timeTillEvent = index_events$age_at_index[match(rownames(pred.data.flat), as.character(index_events$subject_num))] - PRED_AGE
  saveRDS(pred.data.flat, paste('./model_datasets/pred.data.flat.age_',PRED_AGE,'.RDs',sep=''))  
}








## code that was removed from prior analysis
```{r removed_code, echo = FALSE}
# getPredCohortStats <- function(pred.data, follow_up, PRED_AGE)
# {
#   total_included = sum(!duplicated(pred.data$PATIENT_ID))
#   
#   all_ages_stats = data.frame(follow_up, 
#                               PRED_AGE, 
#                               total_included , 
#                               cases = showCountandRate(sum(pred.data$isCase[!duplicated(pred.data$PATIENT_ID)]), total_included), 
#                               controls = showCountandRate(sum(!pred.data$isCase[!duplicated(pred.data$PATIENT_ID)]), total_included))
#   return(all_ages_stats)
# }
# 
# model.data = model.data[model.data$PATIENT_ID %in% included_ids, ]
# pred.data = getPredictionCohort(model.data, NA, NA)
# predCohortStats = getPredCohortStats(pred.data, NA, NA)
# 
# predCohortStats
# saveRDS(pred.data, 'pred.data.RDs')

# # Create prediction cohorts for each age & follow-up time window
# ALL_PRED_AGES = c(1,3,5)
# ALL_FOLLOW_UPS = c(15,20,25,30)
# 
# all_ages_stats = data.frame()
# for(follow_up in ALL_PRED_AGES){
#   for (PRED_AGE in ALL_FOLLOW_UPS){
#     pred.data = getPredictionCohort(model.data, follow_up, PRED_AGE)
#     all_ages_stats = rbind(all_ages_stats, getPredCohortStats(pred.data, follow_up, PRED_AGE))
#     saveRDS(pred.data, paste('pred.data.follow_up_',follow_up,'_age_',PRED_AGE,'.RDs',sep=''))
#   }
# }  
# write.csv(all_ages_stats, paste('all_ages_stats.csv') )



# 
# plot(x = 1:length(survival_cases_mean), y = survival_cases_mean, col = 'red', type  = 'b', cex = 0.5)
# points(x = 1:length(survival_controls_mean), y = survival_controls_mean, col = 'blue', type  = 'b', cex = 0.5)
# 
# plot.survival(rf.fit)
# 
# idx1 = which(y.pred$yvar$isCase == 1)
# idx2 = which(y.pred$yvar$isCase == 0)
# 
# idx1 = idx1[sample(1:length(idx1), 20, replace = F)]
# idx2 = idx2[sample(1:length(idx2), 20, replace = F)]
# 
# par(cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, mar = c(6.0,6,1,1), mgp = c(4, 1, 0))
# 
# plot(round(y.pred$time.interest,2),y.pred$survival[idx1[1],], type="l", xlab="Time (Year)", ylab="Survival", col=1, lty=1, lwd=2, xlim = c(0,10), ylim = c(0,1))
# 
# for(idx in idx1[-1]){
#   lines(round(y.pred$time.interest,2),y.pred$survival[idx,], type="l", xlab="Time (Year)", ylab="Survival", col=1, lty=1, lwd=2, xlim = c(0,10), ylim = c(0,1))
# }
# 
# for(idx in idx2){
#   lines(round(y.pred$time.interest,2), y.pred$survival[idx,], col=2, lty=2, lwd=2)  
# }
# 
# legend("topright", legend=c("Case","Control"), col=c(1:2), lty=1:2, cex=2, lwd=2)
# 
# rf.fit.subsample = subsample(rf.fit)
# print(rf.fit.subsample)
# 
# par(oma = c(0.5, 10, 0.5, 0.5))
# par(cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, mar = c(6.0,17,1,1), mgp = c(4, 1, 0))
# plot(rf.fit.subsample, xlab = "Variable Importance (x 100)", cex = 1.2)
```