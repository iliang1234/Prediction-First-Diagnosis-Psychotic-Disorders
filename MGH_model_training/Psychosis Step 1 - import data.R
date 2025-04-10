RELOAD_DATA <- TRUE
APPLY_CASE_DEFINITION <- TRUE


# Fields that should be modified
# subject_num
# sbirth_date

# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Define Global Functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -

getBirthYear <- function(dates_str){
  if (class(dates_str) == 'character' | class(dates_str) == 'factor'){
    dates_str <- as.character(dates_str)
    res = as.numeric(gsub('.+/.+/','',dates_str))  
  } else {
    dates_str <- as.character(dates_str)
    res = as.numeric(gsub('-.+-.+','',dates_str))
  }
  return(res)
}

getFactorSummaryString <- function(a){
  s = summary(a); s = s[order(s, decreasing = T)]
  res = paste(paste(paste(names(s), round(100*s / length(a), digits = 1), sep = '='), collapse = '%, '), '%',sep='')
  return(res)
}

readAllDatasetSegments <- function(path_name, included_ids = NA){
  all_files = list.files(path_name, pattern = '*.csv', full.names = T)
  res = data.frame()
  count = 1
  for(fn in all_files){
    print(paste('Reading file',count,'/',length(all_files)))
    df = read.delim(fn, sep = '\t'); print(paste('number of rows =', nrow(df), '(total=',nrow(res),')'))
    if(length(included_ids) > 1) df = df[df$PATIENT_ID %in% included_ids, ]
    res = rbind(res, df)
    count = count + 1
  } 
  return(res)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Import txt files into RDs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setwd("/data/js95/PNGU/RPDRml_v2021.1")

time_format = '%m/%d/%Y' # 11/17/1943
# other possible formats: '%Y-%m-%d %H:%M:%OS'

# demographics
demographics = read.delim('./RPDRml__demographics.txt', sep = '\t', stringsAsFactors = F)
demographics$BIRTH_DATE <- as.POSIXct(strptime(demographics$sbirth_date, format = time_format))
demographics$SEX_CD <- factor(demographics$gender)
demographics$MARITAL_STATUS_CD <- factor(demographics$marital_status)

demographics$RACE_CD <- factor(demographics$race)
demographics$VETERAN_CD <- factor(demographics$veteran)
demographics$ETHNICITY_CD <- factor(demographics$ethnicity)

saveRDS(demographics, 'demographics.RDs')

# ethnicity -- at Partners this is part of the demographics dataset
# ethnicity = readAllDatasetSegments('./Ethnicity/')
# ethnicity$ETHNICITY <- factor(ethnicity$ETHNICITY)
# ethinicities_per_subject = tapply(ethnicity$ETHNICITY, ethnicity$PATIENT_ID, function(x) paste(unique(x[order(x)]), collapse = ', '))
# demographics$ethnicity = ethinicities_per_subject[as.character(demographics$PATIENT_ID)]
# saveRDS(demographics, 'demographics.RDs')

# insurance -- at Partners this is part of the demographics dataset 
# insurance = readAllDatasetSegments('./Insurance/')
# saveRDS(insurance, 'insurance.RDs')

# diagnosis
diagnosis = read.delim('./RPDRml__DX.txt', sep = '\t', stringsAsFactors = F)
diagnosis$START_DATE <- as.POSIXct(strptime(diagnosis$START_DATE, format = '%Y-%m-%d %H:%M:%OS'))
saveRDS(diagnosis, 'diagnosis.RDs')

# The remaining datasets are loaded at the end of the script


# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load datasets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
list.files('./sample2/')
demographics = readRDS('./sample2/RPDRml_demographics_sample2.RDs') # 1.3M
diagnosis = readRDS('./sample2/RPDRml__DX_sample2.RDs') # 160M

## Run Inclusion Criteria

# 1. Sub-Sample to meet inclusion criteria
# Include subjects born between 1965-2005 (ages 15-35 between years 2000-2020)
all_subjects = demographics
MIN_BIRTH_DATE = strptime("01/01/1965", format = '%m/%d/%Y')
MAX_BIRTH_DATE = strptime("01/01/2005", format = '%m/%d/%Y')
included_ids = unique(all_subjects$subject_num[all_subjects$sbirth_date >= MIN_BIRTH_DATE & all_subjects$sbirth_date < MAX_BIRTH_DATE])
saveRDS(included_ids, './sample2/included_ids.RDs') 

# diagnosis - get diagnosis data for the selected subset of patients 
# PARTNERS: no need for this. Already ran this part. 
# diagnosis.included = diagnosis[diagnosis$PATIENT_ID %in% included_ids, ] 
# saveRDS(diagnosis.included, 'diagnosis.included.RDs')


# demographics - get demographics data for the selected subset of patients
# PARTNERS: no need for this. Already ran this part. 
# demographics.included = demographics[demographics$PATIENT_ID %in% included_ids, ]
# saveRDS(demographics.included, 'demographics.included.RDs')

# extract encounters from diagnosis data
# PARTNERS: no need for this. Already ran this part. 
# encounters = diagnosis.included[order(diagnosis.included$PATIENT_ID, diagnosis.included$START_DATE), ]
# encounters = encounters[!duplicated(paste(encounters$PATIENT_ID, encounters$START_DATE)), ]
# encounters = encounters[, c('PATIENT_ID', 'ENCOUNTER_NUM',  'START_DATE')]
# saveRDS(encounters, 'encounters.included.RDs')


# - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load Case Definition
# - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Psychosis Primary Inclusion Code Groups (broad-based, narrow-based, or hybrid definition)
#   - Psychosis, unspecified or other
#   - Major depression with psychotic features 
#   - Bipolar disorder with psychotic features 
#   - Schizoaffective disorder
#   - Schizophrenia, unspecified states

# FEP Codes for Case Definition
require(data.table)
broad_narrow_case_def = fread('FEP Case Def - individual codes 2021_0101_v4.txt')
hybrid_case_def = fread('FEP Case Def - hybrid with 95CI.txt')

# Broad-based definition - Use all codes in the preliminary case definition 
case_def_broad = broad_narrow_case_def$`FEP Code`

# Narrow-based definition - Use all codes in the preliminary case definition EXCEPT the 'MDD with psychotic behaviors' group
case_def_narrow = broad_narrow_case_def$`FEP Code`[broad_narrow_case_def$`Code Grouping` != 'Inclusion_DX_MDD_withPsychoticBehavior']

# Hybrid definition - Use codes in the preliminary case definition, EXCLUDING all individual codes whose PPV 95% confidence interval upper limits are below 0.7
case_def_hybrid = hybrid_case_def$`FEP Code`[hybrid_case_def$`95% CI Upper` > 0.7]


# Convert these codes into a regular expression that can be used to filter concept-codes
# [code1, code2, code3] --> '(^code1)|(^code2)|(^code3)'
case_def_broad_regexp = paste('(^', paste(case_def_broad, collapse = ')|(^'), ')',sep='')
case_def_narrow_regexp = paste('(^', paste(case_def_narrow, collapse = ')|(^'), ')',sep='')
case_def_hybrid_regexp = paste('(^', paste(case_def_hybrid, collapse = ')|(^'), ')',sep='')

# additional case-definitions (Tables 2-4 in the 'case definition' document)

# Table2 - Secondary Inclusion Codes
#   - Secondary Inclusion Codes: 295.01, 295.02, 295.03, 295.04, 295.11, 295.12, 295.13, 295.14, 295.21, 295.22, 295.23, 295.24, 295.31, 295.32, 295.33, 295.34, 295.41, 295.42, 295.43, 295.44, 295.51, 295.52, 295.53, 295.54, 295.81, 295.82, 295.83, 295.84
#   - Schizoaffective disorder, chronic or subchronic states: 295.71, 295.72, 295.73, 295.74
table2 = paste('ICD9', c('295.01', '295.02', '295.03', '295.04', '295.11', '295.12', '295.13', '295.14', '295.21', '295.22', '295.23', '295.24', '295.31', '295.32', '295.33', '295.34', '295.41', '295.42', '295.43', '295.44', '295.51', '295.52', '295.53', '295.54', '295.81', '295.82', '295.83', '295.84'), sep=':')
table2 = c(table2, paste('ICD9', c('295.71', '295.72', '295.73', '295.74'), sep=':'))
table2_regexp = paste('(^', paste(table2, collapse = ')|(^'), ')',sep='')

# Table 3 - Prior Psychosis Case Exclusions Codes
#   - Prior Psychosis Case Exclusions Codes: 295.05, 295.15, 295.25, 295.35, 295.45, 295.55, 295.75, 295.85, 295.95
#   - Residual Schizophrenia: 295.6
#   - Shared Psychotic Disorder: 298, F24*
table3 = paste('ICD9', c('295.05', '295.15', '295.25', '295.35', '295.45', '295.55', '295.75', '295.85', '295.95'), sep=':')
table3 = c(table3, paste('ICD9', c('295.6'), sep=':'))
table3 = c(table3, paste('ICD9', c('298'), sep=':'))
table3 = c(table3, paste('ICD10', c('F24.*'), sep=':'))
table3_regexp = paste('(^', paste(table3, collapse = ')|(^'), ')',sep='')

# Table 4 - Control Exclusion Codes
#   - Alcohol-induced psychosis: 291, 291.3, F10.15*, F10.25*, F10.95*
#   - Drug-induced psychosis: 292.1, 292.11, 292.12, F11.15*, F11.25*, F11.95*, F12.15*, F12.25*, F12.95*, F13.15*, F13.25*, F13.95*, F14.15*, F14.25*, F14.95*, F15.15*, F15.25*, F15.95*, F16.15*, F16.25*, F16.95*, F18.15*, F18.25*, F18.95*, F19.15*, F19.25*, F19.95*
table4 = paste('ICD9', c('291', '291.3', '292.1', '292.11', '292.12'), sep=':')
table4 = c(table4, paste('ICD10', c('F10.15.*', 'F10.25.*', 'F10.95.*', 'F11.15.*', 'F11.25.*', 'F11.95.*', 'F12.15.*', 'F12.25.*', 'F12.95.*', 'F13.15.*', 'F13.25.*', 'F13.95.*', 'F14.15.*', 'F14.25.*', 'F14.95.*', 'F15.15.*', 'F15.25.*', 'F15.95.*', 'F16.15.*', 'F16.25.*', 'F16.95.*', 'F18.15.*', 'F18.25.*', 'F18.95.*', 'F19.15.*', 'F19.25.*', 'F19.95.*'), sep=':'))
table4_regexp = paste('(^', paste(table4, collapse = ')|(^'), ')',sep='')


## Apply Case Definition  --------------------------------------------------------------------------

# ```{r apply_case_def, eval = APPLY_CASE_DEFINITION}
# We no longer need the original 'diagnosis' dataset (including all subjects) so for simplicity
# we shall only keep one diagnosis dataset
diagnosis <- diagnosis.included
rm(diagnosis.included)

# first - add a flag to each DX code to see if it meets the case-defs in tables 1-4
diagnosis$table1_narrow = F; diagnosis$table1_broad = F; diagnosis$table1_hybrid = F;  
diagnosis$table2 = F; diagnosis$table3 = F; diagnosis$table4 = F

diagnosis$table1_narrow[grep(case_def_narrow_regexp, diagnosis$DIAGNOSIS_CODE)] <- T
diagnosis$table1_broad[grep(case_def_broad_regexp, diagnosis$DIAGNOSIS_CODE)] <- T
diagnosis$table1_hybrid[grep(case_def_hybrid_regexp, diagnosis$DIAGNOSIS_CODE)] <- T

diagnosis$table2[grep(table2_regexp, diagnosis$DIAGNOSIS_CODE)] <- T
diagnosis$table3[grep(table3_regexp, diagnosis$DIAGNOSIS_CODE)] <- T
diagnosis$table4[grep(table4_regexp, diagnosis$DIAGNOSIS_CODE)] <- T
saveRDS(diagnosis, 'diagnosis.flagged.RDs')


## Apply Case Definition per Subject

# ```{r case_def_per_subject, eval = APPLY_CASE_DEFINITION}
getDxStartDate <- function(case_def, diagnosis){
  res = diagnosis[case_def, ]
  res = res[order(res$PATIENT_ID, res$START_DATE), ]
  res = res[!duplicated(res$PATIENT_ID), c('PATIENT_ID', 'START_DATE')]
  return(res)
}

# second, create a per-subject summary so that the case-definition could be applied
# a. get summary stats per subject
table1_narrow_count = tapply(diagnosis$table1_narrow, diagnosis$PATIENT_ID, sum)
table1_narrow_start = getDxStartDate(diagnosis$table1_narrow, diagnosis)

table1_broad_count = tapply(diagnosis$table1_broad, diagnosis$PATIENT_ID, sum)
table1_broad_start = getDxStartDate(diagnosis$table1_broad, diagnosis)

table1_hybrid_count = tapply(diagnosis$table1_hybrid, diagnosis$PATIENT_ID, sum)
table1_hybrid_start = getDxStartDate(diagnosis$table1_hybrid, diagnosis)

table2_count = tapply(diagnosis$table2, diagnosis$PATIENT_ID, sum)
table2_start = getDxStartDate(diagnosis$table2, diagnosis)

table3_count = tapply(diagnosis$table3, diagnosis$PATIENT_ID, sum)
table3_start = getDxStartDate(diagnosis$table3, diagnosis)

# b. create a patient-based data frame to hold these summary stats
case_def_table = data.frame(PATIENT_ID = demographics$PATIENT_ID, birth_date = demographics$BIRTH_DATE)

# table 1 narrow
case_def_table$table1_narrow_count = table1_narrow_count[as.character(case_def_table$PATIENT_ID)]
case_def_table$table1_narrow_start = table1_narrow_start$START_DATE[match(case_def_table$PATIENT_ID, table1_narrow_start$PATIENT_ID)]
case_def_table$table1_narrow_start_age = round(as.numeric(difftime(case_def_table$table1_narrow_start, case_def_table$birth_date, units = 'days') / 365), digits = 1)
# table 1 broad
case_def_table$table1_broad_count = table1_broad_count[as.character(case_def_table$PATIENT_ID)]
case_def_table$table1_broad_start = table1_broad_start$START_DATE[match(case_def_table$PATIENT_ID, table1_broad_start$PATIENT_ID)]
case_def_table$table1_broad_start_age = round(as.numeric(difftime(case_def_table$table1_broad_start, case_def_table$birth_date, units = 'days') / 365), digits = 1)
# table 1 hybrid
case_def_table$table1_hybrid_count = table1_hybrid_count[as.character(case_def_table$PATIENT_ID)]
case_def_table$table1_hybrid_start = table1_hybrid_start$START_DATE[match(case_def_table$PATIENT_ID, table1_hybrid_start$PATIENT_ID)]
case_def_table$table1_hybrid_start_age = round(as.numeric(difftime(case_def_table$table1_hybrid_start, case_def_table$birth_date, units = 'days') / 365), digits = 1)
# table 2
case_def_table$table2_count = table2_count[as.character(case_def_table$PATIENT_ID)]
case_def_table$table2_start = table2_start$START_DATE[match(case_def_table$PATIENT_ID, table2_start$PATIENT_ID)]
case_def_table$table2_start_age = round(as.numeric(difftime(case_def_table$table2_start, case_def_table$birth_date, units = 'days') / 365), digits = 1)
# table 3
case_def_table$table3_count = table3_count[as.character(case_def_table$PATIENT_ID)]
case_def_table$table3_start = table3_start$START_DATE[match(case_def_table$PATIENT_ID, table3_start$PATIENT_ID)]

saveRDS(case_def_table, 'case_def_table.RDs') # case_def_table = readRDS('case_def_table.RDs')


print(paste('Table 1 narraw:', sum(!duplicated(case_def_table$PATIENT_ID[case_def_table$table1_narrow_count >= 1])) ))
print(paste('Table 1 broad:', sum(!duplicated(case_def_table$PATIENT_ID[case_def_table$table1_broad_count >= 1])) ))
print(paste('Table 1 hybrid:', sum(!duplicated(case_def_table$PATIENT_ID[case_def_table$table1_hybrid_count >= 1])) ))
print(paste('Table 2:', sum(!duplicated(case_def_table$PATIENT_ID[case_def_table$table2_count >= 1])) ))
print(paste('Table 3:', sum(!duplicated(case_def_table$PATIENT_ID[case_def_table$table3_count >= 1])) ))


# ```{r load_case_def, eval = (!APPLY_CASE_DEFINITION)}
diagnosis = readRDS('diagnosis.flagged.RDs')
case_def_table = readRDS('case_def_table.RDs')


## Get Cases and Controls  --------------------------------------------------------------------------
# ```{r get_cases, eval = APPLY_CASE_DEFINITION}

# First, choose which type of outcome are we studying (we started with the 'broad' category)
SELECTED_OUTCOME = 'broad'

if (SELECTED_OUTCOME == 'narrow'){
  case_def_table$table1_count <- case_def_table$table1_narrow_count
  case_def_table$table1_start <- case_def_table$table1_narrow_start
  case_def_table$table1_start_age <- case_def_table$table1_narrow_start_age
}

if (SELECTED_OUTCOME == 'broad'){
  case_def_table$table1_count <- case_def_table$table1_broad_count
  case_def_table$table1_start <- case_def_table$table1_broad_start
  case_def_table$table1_start_age <- case_def_table$table1_broad_start_age
}

if (SELECTED_OUTCOME == 'hybrid'){
  case_def_table$table1_count <- case_def_table$table1_hybrid_count
  case_def_table$table1_start <- case_def_table$table1_hybrid_start
  case_def_table$table1_start_age <- case_def_table$table1_hybrid_start_age
}

# Cases: Patients with at least 2 FEP inclusion codes where first code occurs between age 15-35
cases = case_def_table[ !is.na(case_def_table$table1_start_age) & (case_def_table$table1_start_age >= 15 & case_def_table$table1_start_age < 35), ]
# Patients must have 2 Primary Inclusion Codes (Table 1) OR 1 Primary Inclusion Code and 1+ Secondary Inclusion Code (Table 2)
cases = cases[cases$table1_count >= 2 | (cases$table1_count == 1 & cases$table2_count >= 1), ]

# Patients with a prior exclusion code (Table 3) that occurs before a primary inclusion code are excluded
cases = cases[is.na(cases$table3_start) | (cases$table3_start >= cases$table1_start) , ]

# At least 2 encounters greater than 2yrs prior to index FEP code
cases_encounters = merge(cases[, c('PATIENT_ID', 'table1_start')], encounters, by = 'PATIENT_ID', all.x = T, all.y = F)
cases_prior_encounters = cases_encounters[cases_encounters$START_DATE < (cases_encounters$table1_start - (2*365*3600*24)), ]
cases_prior_encounters_count = tapply(cases_prior_encounters$ENCOUNTER_NUM, cases_prior_encounters$PATIENT_ID, function(x) sum(!duplicated(x)))
included_subjects_1 = names(cases_prior_encounters_count)[cases_prior_encounters_count >= 2]

# At least 1 encounter in the 10 years before index FEP code and index code must occur on or after year 2000
cases_prior_10y = cases_encounters[cases_encounters$START_DATE >= (cases_encounters$table1_start - (10*365*3600*24)), ]
cases_prior_10y_count = tapply(cases_prior_10y$ENCOUNTER_NUM, cases_prior_10y$PATIENT_ID, function(x) sum(!duplicated(x)))
included_subjects_2 = names(cases_prior_10y_count)[cases_prior_10y_count >= 1]
year2000 = strptime('01-01-2000', format = '%d-%m-%Y')
included_subjects_3 = cases$PATIENT_ID[cases$table1_start >= year2000]

# (optional: do not exclude, just flag inclusion/exclusion criteria)
cases$inclusion1_moreThanTwoYears = (cases$PATIENT_ID %in% included_subjects_1)
cases$inclusion2_visitInLast10Years = (cases$PATIENT_ID %in% included_subjects_2)
cases$inclusion3_indexAfterYear2000 = (cases$PATIENT_ID %in% included_subjects_3)
saveRDS(cases, paste('cases_noExclusions_',SELECTED_OUTCOME,'.RDs',sep='')) # cases <- readRDS(paste('cases_',SELECTED_OUTCOME,'.RDs',sep=''))

# ```{r load_cases, eval = (!APPLY_CASE_DEFINITION)}

SELECTED_OUTCOME = 'broad'
cases = readRDS(paste('cases_noExclusions_',SELECTED_OUTCOME,'.RDs',sep=''))



# ```{r load_ccs_codes}
icd9_to_ccs = read.csv('./CCS/ICD9 to CCS.csv')
icd10_to_ccs = read.csv('./CCS/ICD10 to CCS.csv')
icd9_ccs_dict = read.csv('./CCS/ICD9 CCS Dict.csv')

mental_icd9 = icd9_to_ccs$X.ICD.9.CM.CODE.[tolower(icd9_to_ccs$X.CCS.LVL.1.LABEL.) == 'mental illness']
mental_icd10 = icd10_to_ccs$X.ICD.10.CM.CODE.[grep('mental', tolower(icd10_to_ccs$X.CCS.CATEGORY.DESCRIPTION.))]
mental_codes = trimws(gsub("'", "", c(paste('ICD9',mental_icd9,sep = ':'), paste('ICD10',mental_icd10,sep=':')), fixed = T), which = 'both')


# ```{r get_controls}
# Controls: All patients with at least 2 visits between age 15 and 35
encounters$birth_date = demographics$BIRTH_DATE[match(encounters$PATIENT_ID, demographics$PATIENT_ID)]
encounters$age = floor(as.numeric(difftime(encounters$START_DATE, encounters$birth_date, units = 'days') / 365))

encounters_15_35 = encounters[!is.na(encounters$age) & (encounters$age >= 15 & encounters$age <= 35), ]
encounters_15_35_count = tapply(encounters_15_35$ENCOUNTER_NUM, encounters_15_35$PATIENT_ID, function(x) sum(!duplicated(x)))
included_subjects_1 = names(encounters_15_35_count)[encounters_15_35_count >= 2]

# No psychosis related codes (TABLES 1-4)
excluded_ids = unique(diagnosis$PATIENT_ID[diagnosis$table1_broad | diagnosis$table2 | diagnosis$table3 | diagnosis$table4])

# At least 2 visits greater than 2 years apart
first_encounter = encounters[!duplicated(encounters$PATIENT_ID), ]
last_encounter = encounters[rev(!duplicated(rev(encounters$PATIENT_ID))), ]
diff_between_encounters = floor(as.numeric(difftime(last_encounter$START_DATE, first_encounter$START_DATE, units = 'days') / 365))
included_subjects_2 = first_encounter$PATIENT_ID[diff_between_encounters > 2]

controls_ids = demographics$PATIENT_ID[(demographics$PATIENT_ID %in% included_subjects_1) & 
                                         (demographics$PATIENT_ID %in% included_subjects_2) & 
                                         !(demographics$PATIENT_ID %in% excluded_ids) & 
                                         !(demographics$PATIENT_ID %in% cases$PATIENT_ID)]

# NEW: only include controls with some mental DXs
controls_dxs = diagnosis[diagnosis$PATIENT_ID %in% controls_ids, ]
controls_dxs$DIAGNOSIS_CODE_noDot = gsub('.', '', controls_dxs$DIAGNOSIS_CODE, fixed = T) 
controls_dxs$isMentalCode = controls_dxs$DIAGNOSIS_CODE_noDot %in% mental_codes
controls_with_mental_codes = unique(controls_dxs$PATIENT_ID[controls_dxs$isMentalCode])
controls_ids = controls_ids[controls_ids %in% controls_with_mental_codes]

controls = demographics[demographics$PATIENT_ID %in% controls_ids, ]
saveRDS(controls, 'controls.RDs') 

cases_and_controls_ids = c(cases$PATIENT_ID, controls$PATIENT_ID)
saveRDS(cases_and_controls_ids, 'cases_and_controls_ids.RDs')


## Import Txt files into RDs (Part 2)

# 
# Load the following datasets only for included subjects
#
included_ids = readRDS('cases_and_controls_ids.RDs')
# labs
labs = readAllDatasetSegments('./Labs2', included_ids)
labs$START_DATE <- as.POSIXct(strptime(labs$START_DATE, format = '%Y-%m-%d %H:%M:%OS'))
saveRDS(labs, 'labs.RDs')


# medications
meds = readAllDatasetSegments('./Medications2', included_ids)
meds$START_DATE = as.POSIXct(strptime(meds$START_DATE, format = '%Y-%m-%d %H:%M:%OS'))
saveRDS(meds, 'meds.RDs') 

# procedures
procs = readAllDatasetSegments('./Procedures', included_ids)
procs$START_DATE = as.POSIXct(strptime(procs$START_DATE, format = '%Y-%m-%d %H:%M:%OS'))
saveRDS(procs, 'procs.RDs') 


# encounters


