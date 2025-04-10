#
# Notes: 
# The partners dataset includes a total of 2.7M patients. 1.3M patients born between study date range 1965-2005
# Only selecting patients that meet either cases or controls case definition, we get 643,984 subjects. 
# The folder 'sample2' contains these 643,984 subjects that were born between 1965-2005 and meet either case or controls case definition. 
#


.libPaths("~/PNGU/R_PACKAGES/")

install.packages('lubridate')
install.packages('ggplot2')
install.packages('dplyr')
install.packages('caret')
install.packages('reshape2')
install.packages('pillar')
install.packages('devtools')
install.packages('pROC')
devtools::install_version("randomForest", version = "4.6.14", repos = "http://cran.us.r-project.org")


require(lubridate)
require(ggplot2)


## Define Global Functions
# Functions Decleration  --------------------------------------------------------------------------
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
showIQR <- function(x, n_digits = 1){
  a = round(quantile(x, probs = c(0.25,0.75)), digits = n_digits)
  res = paste(a['25%'],'-',a['75%'],sep='')
  return(res)
}
showMedianAndIQR <- function(x, n_digits = 1){
  IQR = showIQR(x, n_digits)
  med = round(median(x), digits = n_digits)
  return(paste(med,' [',IQR,']',sep = ''))
}
showMedianAndRange <- function(x, n_digits = 1){
  range = paste(round(c(min(x),max(x)), n_digits), collapse = '-')
  med = round(median(x), digits = n_digits)
  return(paste(med,' [',range,']',sep = ''))
}

showMeanAnd95CI <- function(x, n_digits = 1){
  tt = t.test(x)
  ci = paste(round(c(tt$conf.int[1],tt$conf.int[2]), n_digits), collapse = '-')
  ave = round(tt$estimate, digits = n_digits)
  return(paste(ave,' [',ci,']',sep = ''))
}


showCountAndPercent <- function(count, total,n_digits = 1){
  prcnt = round(100 * count / total, digits = n_digits)
  res = paste(prcnt, '% (n=',count,')',sep = '')
  return(res)
}
getBirthYear <- function(dates_str){
  if(length(grep('POSIX', class(dates_str))) > 0){
    res = dates_str$year + 1900
  } else {
    res = as.numeric(gsub('.+/.+/','',dates_str))  
  }
  return(res)
}

getFactorSummaryString <- function(a){
  a = factor(as.character(a))
  s = summary(a); s = s[order(s, decreasing = T)]
  res = paste(paste(paste(names(s), round(100*s / length(a), digits = 1), sep = '='), collapse = '%, '), '%',sep='')
  return(res)
}

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

# postfix = 'AllControls' | 'SelectedControls' 
getModelPerformance <- function(postfix, cases, included_ids)
{
  res = data.frame()
  for (PRED_AGE in c(15,20,25,30)){
    print(paste('------------------------ Age =',PRED_AGE,'------------------------'))
    for (FOLLOW_UP in c(1,5,NA)){
      rf_fit_filename = paste('rf.fit2.follow_up_',FOLLOW_UP,'_age_',PRED_AGE,'_',postfix,'.RDs',sep='')
      if(file.exists(rf_fit_filename)){
        print(paste('------------------------ Follow-Up =',FOLLOW_UP,'------------------------'))
        testing <- readRDS(paste('testing.followup_',FOLLOW_UP,'_age_',PRED_AGE,'_','AllControls2','.RDs',sep=''))
        training <- readRDS(paste('training.followup_',FOLLOW_UP,'_age_',PRED_AGE,'_','AllControls2','.RDs',sep=''))
        rf.fit = readRDS(rf_fit_filename) 
        
        if (length(included_ids) > 0){
          training <- training[rownames(training) %in% as.character(included_ids), ]
          testing <- testing[rownames(testing) %in% as.character(included_ids), ]
        }
        
        imp_vars = varImpPlot(rf.fit, main = paste('Prediction age:',PRED_AGE))
        imp_vars <- rownames(imp_vars)[order(imp_vars, decreasing = T)]
        print(paste(imp_vars[1:10], collapse = ', '))
        
        # make sure that the testing set doesn't include information not available before for the model
        missing_lvls = levels(testing$current_zip3)[!(levels(testing$current_zip3) %in% levels(training$current_zip3))]
        if(length(missing_lvls) > 0){
          testing$current_zip3[testing$current_zip3 %in% missing_lvls] <- 'Other'
          testing$current_zip3 <- factor(as.character(testing$current_zip3))
        }
        
        # apply model on testing set
        rf.pred = predict(rf.fit, testing, type = 'prob')
        roc.res = roc(testing$isCase, rf.pred[,2])
        print(roc.res$auc)
        stats_tbl = getSummaryStatsTbl(roc.res)
        
        for(specificity in c(0.9, 0.95, 0.99)){
          # get subjects identified by the model (score > cutoff)
          cutoff_score = quantile(rf.pred[,2][testing$isCase == F], probs = specificity)
          identified_cases = rownames(testing)[which(testing$isCase == T & rf.pred[,2] > cutoff_score)]
          # get the number of years from prediction till index-event
          time_till_index = cases$table1_start_age[match(identified_cases, cases$subject_num)] - PRED_AGE
          
          res = rbind(res, data.frame(ModelType = postfix, 
                                      TimeWindow = FOLLOW_UP,
                                      Age = PRED_AGE,
                                      Specificity = specificity,
                                      Cases = showCountandRate(sum(testing$isCase == 'TRUE'), nrow(testing)),
                                      Controls = showCountandRate(sum(testing$isCase == 'FALSE'), nrow(testing)),
                                      AUC = round(roc.res$auc,digits = 2),
                                      PPV = stats_tbl$ppv[stats_tbl$specificity == specificity],
                                      NPV = stats_tbl$npv[stats_tbl$specificity == specificity],
                                      TimeTillIndex = showMedianAndIQR(time_till_index)))
          
        }        
      } else {
        res = rbind(res, data.frame(ModelType = postfix, TimeWindow = FOLLOW_UP, Age = PRED_AGE,
                                    Specificity = NA,
                                    Cases = NA,
                                    Controls = NA,
                                    AUC = NA,
                                    PPV = NA,
                                    NPV = NA,
                                    TimeTillIndex = NA))
        
      }
    }
  } 
  return(res)
}


getModelPerformance2 <- function(postfix, cases, included_ids)
{
  res = data.frame()
  for (PRED_AGE in c(15,20,25,30)){
    print(paste('------------------------ Age =',PRED_AGE,'------------------------'))
    for (FOLLOW_UP in c(3)){ # 1,5,NA
      print(paste('------------------------ Follow-Up =',FOLLOW_UP,'------------------------'))
      # load files (testing cohort + RF model)
      testing <- readRDS(paste('testing.followup_',FOLLOW_UP,'_age_',PRED_AGE,'_',postfix,'.RDs',sep=''))
      rf.fit = readRDS(paste('rf.fit2.follow_up_',FOLLOW_UP,'_age_',PRED_AGE,'_',postfix,'.RDs',sep='')) 
      testing <- testing[rownames(testing) %in% as.character(included_ids), ]
      
      # apply model on testing set
      rf.pred = predict(rf.fit, testing)
      pred_score = rf.pred$predicted[,'TRUE']
      roc.res = roc(rf.pred$yvar == 'TRUE', pred_score)
      stats_tbl = getSummaryStatsTbl(roc.res)
      
      for(specificity in c(0.9, 0.95, 0.99)){
        # get subjects identified by the model (score > cutoff)
        cutoff_score = quantile(pred_score[testing$isCase == F], probs = specificity)
        identified_cases = rownames(testing)[which(testing$isCase == T & pred_score > cutoff_score)]
        # get the number of years from prediction till index-event
        time_till_index = cases$table1_start_age[match(identified_cases, cases$subject_num)] - PRED_AGE
        
        res = rbind(res, data.frame(ModelType = postfix, 
                                    TimeWindow = FOLLOW_UP,
                                    Age = PRED_AGE,
                                    Specificity = specificity,
                                    Cases = showCountandRate(sum(testing$isCase == 'TRUE'), nrow(testing)),
                                    Controls = showCountandRate(sum(testing$isCase == 'FALSE'), nrow(testing)),
                                    AUC = round(roc.res$auc,digits = 2),
                                    PPV = stats_tbl$ppv[stats_tbl$specificity == specificity],
                                    NPV = stats_tbl$npv[stats_tbl$specificity == specificity],
                                    TimeTillIndex = showMedianAndIQR(time_till_index)))
        
      }        
    }
  } 
  return(res)
}


setwd("/data/js95/PNGU/RPDRml_v2021.1/sample2")
## Load Data

# the subjects meeting the specific controls case definition
cases = readRDS('./cases.RDs')
controls = readRDS('controls_sample2.RDs') 
controls.sample = readRDS('controls_subsample.RDs') # down-sample controls to include x4 the number of cases

# the list of index dates (index event for cases, last documented encounter for controls)
index_events = readRDS('index_events.RDs')

# get all available subjects from demographics dataset
all_subjects = readRDS('../demographics.RDs')
all_subjects = all_subjects[!duplicated(all_subjects$subject_num), ]

# subjects included for models
included_subjects = all_subjects[all_subjects$subject_num %in% c(cases$subject_num, controls$subject_num), ]
included_subjects.sample = all_subjects[all_subjects$subject_num %in% c(cases$subject_num, controls.sample$subject_num), ]
cases_dem = included_subjects[included_subjects$subject_num %in% cases$subject_num, ]

controls_dem = included_subjects[included_subjects$subject_num %in% controls$subject_num, ]



# summary stats for included cohort

# look at the distribution of cases over time
cases_months = as.integer(format(cases$index_date, format = '%m'))
cases_years = as.integer(format(cases$index_date, format = '%Y'))
a = tapply(cases$subject_num, paste(1,cases_months,cases_years,sep='/'), length)
plot_data = data.frame(x = strptime(names(a), format = '%d/%m/%Y'), y = a)
barplot(plot_data$y ~ plot_data$x)

# examine controls birth date vs. cases
require()
plot_data = data.frame(x = c(cases$birth_date, controls$sbirth_date), y = c(rep('Case',nrow(cases)), rep('Control',nrow(controls))))
ggplot(data = plot_data, aes(x = x, fill = y)) + 
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', aes(y=..density..)) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_minimal() +
  labs(fill="")



## Analysis of available prior Encounters
# summary of encounters prior to FEP - 17,886,693 x 3 (subject_num, encounter_num, start_date)
encounters = readRDS('encounters2.RDs')
encounters = encounters[encounters$subject_num %in% included_subjects$subject_num, ]

# create unique encounter identification number for encounters with >24 hours between them
encounters = encounters[order(encounters$subject_num, encounters$sstart_date), ]
encounters$time_from_last_encounter = c(NA, difftime(encounters$sstart_date[-1], encounters$sstart_date[-nrow(encounters)], units = 'hours'))
id_change = c(TRUE, encounters$subject_num[-1] != encounters$subject_num[-nrow(encounters)])
encounters$time_from_last_encounter[which(id_change)] <- NA
encounters$unique_encounter = is.na(encounters$time_from_last_encounter) | (encounters$time_from_last_encounter > 24)


cases_encounters = merge(cases[, c('subject_num', 'table1_start')], encounters, by = 'subject_num', all.x = T, all.y = F)
cases_encounters = cases_encounters[order(cases_encounters$subject_num, cases_encounters$sstart_date), ]
cases_prior_encounters = cases_encounters[cases_encounters$sstart_date < (cases_encounters$table1_start - (2*365*3600*24)), ]
# option 1: look at all encounters, only discriminating by encounter num
cases_prior_encounters_count = tapply(cases_prior_encounters$subject_num, cases_prior_encounters$subject_num, length)
# option 2: treat encounters as the same as long as there is <=24 hours between them
cases_prior_encounters_count = tapply(cases_prior_encounters$unique_encounter, cases_prior_encounters$subject_num, sum)
summary(cases_prior_encounters_count)

# investigate the outliers: 
# ids = names(cases_prior_encounters_count)[which(cases_prior_encounters_count > 100)]
# for(id in sample(ids, 1)){
#   d = cases_encounters[cases_encounters$subject_num == id, ]
#   print(d[order(d$START_DATE), ])
# }


## Create Table-1

# 2. Get Stats on Cases and Controls  --------------------------------------------------------------------------
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# TODO: double check feature names vs. Galina's data
createTable1 <- function(all_subjects, included_subjects, cases_dem, controls_dem)
{
  # filter to not include more than the patients in 'all_subjects'
  included_subjects = included_subjects[included_subjects$subject_num %in% all_subjects$subject_num, ]
  cases_dem = cases_dem[cases_dem$subject_num %in% all_subjects$subject_num, ]
  controls = controls[controls$subject_num %in% all_subjects$subject_num, ]
  
  table_1 = data.frame()
  table_1 = rbind(table_1, data.frame(feature = 'Count', All = nrow(all_subjects), Included = nrow(included_subjects), Cases = nrow(cases_dem), Controls = nrow(controls)))
  table_1 = rbind(table_1, data.frame(feature = 'Sex', All = getFactorSummaryString(all_subjects$SEX_CD), Included = getFactorSummaryString(included_subjects$SEX_CD), Cases = getFactorSummaryString(cases_dem$SEX_CD), Controls = getFactorSummaryString(controls_dem$SEX_CD)))
  table_1 = rbind(table_1, data.frame(feature = 'Birth Year', All = median(getBirthYear(all_subjects$sbirth_date)), Included = median(getBirthYear(included_subjects$sbirth_date)), Cases = median(getBirthYear(cases_dem$sbirth_date)), Controls = median(getBirthYear(controls_dem$sbirth_date))))
  table_1 = rbind(table_1, data.frame(feature = 'Race', All = getFactorSummaryString(all_subjects$race), Included = getFactorSummaryString(included_subjects$race), Cases = getFactorSummaryString(cases_dem$race), Controls = getFactorSummaryString(controls$race)))
  table_1 = rbind(table_1, data.frame(feature = 'Ethnicity', All = getFactorSummaryString(all_subjects$ethnicity), Included = getFactorSummaryString(included_subjects$ethnicity), Cases = getFactorSummaryString(cases_dem$ethnicity), Controls = getFactorSummaryString(controls$ethnicity)))
  table_1 = rbind(table_1, data.frame(feature = 'Marital Status', All = getFactorSummaryString(all_subjects$MARITAL_STATUS_CD), Included = getFactorSummaryString(included_subjects$MARITAL_STATUS_CD), Cases = getFactorSummaryString(cases_dem$MARITAL_STATUS_CD), Controls = getFactorSummaryString(controls_dem$MARITAL_STATUS_CD)))
  #summary(factor(all_subjects$ethnicity))  
  
  return(table_1)
}


table_1 = createTable1(all_subjects, included_subjects, cases_dem, controls_dem)
write.csv(table_1, 'table1.csv')

# extract p-values for table-1
case_def = factor(c(rep('all', nrow(included_subjects)), rep('cases', nrow(cases_dem))))
chisq.test( factor(c(as.character(included_subjects$SEX_CD), as.character(cases_dem$SEX_CD))), case_def )
# chisq.test( factor(c(as.character(included_subjects$race), as.character(cases_dem$race))), case_def )
chisq.test( factor(c(as.character(included_subjects$MARITAL_STATUS_CD), as.character(cases_dem$MARITAL_STATUS_CD))), case_def )

# age of FEP
s = summary(factor(round(cases$table1_start_age, digits = 0)))
round(s / sum(s), digits = 2)

## Table-1 for manuscript: only includes patients that are eligble for our models 
# (within the right age-group of 15-35 and with enough follow-up)

# get first & last encounter for each subject
encounters = encounters[order(encounters$subject_num, encounters$sstart_date), ]
# first, add age to encounters
all_subjects$sbirth_date = strptime(all_subjects$sbirth_date, format = '%m/%d/%Y')
dob = all_subjects$sbirth_date[match(encounters$subject_num, all_subjects$subject_num)]
encounters$age = pmax(round((as.numeric(difftime(encounters$sstart_date, dob, units = 'days'))/365.25), digits = 1), 0)
# now get the 1st and last encounter age
first_encounter = encounters[!duplicated(encounters$subject_num), ]
last_encounter = encounters[rev(!duplicated(rev(encounters$subject_num))), ]
age_range = merge(first_encounter[, c('subject_num', 'age')], last_encounter[, c('subject_num', 'age')], by='subject_num')
colnames(age_range)[2:3] <- c('first_age', 'last_age')
age_range = age_range[!is.na(age_range$first_age) & !is.na(age_range$last_age), ]

saveRDS(age_range, 'range_of_encounters.RDs') # age_range = readRDS('range_of_encounters.RDs')

age_range$isCase = age_range$subject_num %in% cases$subject_num
age_range$indexEvent = cases$index_age[match(age_range$subject_num, cases$subject_num)]

all_potential_subjects = c()
for(pred_age in c(15,20,25,30)){
  potential_subjects = age_range$subject_num[age_range$first_age < pred_age & age_range$last_age > (pred_age+5)]
  # cases can have shorter follow-up if they meet the case definition within the time window 
  # (we need enough follow-up to either rule-in or rule-out cases)
  cases_age_range = age_range[age_range$isCase, ]
  potential_cases = cases_age_range$subject_num[cases_age_range$first_age < pred_age & cases_age_range$indexEvent > (pred_age) & cases_age_range$indexEvent < (pred_age+5)]
    
  # add these cases to all potential subjects (i.e. make sure that we include all cases)
  potential_subjects = unique(c(potential_subjects, potential_cases))
  
  print(paste('Age=',pred_age,', subjects=',length(potential_subjects)))
  table1_sub = createTable1(all_subjects[all_subjects$subject_num %in% potential_subjects, ], included_subjects, cases_dem, controls_dem)
  
  write.csv(table1_sub, paste('table1_age_',pred_age,'.csv',sep=''))
  
  all_potential_subjects = unique(c(all_potential_subjects, potential_subjects))
}
table1_all = createTable1(all_subjects[all_subjects$subject_num %in% all_potential_subjects, ], included_subjects, cases_dem, controls_dem)

write.csv(table1_all, 'table1_all_ages.csv')



## Table 2
# 3. Get number of available patients for each model (Table 2)  -------------------------------------
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# this needs to be re-written - I currently do not save this csv table

# table2 <- read.csv('./all_ages_stats.csv') # file name is misleading but this includes only true-controls
# print(table2)
# print(sum(as.numeric(gsub(' .+', '', table2$cases[table2$follow_up == 5]) )))
```


##  Run Models 

## Naive Bayes Model version 1 - for any time in the future 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
require(pROC, lib.loc = '/data/js95/PNGU/R_PACKAGES/')
require(ggplot2, lib.loc = '/data/js95/PNGU/R_PACKAGES/')

addNBCScores <- function(pred.data, PRED_AGE){
  # load odds-ratios, and replace infinite ORs (patients with 0 controls) with the max OR value 
  odds_ratios = readRDS(paste('./model_datasets/odds_ratio.age_',PRED_AGE,'.RDs',sep=''))
  odds_ratios$OR[is.infinite(odds_ratios$OR)] <- max(odds_ratios$OR[!is.infinite(odds_ratios$OR)], na.rm = T)
  # replace OR=0 with minimum OR which is > 0
  odds_ratios$OR[odds_ratios$OR == 0] <- min(odds_ratios$OR[odds_ratios$OR > 0], na.rm = T)
  
  # keep only one occurrence from each concept code
  pred.data = pred.data[order(pred.data$subject_num, pred.data$date), ]
  pred.data.unique = pred.data[!duplicated(paste(pred.data$subject_num, pred.data$concept_type, pred.data$concept_code,sep = '_')), ]
  
  
  # add NBC score to each row in prediction data
  pred.data.unique$NBC = log(odds_ratios$OR[match(paste(pred.data.unique$concept_type, pred.data.unique$concept_code,sep = '_'), odds_ratios$concept_code)])
  pred.data.unique$NBC[is.na(pred.data.unique$NBC)] <- 0
  
  # calc cumulative NBC over time
  pred.data.unique$sumNBC = unlist(tapply(pred.data.unique$NBC, pred.data.unique$subject_num, cumsum))
  
  # add the rounded age in years (rounded up)
  pred.data.unique$age_rounded = ceiling(pred.data.unique$age)
  
  saveRDS(pred.data.unique, paste('./model_datasets/nbc_scored.age_',PRED_AGE,'.RDs',sep=''))  
  
  return(pred.data.unique)
}


PRED_AGE = 20

for(PRED_AGE in c(15,20,25,30)){
  print(paste('pred age:',PRED_AGE))

  pred.data = readRDS(paste('./model_datasets/pred.data.age_',PRED_AGE,'.RDs',sep=''))
  pred.data.scored = addNBCScores(pred.data, PRED_AGE) # this also saves an RDs with the scored file
  
  # load already scored NBC scores (if not 1st time running uncomment next line and comment above 2 lines)
  # pred.data.scored = readRDS(paste('./nbc_scored.age_',PRED_AGE,'.RDs',sep=''))
  
  total_cases = sum(!duplicated(pred.data.scored$subject_num[pred.data.scored$isCase == T]))
  total_controls = sum(!duplicated(pred.data.scored$subject_num[pred.data.scored$isCase == F]))
  print(paste('Pred age',PRED_AGE, 'Cases=',total_cases,'Controls=',total_controls, 'Prcnt cases=', round(100 * (total_cases / (total_cases + total_controls)), digits = 1),'%' ))
  
  # plot NBC score by age for cases vs. controls
  # get only testing set and only last score for each year of age (rounding up the ages for this)
  plot_data = pred.data.scored[(pred.data.scored$isTraining == F) & !duplicated(paste(pred.data.scored$subject_num, pred.data.scored$age_rounded,sep='_')), ]
  plot_data = plot_data[plot_data$age > 0, ]
  
  # number of patients in testing set
  total_count =sum(!duplicated(plot_data$subject_num))
  p = ggplot(data = plot_data, aes(x = age_rounded, y = sumNBC, color = isCase, fill = isCase)) +
    stat_summary(geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=0.95), alpha = 0.5)+ # fill="lightblue"
    stat_summary(geom="line", fun.y=mean, linetype="dashed") +
    stat_summary(geom="point", fun.y=mean) + # color="red"
    theme_minimal() + xlab('Age (Years)') + ylab('NBC Score') + ggtitle(paste('Age at prediction =',PRED_AGE,'years (n=',total_count,')')) + 
    scale_color_manual(values=c("blue", "red")) + scale_fill_manual(values=c("lightblue", "lightcoral")) + ylim(-100,50)
  print(p)  
  
  last_nbc_score = pred.data.scored[rev(!duplicated(rev(pred.data.scored$subject_num))), ]
  last_nbc_score = last_nbc_score[last_nbc_score$isTraining == F, ]
  roc.res = roc(last_nbc_score$isCase ~ last_nbc_score$sumNBC)
  summary_stats = getSummaryStatsTbl(roc.res)
  write.csv(summary_stats, paste('./model_results/summary_stats_age_',PRED_AGE,'.csv',sep=''))
  print(roc.res)
  print(summary_stats)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



## Naive Bayes Model version 2 - for 5 years time-window
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
require(pROC)


PRED_AGE = 20
age_range = readRDS('range_of_encounters.RDs')

for(PRED_AGE in c(15,20,25,30)){
  print(paste('pred age:',PRED_AGE))
  
  pred.data = readRDS(paste('./model_datasets/pred.data.age_',PRED_AGE,'.RDs',sep=''))
  
  # Only include 5Y time-window for prediction
  potential_subjects = age_range$subject_num[age_range$first_age < PRED_AGE & age_range$last_age > (PRED_AGE+5)]
  cases_to_include = cases$subject_num[!is.na(cases$index_age) & (cases$index_age > PRED_AGE) & (cases$index_age < (PRED_AGE+5))]
  potential_subjects = unique(c(potential_subjects, cases_to_include))
  
  pred.data = pred.data[pred.data$subject_num %in% potential_subjects, ]
  pred.data$isCase = (pred.data$subject_num %in% cases_to_include)
  
  pred.data.scored = addNBCScores(pred.data, PRED_AGE) # this also saves an RDs with the scored file
  
  # load already scored NBC scores (if not 1st time running uncomment next line and comment above 2 lines)
  # pred.data.scored = readRDS(paste('./nbc_scored.age_',PRED_AGE,'.RDs',sep=''))
  
  total_cases = sum(!duplicated(pred.data.scored$subject_num[pred.data.scored$isCase == T]))
  total_controls = sum(!duplicated(pred.data.scored$subject_num[pred.data.scored$isCase == F]))
  print(paste('Pred age',PRED_AGE, 'Cases=',total_cases,'Controls=',total_controls, 'Prcnt cases=', round(100 * (total_cases / (total_cases + total_controls)), digits = 1),'%' ))
  
  last_nbc_score = pred.data.scored[rev(!duplicated(rev(pred.data.scored$subject_num))), ]
  last_nbc_score = last_nbc_score[last_nbc_score$isTraining == F, ]
  roc.res = roc(last_nbc_score$isCase ~ last_nbc_score$sumNBC)
  summary_stats = getSummaryStatsTbl(roc.res)
  write.csv(summary_stats, paste('summary_stats_NBC_timewindow5Y_age_',PRED_AGE,'.csv',sep=''))
  print(roc.res)
  print(summary_stats)
}



# apply score on testing set, this time including all controls - not only down-sampled cohort


## Get top predictors

# # Create dictionaries for DXs, meds, labs
# dx_dict = readRDS('./diagnosis.included.RDs')
# dx_dict = dx_dict[!duplicated(dx_dict$DIAGNOSIS_CODE), ]
# saveRDS(dx_dict, './dx_dict.RDs')
# 
# labs_dict = readRDS('./labs.RDs')
# labs_dict = labs_dict[!duplicated(labs_dict$LAB_CODE), ]
# saveRDS(labs_dict, './labs_dict.RDs')
# 
# meds_dict = readRDS('./meds.RDs')
# meds_dict = meds_dict[!duplicated(meds_dict$MED_CODE), ]
# saveRDS(meds_dict, './meds_dict.RDs')
# 
# library(data.table)
# dx_dict = readRDS('../RPDRml__DX_dictionary.RDs')
# labs_dict = fread('../RPDRml__LAB_dictionary.txt', select = c('concept_code', 'concept_name'))
# meds_dict = fread('../RPDRml__MED_dictionary.txt', select = c('concept_code', 'concept_name'))
# labs_dict = fread('../RPDRml__LAB_dictionary.txt', select = c('concept_code', 'concept_name'))
# meds_dict = fread('../RPDRml__MED_dictionary.txt', select = c('concept_code', 'concept_name'))


dx_dict = readRDS('../RPDRml__DX_dictionary.RDs')
procs_dict = read.delim('../RPDRml__PROC_dictionary.txt', stringsAsFactors = F)
labs_dict = read.delim('../RPDRml__LAB_dictionary.txt', stringsAsFactors = F)
meds_dict = read.delim('../RPDRml__MED_dictionary.txt', stringsAsFactors = F)


# add concept name to the odds ratio table
for(PRED_AGE in c(15,20,25,30)){
  print(paste('PRED AGE', PRED_AGE))
  odds_ratios = readRDS(paste('./model_datasets/odds_ratio.age_',PRED_AGE,'.RDs',sep=''))
  
  dx_idxs = grep('dx_', odds_ratios$concept_code)
  labs_idxs = grep('labs_', odds_ratios$concept_code)
  meds_idxs = grep('meds_', odds_ratios$concept_code)
  
  odds_ratios$description = ''
  odds_ratios$description[dx_idxs] = dx_dict$concept_name[match(gsub('dx_','',odds_ratios$concept_code[dx_idxs]), dx_dict$concept_code)]
  odds_ratios$description[labs_idxs] = labs_dict$concept_name[match(gsub('labs_','',odds_ratios$concept_code[labs_idxs]), labs_dict$concept_code)]
  odds_ratios$description[meds_idxs] = meds_dict$concept_name[match(gsub('meds_','',odds_ratios$concept_code[meds_idxs]), meds_dict$concept_code)]
  saveRDS(odds_ratios, paste('./model_datasets/odds_ratio.age_',PRED_AGE,'.RDs',sep=''))
  
  print(odds_ratios[1:20, ])
  write.csv(odds_ratios[1:20, ], paste('./model_results/top_risk_factors.age_',PRED_AGE,'.csv',sep=''))
}

# debug cheat codes
pred.data.scored = readRDS(paste('./nbc_scored.age_',25,'.RDs',sep=''))

ids = unique(pred.data.scored$subject_num[pred.data.scored$concept_code == 'ICD9:298.9'])

View(pred.data.scored[pred.data.scored$subject_num %in% ids[1:5], ])






#
# Run Random Forest classification algorithm
#
reduceNumberOfFactorLevels <- function(fctr, num_levels){
  s = summary(fctr)
  lvls_to_include = names(s)[1:num_levels]
  fctr = as.character(fctr)
  fctr[!(fctr %in% lvls_to_include)] <- 'Other'
  fctr = factor(fctr)
  
  return(fctr)
}

library(caret)
library(randomForest)

# Load relevant files
PRED_AGE = 20
setwd("/data/js95/PNGU/RPDRml_v2021.1/sample2/")
pred.data.flat = readRDS(paste('./model_datasets/pred.data.flat.age_',PRED_AGE,'.RDs',sep=''))  

# some clean-up of the pred dataset
colnames(pred.data.flat) <- make.names(colnames(pred.data.flat))
pred.data.flat = pred.data.flat[, !duplicated(colnames(pred.data.flat))]

training = pred.data.flat[pred.data.flat$isTraining == 1, ]
testing = pred.data.flat[pred.data.flat$isTraining == 0, ]

cases_idxs = which(training$isCase == 1)
controls_idxs = which(training$isCase == 0)
controls_idxs_sample = sample(controls_idxs, length(cases_idxs) * 4, replace = F)
training.sample = training[c(cases_idxs, controls_idxs_sample), ]


# Assuming that "target" is the name of your outcome variable

# Define the control method for train() function
# This will use 5-fold cross-validation and down-sampling to handle class imbalance
train_control <- trainControl(method="cv", 
                              number=5, 
                              verboseIter=TRUE, 
                              returnData=FALSE, 
                              returnResamp="all", 
                              classProbs = TRUE,
                              sampling = "down",
                              search = "random",
                              summaryFunction=twoClassSummary)

# Set up a simple grid of tuning parameters to try
# Just using mtry (number of variables randomly sampled at each split)
grid <- expand.grid(.mtry = c(2, 4, 6, 10))

# Fit the model on your data
# Note: replace `your_data` with your actual data frame, and `target` with your actual target variable
set.seed(123)
levels(training.sample$isCase) <- c('Control','Case')
rf_model <- train(isCase ~ ., 
                  data = training.sample, 
                  method = "rf", 
                  metric = "ROC",
                  tuneGrid = grid, 
                  trControl = train_control)

# Get feature importance
importance <- varImp(rf_model, scale = FALSE)
# Get the names of the top 500 variables
topvars = rownames(importance$importance)[order(importance$importance, decreasing = T)[1:500]]
# remove 'cheat codes' 
topvars = c(topvars[!(topvars %in% 'timeTillEvent')], 'isCase')

# This will give you the names of the factor columns
factor_columns_names <- names(training.sample)[sapply(training.sample, is.factor)]

# for categorical values - we don't want to include specific values, but the entire variable
# [1] "SEX"            "MARITAL_STATUS" "race"           "ethnicity"      "veteran"        "public_payer"  
topvars <- sub("^SEX.*", "SEX", topvars)
topvars <- sub("^MARITAL_STATUS.*", "MARITAL_STATUS", topvars)
topvars <- sub("^race.*", "race", topvars)
topvars <- sub("^ethnicity.*", "ethnicity", topvars)
topvars <- sub("^veteran.*", "veteran", topvars)
topvars <- sub("^public_payer.*", "public_payer", topvars)
topvars = unique(topvars)

# Create a new data frame that only includes the top variables
new_data.training <- training.sample[, topvars]
new_data.testing <- testing[, topvars]

# test the model using only top variables
rf_model <- train(isCase ~ ., 
                  data = new_data.training, 
                  method = "rf", 
                  metric = "ROC",
                  tuneGrid = grid, 
                  trControl = train_control)

y.pred <- predict(rf_model, newdata = new_data.testing, type = 'prob')


roc.res = roc(new_data.testing$isCase, y.pred[,2])
View(getSummaryStatsTbl(roc.res))























## Random Forest Survival Model
```{r random_forest_survival_model, echo = FALSE}
require(randomForestSRC)

PRED_AGE = 25

for (PRED_AGE in c(20,25,30)) {
  print(paste('PRED_AGE',PRED_AGE, sep = '='))
  pred.data.ready = readRDS(paste('pred.data.ready.age_',PRED_AGE,'.RDs',sep=''))
  pred.data.ready = pred.data.ready[pred.data.ready$timeTillEvent >= 0, ]
  training_ids = readRDS(paste('training_ids.age_',PRED_AGE,'.RDs',sep=''))
  training_idxs = which(rownames(pred.data.ready) %in% as.character(training_ids))
  
  pred.data.ready$isCase = as.integer(pred.data.ready$isCase) - 1
  
  training = pred.data.ready[training_idxs, ]
  testing = pred.data.ready[-training_idxs, ]
  
  rf.fit = rfsrc(Surv(timeTillEvent,isCase)~., data = training, ntree = 128, nodesize = 5, nsplit = 50, importance = TRUE)
  saveRDS(rf.fit, paste('rf.fit.age_',PRED_AGE,'.RDs',sep=''))
}

for (PRED_AGE in c(15, 20,25,30)) {
  print(paste('PRED_AGE',PRED_AGE, sep = '='))
  pred.data.ready = readRDS(paste('pred.data.ready.age_',PRED_AGE,'.RDs',sep=''))
  pred.data.ready = pred.data.ready[pred.data.ready$timeTillEvent >= 0, ]
  pred.data.ready$isCase = as.integer(pred.data.ready$isCase) - 1
  training_ids = readRDS(paste('training_ids.age_',PRED_AGE,'.RDs',sep=''))
  training_idxs = which(rownames(pred.data.ready) %in% as.character(training_ids))
  testing = pred.data.ready[-training_idxs, ]
  
  rf.fit = readRDS(paste('rf.fit.age_',PRED_AGE,'.RDs',sep=''))
  y.pred <- predict(rf.fit,newdata = pred.data.ready)
  
  survival_cases = y.pred$survival[which(y.pred$yvar$isCase == 1), ]
  survival_controls = y.pred$survival[which(y.pred$yvar$isCase == 0), ]
  
  survival_cases_mean = apply(survival_cases, 2, mean)
  survival_cases_sd = apply(survival_cases, 2, sd)
  survival_controls_mean = apply(survival_controls, 2, mean)
  survival_controls_sd = apply(survival_controls, 2, sd)
  
  plot_data = rbind(data.frame(year = 1:length(survival_cases_mean), survival = survival_cases_mean, std = survival_cases_sd, is_case = TRUE), 
                    data.frame(year = 1:length(survival_controls_mean), survival = survival_controls_mean, std = survival_controls_sd, is_case = FALSE))
  
  p <- ggplot(data = plot_data, aes(x = year, y = survival, fill = is_case)) + geom_ribbon(aes(ymin = survival - std, ymax = survival + std, group=is_case), alpha = 0.3) + geom_line(color = "firebrick", size = 1) + theme_minimal() + ggtitle(paste('Pred age = ', PRED_AGE))  
  print(p)
}




























```{r}
require(xgboost)
require(pROC)
require(ggplot2)

pred.data = readRDS('pred.data.ready.RDs')
pred.data$isCase = rownames(pred.data) %in% as.character(cases$subject_num)
training_ids = readRDS('training_ids.RDs')

getXGboostModelPerformance <- function(xgb_model, testing){
  xgb.pred =  predict (xgb_model,testing)
  
  true_outcome = getinfo(testing, 'label')
  xgb.roc.res = roc(true_outcome, xgb.pred)
  return(getSummaryStatsTbl(xgb.roc.res))
}

#
# Return the top 'count' significant variables from a XGBoost model
#
getBestParams <- function(xgb.fit, count = 50){
  mat <- xgb.importance(model = xgb.fit)
  count = min(nrow(mat), count)
  return(unlist(mat[1:count, 'Feature']))
}

convertToXGBoostMatrix <- function(df)
{
  outcome_col = grep('isCase', colnames(df))
  for (col in colnames(df)){ if (class(df[, col]) == 'factor') df[, col] <- as.numeric(df[, col]) - 1 }
  df.xgb = xgb.DMatrix(data = as.matrix(df[,-outcome_col]), label = as.numeric(df$isCase) )
  return(df.xgb)
}

buildAndValidateXGBoost <- function(training, testing, params)
{
  # convert dataset to XGBoost format
  training <- convertToXGBoostMatrix(training)
  testing <- convertToXGBoostMatrix(testing)
  
  # find the nrounds parameter (equivalent to number of trees to grows)
  xgbcv <- xgb.cv( params = params, data = training, nrounds = 100, nfold = 5, showsd = T, stratified = T, print_every_n = 10, early_stopping_rounds = 40, maximize = F)
  best_nrounds = xgbcv$best_iteration; print(paste('selected nrounds:', best_nrounds))
  xgb.fit <- xgb.train (params = params, data = training, nrounds = best_nrounds, watchlist = list(train=training, val=testing), print_every_n = 10, early_stopping_rounds = 10, maximize = F , eval_metric = "error") # optional: you can show simultaneous performance over the testing set using: list(val=testing,train=training)
  
  # validate performance on the testing set
  xgb.pred <- predict (xgb.fit,testing)
  true_outcome = getinfo(testing, 'label')
  xgb.roc.res = roc(true_outcome, xgb.pred)
  print(getSummaryStatsTbl(xgb.roc.res)[2, ])
  
  return(list(training, testing, xgb.fit, xgb.pred))
}

# PRED_AGE = NA; FOLLOW_UP = NA
builModel <- function(pred.data, training_ids, PRED_AGE, FOLLOW_UP)
{
  colnames(pred.data) <- make.names(colnames(pred.data))
  pred.data$ethnicity = factor(as.character(pred.data$ethnicity))
  training_idxs = which(rownames(pred.data) %in% as.character(training_ids))
  training <- pred.data[training_idxs, ]
  testing <- pred.data[-training_idxs, ]
  
  params = list(booster = "gbtree", objective = "binary:logistic", eta=0.3, gamma=0, max_depth=6, min_child_weight=1, subsample=1, colsample_bytree=1)
  
  xgboost.res = buildAndValidateXGBoost(training, testing, params)
  
  return(xgboost.res)
}

xgboost.res = builModel(pred.data, training_ids, NA, NA) # returns: list(training, testing, xgb.fit, xgb.pred)

getXGboostModelPerformance(xgboost.res[[3]], xgboost.res[[2]])

# Show the different variables contributing to the model
xgboost.fit = xgboost.res[[3]] 
mat <- xgb.importance(model = xgboost.fit)
print(xgb.ggplot.importance(importance_matrix = mat[1:20], rel_to_first = F) + ggtitle("") + theme(legend.position = 'none') + theme_minimal())


```


<!-- for (PRED_AGE in c(15,20,25,30)){ -->
    <!--   print(paste('------------------------ Age =',PRED_AGE,'------------------------')) -->
    <!--   for (FOLLOW_UP in c(1,3,5)){ -->
        
        <!--     roc.res = roc(testing$isCase, rf.pred[,2]) -->
          <!--     print(paste('AUC=',round(roc.res$auc, digits = 2))) -->
          <!--     print(getSummaryStatsTbl(roc.res)) -->
          <!--   } -->
    <!-- }   -->
  
  
  
  
  
  
  
  
  
  
  






















# 4. Get models' performance (Table 3)  -------------------------------------------------------------
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(randomForest)
library(pROC)


demographics = readRDS('RPDRml_demographics_sample2.RDs')
all_sujects_ids = demographics$subject_num

cases = readRDS('./cases_broad_sample2.RDs')
controls = readRDS('./controls_sample2.RDs')
# include only true cases and true controls
included_ids = unique(c(cases$subject_num, controls$subject_num))


# postfix = 'AllControls' | 'SelectedControls' --> original cohort, including 80% of the data
# postfix = 'AllControls2' | 'SelectedControls2' --> new cohort including 100% of the data

table3.all_controls <- getModelPerformance2('AllControls2', cases, all_sujects_ids)
table3.selected_controls <- getModelPerformance2('SelectedControls2', cases, included_ids)


table3 = cbind(table3.all_controls, table3.selected_controls)
View(table3[order(table3$Age, table3$TimeWindow), ])

View(table3.all_controls[order(table3.all_controls$TimeWindow), ])
View(table3.selected_controls[order(table3.selected_controls$TimeWindow), ])


res$AUC <- round(res$AUC, digits = 2)
View(res)
write.csv(res, paste('al_models_aucs_',postfix, '.csv') ) # all_ages_stats = read.csv('all_ages_stats.csv')

top_vars = (rownames(imp_vars)[order(imp_vars, decreasing = T)])[1:20]

top_vars_dx = gsub('(dx_)|(0$)','',top_vars[grep('dx_', top_vars)])
top_vars_dx = gsub('ICD9.', 'ICD9:', top_vars_dx)
top_vars_dx_desc = dx_dict$concept_name[match(top_vars_dx, dx_dict$concept_code)]
print(paste(paste(top_vars_dx, top_vars_dx_desc,sep='='), collapse = ', '))

top_vars_procs = gsub('(procs_)|(0$)','',top_vars[grep('procs_', top_vars)])
top_vars_procs = gsub('CPT4.', 'CPT4:', top_vars_procs)
top_vars_procs_desc = procs_dict$concept_name[match(top_vars_procs, procs_dict$concept_code)]
print(paste(paste(top_vars_procs, paste(substr(top_vars_procs_desc,1,200),'...'),sep='='), collapse = ', '))
print(paste(paste(top_vars_procs_desc, paste('(',top_vars_procs,')',sep=''), collapse = ', ')))



#
# Study specific risk-factors
#

model.data.selected_controls = readRDS('modelDataSelectedPatients2.RDs')

all_ids = unique(model.data.selected_controls$subject_num)
cases_ids = unique(model.data.selected_controls$subject_num[model.data.selected_controls$isCase])
controls_ids = unique(model.data.selected_controls$subject_num[!model.data.selected_controls$isCase])

cpt90801 <- model.data.selected_controls[grep('CPT4:90801', model.data.selected_controls$concept_code), ]
cpt90801 <- cpt90801[is.na(cpt90801$time_till_index) | ( !is.na(is.na(cpt90801$time_till_index)) & cpt90801$time_till_index>0 ), ] # keep all controls and all cases with code before FEP
ids_CPT90801 <- unique(cpt90801$subject_num)

showCountandRate(sum(ids_CPT90801 %in% all_ids), length(all_ids))
showCountandRate(sum(ids_CPT90801 %in% cases_ids), length(cases_ids))
showCountandRate(sum(ids_CPT90801 %in% controls_ids), length(controls_ids))


cpt90801 = cpt90801[order(cpt90801$subject_num, cpt90801$date), ]
cpt90801_first = cpt90801[!duplicated(cpt90801$subject_num), ]
summary(cpt90801_first$time_till_index[!is.na(cpt90801_first$time_till_index)])
