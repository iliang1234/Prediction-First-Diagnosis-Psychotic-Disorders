# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                         Functions Decleration
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

reduceNumberOfFactorLevels <- function(fctr, num_levels){
  s = summary(fctr)
  lvls_to_include = names(s)[1:num_levels]
  fctr = as.character(fctr)
  fctr[!(fctr %in% lvls_to_include)] <- 'Other'
  fctr = factor(fctr)
  
  return(fctr)
}

getSummaryStatsTbl <- function(roc.res){
  res = data.frame()
  for(spec in c(0.8,0.9,0.95,0.99)){
    res = rbind(res, round(coords(roc.res, x = spec, input = 'specificity', ret = c('spec','sen','ppv','npv', 'tp','fp','tn','fn'), transpose = F), digits = 2))
  }
  res$AUC = round(roc.res$auc, digits = 2)
  for(col in c('tp','fp','tn','fn')) res[,col] <- round(res[,col], digits = 0)
  return(res)
}

prepareModelData <- function(df, demographics, min_count = 20)
{
  # get the ids of all subjects initially in dataset (before we do any filtering)
  all_ids = unique(df$subject_num)
  
  # get list of concepts that appear for more than min_count subjects
  concept_count = tapply(df$subject_num, df$concept_code, function(x) sum(!duplicated(x)))
  selected_concepts = unique(names(concept_count)[concept_count > min_count])
  
  df = df[df$concept_code %in% selected_concepts, ]
  df$concept_type_code = paste(df$concept_type, df$concept_code,sep='_')
  df = df[, c('subject_num', 'concept_type_code')]
  
  # to make sure that all patients appear in final table, create empty (stub) rows for patients with missing data
  missing_ids = all_ids[!(all_ids %in% df$subject_num)]
  if(length(missing_ids) > 0){
    df = rbind(df, data.frame(subject_num = missing_ids, concept_type_code = 'stub_stub'))  
  }
  
  tbl = as.data.frame.matrix(table(df$subject_num, df$concept_type_code))
  # remove the 'stub' column
  if(length(grep('stub_stub', colnames(tbl))) > 0){
    tbl = tbl[, -grep('stub_stub', colnames(tbl))]  
  }
  
  # add demographics to the model data 
  tbl$SEX = demographics$gender[match(rownames(tbl), as.character(demographics$subject_num))]
  tbl$MARITAL_STATUS = demographics$marital_status[match(rownames(tbl), as.character(demographics$subject_num))]
  tbl$race = factor(as.character(demographics$race[match(rownames(tbl), as.character(demographics$subject_num))]))
  tbl$ethnicity = factor(as.character(demographics$ethnicity[match(rownames(tbl), as.character(demographics$subject_num))]))
  tbl$veteran = factor(as.character(demographics$veteran[match(rownames(tbl), as.character(demographics$subject_num))]))
  tbl$public_payer = factor(as.character(demographics$public_payer[match(rownames(tbl), as.character(demographics$subject_num))]))
  
  return(tbl)
}


addDescriptionToPredictors <- function(predictors){
  dx_idxs = grep('dx_', predictors)
  procs_idxs = grep('procs_', predictors)
  labs_idxs = grep('labs_', predictors)
  meds_idxs = grep('meds_', predictors)
  
  res = data.frame(pred = predictors)
  
  res$description = ''
  
  if (length(dx_idxs) > 0){
    res$description[dx_idxs] = dx_dict$concept_name[match(gsub('dx_','',predictors[dx_idxs]), make.names(dx_dict$concept_code))]
  }
  if (length(procs_idxs) > 0){
    res$description[procs_idxs] = procs_dict$concept_name[match(gsub('procs_','',predictors[procs_idxs]), make.names(procs_dict$concept_code))]
  }
  
  if (length(labs_idxs) > 0){
    res$description[labs_idxs] = labs_dict$concept_name[match(gsub('labs_','',predictors[labs_idxs]), make.names(labs_dict$concept_code))]
  }
  if (length(meds_idxs) > 0){
    res$description[meds_idxs] = meds_dict$concept_name[match(gsub('meds_','',predictors[meds_idxs]), make.names(meds_dict$concept_code))]
  }
  rownames(res) = 1:nrow(res)
  return(res)
}

print_count_and_percentage <- function(count, total) {
  percentage <- (count / total) * 100
  message <- paste("n=", count, " (", round(percentage, 2), "%)", sep = "")
  return(message)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                         Code Starts Here
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#
# Load required libraries
# ----------------------------------------------------------------------------
.libPaths("~/PNGU/R_PACKAGES/")

# install.packages('ggplot2', dependencies = TRUE)
install.packages('randomForest', dependencies = TRUE)
# remotes::install_github("NorskRegnesentral/shapr", dependencies = TRUE)

# install.packages('DALEX')
# install.packages('ggplot2')
# install.packages("shapper")
# install.packages('fastshap')
install.packages('shapr', force = TRUE)
install.packages('treeshap')
install.packages('shapviz')
install.packages('vip')

install.packages("DALEX")
shapper::install_shap()

install.packages("kernelshap", dependencies = TRUE)

# install.packages('caret', dependencies = TRUE)
# install.packages('dplyr', dependencies = TRUE)
# install.packages('pROC', dependencies = TRUE)

library(ggplot2)
library(caret)
library(randomForest, lib = "/data/js95/PNGU/R_PACKAGES")
library(dplyr, lib = "/data/js95/PNGU/R_PACKAGES")
library(pROC)
library(shapr)
library(DALEX)
library(ggplot2)
library(shapper)
library(fastshap)
library(treeshap)
library(kernelshap, lib = "/tmp/Rtmp6wTf3q/downloaded_packages")
library(vip)

#
# Load data
# ----------------------------------------------------------------------------
setwd("/data/js95/PNGU/RPDRml_v2021.1/sample2/")
demographics <- readRDS('/data/js95/PNGU/RPDRml_v2021.1/sample2/RPDRml_demographics_sample2.RDs')
cases = readRDS('./cases.RDs')
controls = readRDS('controls_sample2.RDs') 
# the list of index dates (index event for cases, last documented encounter for controls)
index_events = readRDS('index_events.RDs')

# dictionaries to name top predictors
dx_dict = readRDS('../RPDRml__DX_dictionary.RDs')
procs_dict = read.delim('../RPDRml__PROC_dictionary.txt', stringsAsFactors = F)
labs_dict = read.delim('../RPDRml__LAB_dictionary.txt', stringsAsFactors = F)
meds_dict = read.delim('../RPDRml__MED_dictionary.txt', stringsAsFactors = F)



#
# Convert to flat (wide) format
# ----------------------------------------------------------------------------

print('Create wide/flat version of prediction dataset')
for(PRED_AGE in c(15, 20, 25,30)){
  print(paste0('PRED_AGE = ', PRED_AGE))
  pred.data = readRDS(paste('./model_datasets/pred.data.age_',PRED_AGE,'.RDs',sep=''))
  
  print(paste('Number of patients=', sum(!duplicated(pred.data$subject_num))))
  print(paste('Cases', print_count_and_percentage(sum( unique(pred.data$subject_num) %in% cases$subject_num), sum(!duplicated(pred.data$subject_num)))) )
  
  # the dataset is too large so we will first down-sample it for this part
  controls_sample = unique(sample(controls$subject_num, length(cases$subject_num) * 2, replace = F))
  pred.data.sample = pred.data[pred.data$subject_num %in% c(cases$subject_num, controls_sample), ]
  
  pred.data.flat = prepareModelData(pred.data.sample, demographics, min_count = 20)
  print(paste('Number of features:', ncol(pred.data.flat)))
  
  # add case definition  
  pred.data.flat$isCase <- factor(ifelse(rownames(pred.data.flat) %in% as.character(cases$subject_num), 1, 0))
  
  # add trainin/testing flag
  training_ids = unique(pred.data$subject_num[pred.data$isTraining == TRUE])
  pred.data.flat$isTraining <- ifelse(rownames(pred.data.flat) %in% as.character(training_ids), 1, 0)
  
  # pred.data.flat$timeTillEvent = index_events$age_at_index[match(rownames(pred.data.flat), as.character(index_events$subject_num))] - PRED_AGE
  # saveRDS(pred.data.flat, paste('./model_datasets/pred.data.flat.age_',PRED_AGE,'.RDs',sep=''))  
  saveRDS(pred.data.flat, paste('~/Psychosis_analysis/No_Race_Data/pred.data.flat.age_',PRED_AGE,'.RDs',sep=''))  
}


# pred.data.flat = readRDS(paste('./model_datasets/pred.data.flat.age_',PRED_AGE,'.RDs',sep=''))  
# colnames(pred.data.flat) <- make.names(colnames(pred.data.flat))
# pred.data.flat = pred.data.flat[, !duplicated(colnames(pred.data.flat))]


#
# Run feature selection
# ----------------------------------------------------------------------------

print('Run feature selection')
for(PRED_AGE in c(20, 25,30)){
  print(paste0('PRED_AGE = ', PRED_AGE))
  # read the flat/wide dataset
  pred.data.flat = readRDS(paste('./model_datasets/pred.data.flat.age_',PRED_AGE,'.RDs',sep=''))
  # extract training set out of it
  training = pred.data.flat[pred.data.flat$isTraining == 1, ]
  colnames(training) <- make.names(colnames(training))
  training = training[, !duplicated(colnames(training))]
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
  levels(training$isCase) <- c('Control','Case')
  rf_model <- train(isCase ~ ., 
                    data = training, 
                    method = "rf", 
                    metric = "ROC",
                    tuneGrid = grid, 
                    trControl = train_control)
  
  # Get feature importance
  importance <- varImp(rf_model, scale = FALSE)
  # Get the names of the top 500 variables
  topvars = rownames(importance$importance)[order(importance$importance, decreasing = T)[1:500]]
  # remove 'cheat codes' (currently it's only the timeTillEvent code)
  topvars = c(topvars[!(topvars %in% 'timeTillEvent')], 'isCase')
  
  # for categorical values - we don't want to include specific values, but the entire variable:
  # This will give you the names of the factor columns
  factor_columns_names <- names(training)[sapply(training, is.factor)]
  
  # For each categorical variable, convert from var+value to just var: 
  topvars <- sub("^SEX.*", "SEX", topvars)
  topvars <- sub("^MARITAL_STATUS.*", "MARITAL_STATUS", topvars)
  topvars <- sub("^race.*", "race", topvars)
  topvars <- sub("^ethnicity.*", "ethnicity", topvars)
  topvars <- sub("^veteran.*", "veteran", topvars)
  topvars <- sub("^public_payer.*", "public_payer", topvars)
  topvars = unique(topvars)
  
  topvars <- make.names(gsub("`", "", topvars))
  
  # look for possible cheat codes
  topvars_description = addDescriptionToPredictors(topvars)
  cheat_codes = grep('psychosis', topvars_description$description, ignore.case = T)
  print(paste('number  of potential cheat codes:', length(cheat_codes)))
  if(length(cheat_codes) > 0){
    print(topvars_description[cheat_codes, ])
    # remove these possible cheat codes from the model
    topvars = topvars[-cheat_codes]
  }
  
  saveRDS(topvars, paste('~/Psychosis_analysis/No_Race_Data/top_vars.age_',PRED_AGE,'.RDs',sep=''))
}


# topvars = readRDS(paste('./model_datasets/top_vars.age_',PRED_AGE,'.RDs',sep=''))


#
# Fine-tune hyper parameters
# ----------------------------------------------------------------------------
FIND_HYPER_PARAMS = FALSE
# For computational purposes, only do hyper param tuning once, using the age=20 group

if (FIND_HYPER_PARAMS) {
  PRED_AGE = 20
  pred.data.flat = readRDS(paste('./model_datasets/pred.data.flat.age_',PRED_AGE,'.RDs',sep=''))
  topvars = readRDS(paste('./model_datasets/top_vars.age_',PRED_AGE,'.RDs',sep=''))
  # extract training set out of it
  training = pred.data.flat[pred.data.flat$isTraining == 1, ]
  colnames(training) <- make.names(colnames(training))
  training = training[, topvars]
  
  # 2) hyperparameter tuning using k-fold cross validation
  hyper_grid <- expand.grid(
    ntree = c(128, 256, 512, 1024),
    mtry = seq(2, 10, 2),
    nodesize = seq(1, 5, 1),
    nsplit = c(1, 3, 5, 10, 15)
  )
  
  k <- 3
  folds <- createFolds(training$isCase, k = k, list = TRUE, returnTrain = TRUE)
  
  results <- matrix(NA, nrow = nrow(hyper_grid), ncol = k)
  colnames(results) <- paste0("Fold", 1:k)
  
  for (i in 1:nrow(hyper_grid)) {
    print(paste('hyper param: row',i,'/',nrow(hyper_grid)))
    ntree <- hyper_grid$ntree[i]
    mtry <- hyper_grid$mtry[i]
    nodesize <- hyper_grid$nodesize[i]
    nsplit <- hyper_grid$nsplit[i]
    
    for (j in 1:k) {
      training_set <- training[folds[[j]], ]
      validation_set <- training[-folds[[j]], ]
      
      model <- randomForest(isCase ~ . ,
                            data = training_set,
                            ntree = ntree,
                            mtry = mtry,
                            nodesize = nodesize,
                            nsplit = nsplit
      )
      
      predictions <- predict(model, validation_set, type = 'prob') 
      
      roc.res = roc(validation_set$isCase, predictions[,'Case'])
      
      results[i, j] <- roc.res$auc
    }
    print(sprintf("ntree=%i, mtry=%i, nodesize=%i, nsplit=%i --> AUC=%.2f", ntree, mtry, nodesize, nsplit, mean(results[i,], na.rm = T)))
    
  }
  
  # determine the best configuration
  
  mean_auc <- rowMeans(results, na.rm = TRUE)
  
  hyper_grid$mean_auc <- mean_auc
  max_auc =  max(hyper_grid$mean_auc)
  print(max_auc)
  
  best_hyperparameters <- hyper_grid[which.max(hyper_grid$mean_auc), ]
  print(best_hyperparameters)
  saveRDS(best_hyperparameters, paste('~/Psychosis_analysis/No_Race_Data/hyper_parameters.age_',PRED_AGE,'.RDs',sep=''))  
}


# 
# execute model
# ----------------------------------------------------------------------------
best_hyperparameters = readRDS('./model_datasets/hyper_parameters.age_20.RDs')
optimal_ntree <- best_hyperparameters$ntree  # 512
optimal_mtry <- best_hyperparameters$mtry # 10
optimal_nodesize <- best_hyperparameters$nodesize # 5 (2)
optimal_nsplit <- best_hyperparameters$nsplit # 10 (15)

EXCLUDE_MCL = F
if (EXCLUDE_MCL){
  # MGB - filter patients all patients who had index event at MCL before 2018. Compare results between models. 
  encounters = readRDS('./encounters2.RDs') # derived from diagnosis dataset
  MCL_before_2018 = encounters[encounters$site == 'MCL' & encounters$sstart_date < as.Date('2018-01-01'), ]
  ids_to_exclude = unique(MCL_before_2018$subject_num)  
} else {
  ids_to_exclude = c()
}

print('Train and validate the model')

RETRAIN_MODEL = FALSE
GET_PREDICTORS = TRUE
CALC_SUMMARY_STATS = FALSE

age_range = readRDS('range_of_encounters.RDs')

age_range$isCase = age_range$subject_num %in% cases$subject_num
age_range$indexEvent = cases$index_age[match(age_range$subject_num, cases$subject_num)]

for(PRED_AGE in c(15, 20, 25,30))
{
  print(paste0('PRED_AGE = ', PRED_AGE))
  # read the long-format dataset (this dataset already excludes all subjects who met the case definition before PRED_AGE)
  pred.data = readRDS(paste('./model_datasets/pred.data.age_',PRED_AGE,'.RDs',sep=''))
  
  # filter the dataset to only include selected features
  topvars = readRDS(paste('./model_datasets/top_vars.age_',PRED_AGE,'.RDs',sep=''))
  pred.concepts = make.names(paste(pred.data$concept_type, pred.data$concept_code, sep = '_'))
  pred.data = pred.data[pred.concepts %in% topvars, ]
  
  # convert to long/flat format
  pred.data.flat = prepareModelData(pred.data, demographics, min_count = 20)
  colnames(pred.data.flat) <- make.names(colnames(pred.data.flat))
  
  # Optional: only use a time-window of 5 years for prediction
  potential_subjects = age_range$subject_num[age_range$first_age < PRED_AGE & age_range$last_age > (PRED_AGE+5)]
  cases_to_include = cases$subject_num[!is.na(cases$index_age) & (cases$index_age > PRED_AGE) & (cases$index_age < (PRED_AGE+5))]
  potential_subjects = as.character(unique(c(potential_subjects, cases_to_include)))
  
  pred.data.flat = pred.data.flat[rownames(pred.data.flat) %in% potential_subjects, ]
  pred.data.flat$isCase <- factor(ifelse(rownames(pred.data.flat) %in% as.character(cases_to_include), 1, 0))
  
  # extract training & testing sets out of it
  training_ids = unique(pred.data$subject_num[pred.data$isTraining == TRUE])
  pred.data.flat$isTraining <- ifelse(rownames(pred.data.flat) %in% as.character(training_ids), 1, 0)
  training = pred.data.flat[pred.data.flat$isTraining == 1, ]
  testing = pred.data.flat[pred.data.flat$isTraining == 0, ]
  
  if (RETRAIN_MODEL){
    print('train model using training set...')
    # first, train the tuned Random Forest model using the training data
    # rf.fit = randomForest(isCase ~ . - race, data = training[!(rownames(training) %in% as.character(ids_to_exclude)), ], 
                          # ntree = optimal_ntree, 
                          # mtry = optimal_mtry, 
                          # nodesize = optimal_nodesize, 
                          # nsplit = optimal_nsplit)
    rf.fit = randomForest(isCase ~ ., data = training[!(rownames(training) %in% as.character(ids_to_exclude)), ], 
                          ntree = optimal_ntree, 
                          mtry = optimal_mtry, 
                          nodesize = optimal_nodesize, 
                          nsplit = optimal_nsplit)
    # saveRDS(rf.fit, paste('rf.fit.window5y.age_',PRED_AGE,'.RDs',sep=''))    
  } else {
    rf.fit = readRDS(paste('rf.fit.window5y.age_',PRED_AGE,'.RDs',sep=''))
  }

  print('validate model on testing set...')
  # apply the fitted/trained model on the testing set
  y.pred <- predict(rf.fit,newdata = testing, type = 'prob')
  
  demographics_with_pred = index_events[match(rownames(testing), index_events$subject_num), ]
  demographics_with_pred = cbind(demographics_with_pred, y.pred)
  
  # saveRDS(demographics_with_pred, paste0('~/Psychosis_analysis/No_Race_Data/rf_model_pred_res_age_', PRED_AGE, '.RDs'))
  
  if (CALC_SUMMARY_STATS){
    # Calculate ROC and summary statistics
    roc.res = roc(factor(testing$isCase), y.pred[,2])
    summaryStats = getSummaryStatsTbl(roc.res)
    # write.csv(summaryStats, paste0('~/Psychosis_analysis/No_Race_Data/model_results_window5y_age',PRED_AGE,'.csv'))
    print(summaryStats)
  }
  
  library(shapviz)
  library(kernelshap)
  library(randomForest)
  GET_PREDICTORS = TRUE
  #
  # Get top predictors
  # ----------------------------------------------------------------------------
  if (GET_PREDICTORS){
    # print('Get top predictors...')
    # var_imp = importance(rf.fit)
    # NUM_PREDS = 20
    # var_imp = sort(var_imp[, "MeanDecreaseGini"], decreasing = T)
    # var_imp = var_imp[1:NUM_PREDS]
    # top_predictors = addDescriptionToPredictors(names(var_imp))
    # top_predictors$Gini = var_imp
    
    # write.csv(top_predictors, paste0('~/Psychosis_analysis/top_predictors_window5y_age',PRED_AGE,'.csv'))
    # View(top_predictors) 
    
    # NUM_PREDS = 20
    # predict_function <- function(rf, testing) predict(rf, testing, type = "prob")
    
    # Calculate SHAP values for the test set
    # shap_values <- fastshap::explain(rf, X = testing[, -5], pred_wrapper = predict_function, nsim = 100)
    
    # Calculate mean absolute SHAP values for each feature
    # shap_importance <- colMeans(abs(shap_values))
    
    # Create a data frame for feature importances
    # feature_importance_df <- data.frame(
      # Feature = names(shap_importance),
      # SHAP_Importance = shap_importance
    # )
    # Sort features by importance
    # feature_importance_df <- feature_importance_df %>% arrange(desc(SHAP_Importance))
    
    # Print the top features
    # print(feature_importance_df)
    
    # Prepare the data for explanation
    training_data <- training[!(rownames(training) %in% as.character(ids_to_exclude)), ]
    
    x_train <- training[, names(training)!='isCase']
    y_train <- training[, 'isCase']
    xvars <- row.names(x_train)
    
    # 1) Sample rows to be explained (recommendation is 500 to 2k rows)
    set.seed(10)
    # x <- x_train[sample(nrow(x_train), 1000), xvars]
    
    # 2) Select background data (recommendation is 50 to 500 rows)
    bg_X <- x_train[sample(nrow(x_train), 200), ]
    
    pred_rf <- function(rf_model, df) {
      predict(rf_model, newdata = df, type = 'vote')[,2]
    }
    
    # 3) Kernel SHAP 
    selected_features = vi(rf.fit)$Variable[1:200]
    x_train_subset <- x_train[sample(nrow(x_train), 500), ]
    system.time( 
      ks <- kernelshap(rf.fit, x_train_subset, bg_X = bg_X, pred_fun = pred_rf, feature_names = selected_features)
    )
    ks
    
    # 4) Analyze with our sister package {shapviz}
    viz <- shapviz(ks)
    sv_importance(viz)
    
    save(ks, file = shaps.RData)
  }
}