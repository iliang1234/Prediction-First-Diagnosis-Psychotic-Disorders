# Prediction of First-Diagnosis Psychotic Disorders

This repository contains the analysis pipeline for predicting first-diagnosis psychotic disorders using R-based statistical modeling.

## Project Overview

This project aims to develop predictive models for identifying individuals at risk of receiving their first diagnosis of psychotic disorders. The analysis includes model development using data from two healthcare systems: Boston Children's Hospital (BCH) and Massachusetts General Hospital (MGH). The analysis is structured in three main steps:

1. Data Import and Case Definition
2. Data Preparation for Modeling
3. Model Development and Evaluation

## Repository Structure

### `BCH_model_training/`
- `Psychosis Step 1 - import files and apply case def.Rmd`: Initial data processing and case definition application
- `Psychosis Step 2 - Prepare data for modelling v2.Rmd`: Feature engineering and data preparation for modeling
- `Psychosis Step 3 - Run Models v2.Rmd`: Implementation of predictive models and analysis

### `MGH_model_training/`
- `Psychosis Step 1 - import data.R`: Initial data processing for MGH dataset
- `Psychosis Step 2 - prepare data.R`: Feature engineering and data preparation
- `Psychosis Step 3 - run analysis.R`: Model implementation and analysis
- `Train_Without_Race.R`: Additional model training excluding race variables

### `model_result_analysis/`
Directory containing model evaluation results and analysis outputs at MGH 

## Getting Started

### Prerequisites

- R (recommended version >= 4.0)
- RStudio (for working with .Rmd files)
- Python version >= 3.0 (analysis code ran with Python 3.6.15)
- Jupyter Notebook (for working with the analysis notebooks)

### Workflow

1. Execute the R Markdown files in sequential order:
   - Start with Step 1 for data import and case definitions
   - Proceed to Step 2 for data preparation
   - Finally, run Step 3 for model development and evaluation

## Analysis Pipeline

Each step in the analysis pipeline is documented in detail within the respective R Markdown files, ensuring reproducibility and transparency of the research process.

## Results

Model results and detailed analysis can be found in the `model_result_analysis` directory.


## Contact

For questions, please contact Ivy Liang (Ivy.Liang@childrens.harvard.edu) and Ben Reis (Ben.Reis@childrens.harvard.edu).