# nowcasting_the_2022_mpox_outbreak_in_england
Code and data associated with the paper "Nowcasting the 2022 mpox outbreak in England"

arXiv: https://arxiv.org/abs/2302.09076

Code for running parametric, non-parametric, and regional nowcasts are provided. 

A subset of the data used in the paper is provided. Individual level data could not be shared, but we have provided the estimates for the shape and scale parameters of the right-truncation corrected Weibull distribution fitted to the individual level data, from which the parametric model can be used. The modified versions of the data are provided, since these showed the strongest model performance. Dates have been removed, and data are only provided for a subset of the data range considered in the paper. 

Note that since submission of the paper, improvements have been identified in the model structure. These improvements will be added to the code/data here as they are implemented during revisions to the paper, so this repo will change over time. 

# Structure of the repo.

## /publication_data/
This folder contains all the data needed to run the nowcasting models, which can easily be adapted to meet the formats of other nowcasting models. There are 4 types of data file included: nowcasting_input_specimen_(date), nowcasting_input_symptom_(date), full_data_specimen, and full_data_symptom. The nowcasting input files are the inputs required for the nowcasting models, representing the state of the reporting triangle as of (date). The structure of these files follows the structure of the file inputs used in the ChainLadder R package. The full_data files contain the full data as they looked 2 months after the final date, representing the actual observations once all backfilling had taken place. These should be considered the ground truth when assessing model performance. 

The "nowcasting_input" files have the following structure:
### origin (numeric variable to represent a date)
The date the event occurred.
### dev (numeric variable)
The delay between the event occurring and being reported. Note that a value of 1 represents occurring on the day of the event (so a delay between 0 and 1 days) and so on. 
### value (numeric variable)
The number of events according to that origin and dev pair.
### dow_origin (factor variable)
Day of week for the origin.
### dow_dev (factor variable)
Day of week for the dev.
### weekend_reporting (factor variable)
Whether the origin date was affected by weekend processing changes, as described in the paper.
### var1 (numeric variable)
The shape parameter of the Weibull distribution fitted to the individual level data. Note that for specimen date, the Weibull distribution is fitted to the observed delay -1, since no events were reported on the day of occurrence. For symptom onset, this was not the case, but with the 2 day modification we will make this change during revisions as it should improve the model performance. 
### var2 (numeric variable)
The scale parameter of the Weibull distribution fitted to the individual level data. Note that for specimen date, the Weibull distribution is fitted to the observed delay -1, since no events were reported on the day of occurrence. For symptom onset, this was not the case, but with the 2 day modification we will make this change during revisions as it should improve the model performance. 

The "full_data" files have the following structure:
### origin (numeric variable to represent a date)
The date the event occurred.
### n (numeric variable)
The total number of events that occurred on that day. 

## /src/
This folder contains the R script "nowcast_functions.R" which contains the functions used for running the different nowcasting models. Note that code for nonparametric regional models is included, but no data with a regional breakdown is included, so that functions cannot be used. 

## /scripts/
This folder contains the R script "call_nowcasting_models" which loads the required input data file and calls the corresponding nowcasting function. This function loops over the full range of available dates for each scenario. 

# Dependencies
These functions require other R packages to run. The dependencies, and versions used in the paper are,

dplyr - version 1.0.9

ChainLadder - version 0.2.15

mgcv - version 1.8.40 (1.8.38 or higher required)

matrixStats - version 0.62.0
