library(dplyr)
library(ChainLadder)
library(mgcv)
library(matrixStats)

## Ensure working directory is the git root directory
source("./src/nowcast_functions.R")

## Model specification parameters_________________________________________________
denom <- 14 # average number of days per knot
model <- "nonpara" # nonpara or para
outcome <- "symptom" # symptom or specimen

if (outcome == "symptom"){
  model_run_range <- 105:160 # date range for symptom onset data
}
if (outcome == "specimen"){
  model_run_range <- 75:160 # date range for specimen date data
}

## Define lead times we are interested in_______________________________________
lead_times <- c(1,2,3,4,5,6,7,8,9)
output_scoring <- list()
for (j in lead_times){
  output_scoring[[j]] <- data.frame()
}

for (d in model_run_range){
  print(d)
  max_date = d
  
  ## Nowcast data pre-processing__________________________________________________
  lda <- read.csv(paste0("C:/Users/christopher.overton/OneDrive - UK Health Security Agency/Documents/Projects/Monkeypox/publication_data/data/data/nowcasting_input_",outcome,"_",d,".csv"))
  date_list <- unique(lda$origin)
  lda <- lda %>% mutate(origin = origin - min(lda$origin) + 1)
  lda$dow_dev <- as.factor(lda$dow_dev)
  lda$dow_origin <- as.factor(lda$dow_origin)
  counts_observed <- read.csv(paste0("C:/Users/christopher.overton/OneDrive - UK Health Security Agency/Documents/Projects/Monkeypox/publication_data/data/data/full_data_",outcome,".csv"))

  ldaFit <- subset(lda, !is.na(lda$value))
  ldaOut <- subset(lda, is.na(lda$value))
  ldaFit$value <- round(ldaFit$value)

  ###_____________________________________________________________________________
  
  ## Nowcast calculation__________________________________________________________
  output <- try(nowcasting_mpox(max_date = max_date,ldaFit = ldaFit,ldaOut = ldaOut,lda = lda,n_samples = 10000,date_list = date_list,denom = denom,model = model,target = outcome))
  counts <- try(output[[1]])

  ### Append earlier observations (before model fitting window) for fitting growth rates to the full time series
  counts_old <- try(counts_observed %>% filter(origin < min(counts$date_list)) %>% rename("date_list"="origin"))
  try(counts_old[,names(counts %>% dplyr::select(-n,-date_list))] <- counts_old$n)
  counts <- try(bind_rows(counts_old,counts))
  
  ## If the model failed, add NA________________________________________________
  if (class(output) == "try-error"){
    counts <- as.data.frame(date_list)
    counts <- counts %>% mutate(
      gam_preds = NA,
      gam_ci_5 = NA,
      gam_ci_25 = NA,
      gam_ci_75 = NA,
      gam_ci_95 = NA,
      central_pred = NA,
      lower_pred = NA,
      upper_pred = NA,
      n = NA
    )
  }
  
  ## Output results for different lead times of interest________________________
  for (j in lead_times){
    output_scoring[[j]] <- rbind(output_scoring[[j]],
                                 counts[(nrow(counts)-j+1),c("date_list","gam_preds","gam_ci_5","gam_ci_25","gam_ci_75","gam_ci_95",
                                                             "central_pred","lower_pred","upper_pred","n")])
  }
}

## Combined results across all model runs_______________________________________
combined = bind_rows(output_scoring, .id = "id") %>% dplyr::select(-n) %>%
  left_join(counts_observed %>% dplyr::select(origin,n), by = c("date_list" = "origin") ) %>%
  mutate(res = n -gam_preds,rel_res = (n-gam_preds)/n) %>%
  mutate(id = paste0("Lead time = ",as.numeric(id)-1))


#################
# OUTPUT
################

output_path <- paste0("./outputs/",
                      "model_scoring_",outcome,"_",max_date,"_",model,"_",denom,"no05.csv")

write.csv(combined, file = output_path)

