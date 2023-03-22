## These functions are based on the ChainLadder package in R. 
## Markus Gesmann, Daniel Murphy, Yanwei (Wayne) Zhang, Alessandro Carrato, Mario Wuthrich, Fabio Concina and Eric Dal Moro (2022). ChainLadder: Statistical
## Methods and Models for Claims Reserving in General Insurance. R package version 0.2.15. https://CRAN.R-project.org/package=ChainLadder

nowcasting_mpox <- function(max_date, # date or integer representing date
                            reporting_triangle_fit, # dataframe
                            reporting_triangle_out, # dataframe
                            reporting_triangle, # dataframe
                            n_samples = 10000, # integer
                            date_list, # vector
                            denom = 14, # integer
                            model = "nonpara", # string
                            target = "specimen" # string
) {
  ## This function is the national nowcast, and can run both nonpara and para, for
  ## specimen date and symptom onset date


  k <- round(length(unique(reporting_triangle_fit$origin)) / denom) # number of knots

  ## Fit the GAM________________________________________________________________
  if (target == "symptom") {
    if (model == "nonpara") {
      gam_fit <- gam(
        value ~ s(as.numeric(origin), k = k, bs = "cr") +
          s(as.numeric(dev), k = 10, bs = "cr") +
          s((dow_dev), bs = "re") +
          s((dow_origin), bs = "re"),
        data = reporting_triangle_fit, family = "nb"
      )
    } else if (model == "para") {
      ## Time delays for symptom date should be shifted as no events were reported on the first two days, to do in revisions
      gam_fit <- gam(
        value ~ s(as.numeric(origin), k = k, bs = "cr") +
          offset(log((pweibull(as.numeric(dev), shape = var1, scale = var2)) - (pweibull(as.numeric(dev - 1), shape = var1, scale = var2)) + 0.00000000000001)) + # without offset
          s((dow_dev), bs = "re") +
          s((dow_origin), bs = "re"),
        data = reporting_triangle_fit, family = "nb"
      )
    } else {
      print("unknown model structure")
    }
  } else if (target == "specimen") {
    if (model == "nonpara") {
      gam_fit <- gam(
        value ~ s(as.numeric(origin), k = k, bs = "cr") +
          s(as.numeric(dev), by = as.factor(weekend_reporting), k = 5, bs = "cr") +
          s((dow_dev), bs = "re") +
          s((dow_origin), bs = "re"),
        data = reporting_triangle_fit, family = "nb"
      )
    } else if (model == "para") {
      ## Time delays for specimen date were shifted as no events were reported on the first day
      gam_fit <- gam(
        value ~ s(as.numeric(origin), k = k, bs = "cr") +
          offset(log((pweibull(as.numeric(dev + 1 - 2), shape = var1, scale = var2)) - (pweibull(as.numeric(dev - 2), shape = var1, scale = var2)) + 0.00000000000001)) + # without offset
          s((dow_dev), bs = "re") +
          s((dow_origin), bs = "re"),
        data = reporting_triangle_fit, family = "nb"
      )
    } else {
      print("unknown model structure")
    }
  } else {
    print("unknown target variable")
  }

  ## Generate posterior sample of gam coefficients______________________________
  mh_sample <- gam.mh(gam_fit, ns = n_samples, burn = 1000, t.df = 40, rw.scale = .25, thin = 1)$bs
  linv <- gam_fit$family$linkinv # inverse link function

  ## Generate model matrix______________________________________________________
  model_matrix <- predict(gam_fit, reporting_triangle_out, type = "lpmatrix")

  ## Preallocate matrices_______________________________________________________
  updated_reports <- matrix(0, length(date_list), n_samples)
  predicted_total_events <- matrix(0, length(date_list) + 1, n_samples)
  updated_reports_total <- matrix(0, 1, n_samples)

  ## Convert reporting triangle to cumulative counts____________________________
  latest <- getLatestCumulative(incr2cum(as.triangle(reporting_triangle)))
  latest <- c(latest, total = sum(latest))

  ## Generate posterior sample of predictions___________________________________
  for (i in 1:n_samples)
  {
    if (model == "nonpara") {
      fits <- linv(model_matrix %*% mh_sample[i, ])
    } else if (model == "para") {
      if (target == "symptom") {
        fits <- linv(model_matrix %*% mh_sample[i, ] + log(((pweibull(as.numeric(reporting_triangle_out$dev + 0.5 - 1), shape = reporting_triangle_out$var1, scale = reporting_triangle_out$var2)) - (pweibull(as.numeric(reporting_triangle_out$dev - 0.5 - 1), shape = reporting_triangle_out$var1, scale = reporting_triangle_out$var2)))))
      } else if (target == "specimen") {
        fits <- linv(model_matrix %*% mh_sample[i, ] + log(((pweibull(as.numeric(reporting_triangle_out$dev + 1 - 2), shape = reporting_triangle_out$var1, scale = reporting_triangle_out$var2)) - (pweibull(as.numeric(reporting_triangle_out$dev - 2), shape = reporting_triangle_out$var1, scale = reporting_triangle_out$var2)))))
      } else {
        print("unknown target variable")
      }
    } else {
      print("unknown model structure")
    }
    fits <- rnbinom(n = length(fits), mu = fits, size = gam_fit$family$getTheta(TRUE)) # generate random prediction sample
    updated_reports[sort(unique(reporting_triangle_out$origin)), i] <- tapply(fits, factor(reporting_triangle_out$origin), sum)
    if (nrow(updated_reports) > length(unique(reporting_triangle_fit$origin))) {
      updated_reports[1:(nrow(updated_reports) - length(unique(reporting_triangle_fit$origin))), i] <- 0
    }
    updated_reports_total[i] <- sum(updated_reports[, i])
    updated_reports_rounded <- round(c(updated_reports[, i], updated_reports_total[i]))
    predicted_total_events[, i] <- latest + updated_reports_rounded
  }

  ## Calculate quantiles and add to data frame__________________________________
  central_pred <- rowQuantiles(predicted_total_events, p = 0.5)
  upper_pred <- rowQuantiles(predicted_total_events, p = 0.975)
  lower_pred <- rowQuantiles(predicted_total_events, p = 0.025)

  predictions <- rowQuantiles(predicted_total_events, p = 0.5)
  predictions_5 <- rowQuantiles(predicted_total_events, p = 0.025)
  predictions_95 <- rowQuantiles(predicted_total_events, p = 0.975)
  predictions_25 <- rowQuantiles(predicted_total_events, p = 0.25)
  predictions_75 <- rowQuantiles(predicted_total_events, p = 0.75)

  date_list <- date_list[1:length(date_list)]
  counts <- as.data.frame(date_list) %>% mutate(n = latest[1:(length(latest) - 1)])
  counts$gam_preds <- predictions[1:(length(predictions) - 1)]
  counts$gam_ci_5 <- predictions_5[1:(length(predictions) - 1)]
  counts$gam_ci_95 <- predictions_95[1:(length(predictions) - 1)]
  counts$gam_ci_25 <- predictions_25[1:(length(predictions) - 1)]
  counts$gam_ci_75 <- predictions_75[1:(length(predictions) - 1)]
  counts$central_pred <- central_pred[1:(length(central_pred) - 1)]
  counts$upper_pred <- upper_pred[1:(length(upper_pred) - 1)]
  counts$lower_pred <- lower_pred[1:(length(lower_pred) - 1)]

  counts <- cbind(counts, predicted_total_events[1:(nrow(predicted_total_events) - 1), ])
  return(list(counts))
}


symptom_nowcast_region <- function(max_date, # date or integer representing date
                                   reporting_triangle_fit, # dataframe
                                   reporting_triangle_out, # dataframe
                                   reporting_triangle, # dataframe
                                   n_samples = 1000, # integer 
                                   date_list_full, # vector
                                   denom = 7 # integer
                                   ) {
  ## this function is the nonpara regional nowcast by symptom onset date

  k <- round(length(unique(reporting_triangle_fit$origin)) / denom)

  ## Fit the GAM________________________________________________________________
  gam_fit <- gam(
    value ~
      s(as.numeric(origin), k = k, bs = "cr") +
      s(as.numeric(origin), as.factor(region), k = k, bs = "fs") +
      s(as.numeric(dev), k = 10, bs = "cr") +
      s((dow_dev), bs = "re") +
      s(as.factor(region), bs = "re") +
      s((dow_origin), bs = "re"),
    data = reporting_triangle_fit, family = "nb"
  )

  ## Generate posterior sample of gam coefficients______________________________
  mh_sample <- gam.mh(gam_fit, ns = 10000, burn = 1000, t.df = 40, rw.scale = .25, thin = 1)$bs
  linv <- gam_fit$family$linkinv

  counts_out <- data.frame()
  region_list <- c("London", "Not-London")
  for (region_ in region_list) {
    ## Filter to region of interest_____________________________________________
    temp <- reporting_triangle_out %>% filter(region == region_)
    tempreporting_triangle <- reporting_triangle %>% filter(region == region_)

    ## Generate model matrix______________________________________________________
    model_matrix <- predict(gam_fit, temp, type = "lpmatrix")

    ## Preallocate matrices_______________________________________________________
    updated_reports <- matrix(0, length(unique(reporting_triangle_out$origin)), n_samples)
    predicted_total_events <- matrix(0, length(unique(reporting_triangle_out$origin)) + 1, n_samples)
    updated_reports_total <- matrix(0, 1, n_samples)

    ## Convert reporting triangle to cumulative counts____________________________
    latest <- getLatestCumulative(incr2cum(as.triangle(tempreporting_triangle)))
    latest <- latest[-(1:(length(latest) - length(unique(reporting_triangle_out$origin))))]
    latest <- c(latest, total = sum(latest))

    ## Generate posterior sample of predictions___________________________________
    for (i in 1:n_samples)
    {
      fits <- linv(model_matrix %*% mh_sample[i, ])
      fits <- rnbinom(n = length(fits), mu = fits, size = gam_fit$family$getTheta(TRUE))
      updated_reports[, i] <- tapply(fits, factor(temp$origin), sum)
      updated_reports_total[i] <- sum(updated_reports[, i])
      updated_reports_rounded <- round(c(updated_reports[, i], updated_reports_total[i]))
      predicted_total_events[, i] <- latest + updated_reports_rounded
    }

    ## Calculate quantiles and add to data frame__________________________________
    predictions <- rowQuantiles(predicted_total_events, p = 0.5)
    predictions_5 <- rowQuantiles(predicted_total_events, p = 0.10)
    predictions_95 <- rowQuantiles(predicted_total_events, p = 0.90)
    predictions_25 <- rowQuantiles(predicted_total_events, p = 0.25)
    predictions_75 <- rowQuantiles(predicted_total_events, p = 0.75)

    date_list <- date_list_full[2:length(date_list_full)]
    counts <- as.data.frame(date_list) %>% mutate(n = latest[1:(length(latest) - 1)])
    counts$gam_preds <- predictions[1:length(predictions) - 1]
    counts$gam_ci_5 <- predictions_5[1:length(predictions) - 1]
    counts$gam_ci_95 <- predictions_95[1:length(predictions) - 1]
    counts$gam_ci_25 <- predictions_25[1:length(predictions) - 1]
    counts$gam_ci_75 <- predictions_75[1:length(predictions) - 1]

    counts <- cbind(counts, predicted_total_events[1:nrow(predicted_total_events) - 1, ])
    counts$region <- region_
    counts_out <- rbind(counts, counts_out)
  }
  return(list(counts_out))
}


nowcast_specimen_region <- function(max_date, # date or integer representing date
                                    reporting_triangle_fit, # dataframe
                                    reporting_triangle_out, #dataframe
                                    reporting_triangle, # dataframr
                                    n_samples = 1000, # integer 
                                    date_list_full, # vector
                                    denom = 7 # integer
                                    ) {
  ## this function is the nonpara regional nowcast by specimen date

  k <- round(length(unique(reporting_triangle_fit$origin)) / denom)

  ## Fit the GAM________________________________________________________________
  gam_fit <- gam(
    value ~
      s(as.numeric(origin), k = k, bs = "cr") +
      s(as.numeric(origin), as.factor(region), k = k, bs = "fs") +
      s(as.numeric(dev), by = as.factor(weekend_reporting), k = 10, bs = "cr") +
      s(as.factor(region), bs = "re") +
      s((dow_dev), bs = "re") +
      s((dow_origin), bs = "re"),
    data = reporting_triangle_fit, family = "nb"
  )

  ## Generate posterior sample of gam coefficients______________________________
  mh_sample <- gam.mh(gam_fit, ns = 10000, burn = 1000, t.df = 40, rw.scale = .25, thin = 1)$bs
  linv <- gam_fit$family$linkinv

  counts_out <- data.frame()
  region_list <- c("London", "Not-London")
  for (region_ in region_list) {
    ## Filter to region of interest_____________________________________________
    temp <- reporting_triangle_out %>% filter(region == region_)
    tempreporting_triangle <- reporting_triangle %>% filter(region == region_)

    ## Generate model matrix______________________________________________________
    model_matrix <- predict(gam_fit, temp, type = "lpmatrix")

    ## Preallocate matrices_______________________________________________________
    updated_reports <- matrix(0, length(unique(reporting_triangle_out$origin)), n_samples)
    predicted_total_events <- matrix(0, length(unique(reporting_triangle_out$origin)) + 1, n_samples)
    updated_reports_total <- matrix(0, 1, n_samples)

    ## Convert reporting triangle to cumulative counts____________________________
    latest <- getLatestCumulative(incr2cum(as.triangle(tempreporting_triangle)))
    latest <- latest[-(1:(length(latest) - length(unique(reporting_triangle_out$origin))))]
    latest <- c(latest, total = sum(latest))

    ## Generate posterior sample of predictions___________________________________
    for (i in 1:n_samples)
    {
      fits <- linv(model_matrix %*% mh_sample[i, ])
      fits <- rnbinom(n = length(fits), mu = fits, size = gam_fit$family$getTheta(TRUE))
      updated_reports[, i] <- tapply(fits, factor(temp$origin), sum)
      updated_reports_total[i] <- sum(updated_reports[, i])
      updated_reports_rounded <- round(c(updated_reports[, i], updated_reports_total[i]))
      predicted_total_events[, i] <- latest + updated_reports_rounded
    }

    ## Calculate quantiles and add to data frame__________________________________
    predictions <- rowQuantiles(predicted_total_events, p = 0.5)
    predictions_5 <- rowQuantiles(predicted_total_events, p = 0.10)
    predictions_95 <- rowQuantiles(predicted_total_events, p = 0.90)
    predictions_25 <- rowQuantiles(predicted_total_events, p = 0.25)
    predictions_75 <- rowQuantiles(predicted_total_events, p = 0.75)

    date_list <- date_list_full[2:length(date_list_full)]
    counts <- as.data.frame(date_list) %>% mutate(n = latest[1:(length(latest) - 1)])
    counts$gam_preds <- predictions[1:length(predictions) - 1]
    counts$gam_ci_5 <- predictions_5[1:length(predictions) - 1]
    counts$gam_ci_95 <- predictions_95[1:length(predictions) - 1]
    counts$gam_ci_25 <- predictions_25[1:length(predictions) - 1]
    counts$gam_ci_75 <- predictions_75[1:length(predictions) - 1]

    counts <- cbind(counts, predicted_total_events[1:nrow(predicted_total_events) - 1, ])
    counts$region <- region_
    counts_out <- rbind(counts, counts_out)
  }
  return(list(counts_out))
}
