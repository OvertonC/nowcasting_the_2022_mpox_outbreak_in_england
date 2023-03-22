nowcasting_mpox <- function(max_date,
                            df_dates,
                            ldaFit,
                            ldaOut,
                            lda,
                            n_samples = 10000, 
                            date_list,
                            denom = 14,
                            model = "nonpara",
                            target = "specimen"
                            ){
## This function is the national nowcast, and can run both nonpara and para, for
## specimen date and symptom onset date


  k = round(length(unique(ldaFit$origin))/denom) # number of knots
  
  ## Fit the GAM________________________________________________________________
  if (target == "symptom"){
    if (model == "nonpara"){
      gamFit <- gam(value ~ s(as.numeric(origin),k = k ,bs = "cr") + 
                    s(as.numeric(dev), k = 10, bs = "cr") +
                    s((dow_dev),bs = "re") +
                    s((dow_origin),bs = "re"),
                  data = ldaFit, family = "nb")
    } else if (model == "para"){
      ## Time delays for symptom date should be shifted as no events were reported on the first two days, to do in revisions
      gamFit <- gam(value ~ s(as.numeric(origin),k= k, bs = "cr") + 
                    offset(log((pweibull(as.numeric(dev),shape = var1, scale = var2))-(pweibull(as.numeric(dev-1),shape = var1, scale = var2))+0.00000000000001)) + # without offset
                    s((dow_dev),bs = "re") +
                    s((dow_origin),bs = "re"),
                  data = ldaFit, family = "nb")
    } else {
      print("unknown model structure")
    }
  } else if (target == "specimen"){
    if (model == "nonpara"){
      gamFit <- gam(value ~ s(as.numeric(origin),k = k ,bs = "cr") + 
                      s(as.numeric(dev), by = as.factor(weekend_reporting), k = 5, bs = "cr") +
                      s((dow_dev),bs = "re") +
                      s((dow_origin),bs = "re"),
                    data = ldaFit, family = "nb")
    } else if (model == "para"){
      ## Time delays for specimen date were shifted as no events were reported on the first day 
      gamFit <- gam(value ~ s(as.numeric(origin),k= k, bs = "cr") + 
                      offset(log((pweibull(as.numeric(dev+1-2),shape = var1, scale = var2))-(pweibull(as.numeric(dev-2),shape = var1, scale = var2))+0.00000000000001)) + # without offset
                      # offset(log((pweibull(as.numeric(dev+0.5-2),shape = var1, scale = var2))-(pweibull(as.numeric(dev-0.5-2),shape = var1, scale = var2))+0.00000000000001)) + # without offset
                      s((dow_dev),bs = "re") +
                      s((dow_origin),bs = "re"),
                    data = ldaFit, family = "nb")
    } else {
      print("unknown model structure")
    }  
  } else {
    print("unknown target variable")
  }
  
  ## Generate posterior sample of gam coefficients______________________________
  mh_sample <- gam.mh(gamFit,ns=n_samples,burn=1000,t.df=40,rw.scale=.25,thin=1)$bs
  linv <- gamFit$family$linkinv # inverse link function
  
  ## Generate model matrix______________________________________________________
  Xp <- predict(gamFit, ldaOut, type = "lpmatrix")
  
  ## Preallocate matrices_______________________________________________________
  resMeanAy_out <- matrix(0,length(date_list),n_samples)
  resMeanAy_pred <- matrix(0,length(date_list),n_samples)
  Ultimate_out <- matrix(0,length(date_list)+1,n_samples)
  Ultimate_pred <- matrix(0,length(date_list)+1,n_samples)
  resMeanTot_out <- matrix(0,1,n_samples)
  resMeanTot_pred <- matrix(0,1,n_samples)
  
  ## Convert reporting triangle to cumulative counts____________________________
  Latest <- getLatestCumulative(incr2cum(as.triangle(lda)))
  Latest <- c(Latest, total = sum(Latest))
  
  ## Generate posterior sample of predictions___________________________________
  for (i in 1:n_samples)
  { 
    if (model == "nonpara"){
      fits <- linv(Xp %*% mh_sample[i,])
    } else if (model == "para"){
      if (target == "symptom"){
        fits <- linv(Xp %*% mh_sample[i,] + log(((pweibull(as.numeric(ldaOut$dev+0.5-1),shape = ldaOut$var1, scale = ldaOut$var2))-(pweibull(as.numeric(ldaOut$dev-0.5-1),shape = ldaOut$var1, scale = ldaOut$var2)))))
      } else if (target == "specimen"){
        # fits <- linv(Xp %*% mh_sample[i,] + log(((pweibull(as.numeric(ldaOut$dev+0.5-2),shape = ldaOut$var1, scale = ldaOut$var2))-(pweibull(as.numeric(ldaOut$dev-0.5-2),shape = ldaOut$var1, scale = ldaOut$var2)))))
        fits <- linv(Xp %*% mh_sample[i,] + log(((pweibull(as.numeric(ldaOut$dev+1-2),shape = ldaOut$var1, scale = ldaOut$var2))-(pweibull(as.numeric(ldaOut$dev-2),shape = ldaOut$var1, scale = ldaOut$var2)))))
      } else {
        print("unknown target variable")
      }
    } else {
      print("unknown model structure")
    }
    resMeanAy_out[sort(unique(ldaOut$origin)),i] <- tapply(fits,factor(ldaOut$origin),sum)
    fits <- rnbinom(n = length(fits),mu = fits,size = gamFit$family$getTheta(TRUE)) # generate random prediction sample
    resMeanAy_pred[sort(unique(ldaOut$origin)),i] <- tapply(fits,factor(ldaOut$origin),sum)
    # if (nrow(resMeanAy_out) > length(unique(ldaFit$origin))){
    if (nrow(resMeanAy_out) > 30){
        resMeanAy_out[1:(nrow(resMeanAy_out)-30),i] <- 0
    }
    resMeanTot_out[i] <- sum(resMeanAy_out[,i])
    IBNR <- round(c(resMeanAy_out[,i], resMeanTot_out[i]))
    Ultimate_out[,i] <- Latest + IBNR
    
    if (nrow(resMeanAy_pred) > 30){
      resMeanAy_pred[1:(nrow(resMeanAy_pred)-30),i] <- 0
    }
    resMeanTot_pred[i] <- sum(resMeanAy_pred[,i])
    IBNR <- round(c(resMeanAy_pred[,i], resMeanTot_pred[i]))
    Ultimate_pred[,i] <- Latest + IBNR
  }
  
  ## Calculate quantiles and add to data frame__________________________________
  central_pred <- rowQuantiles(Ultimate_pred,p=0.5)
  upper_pred <- rowQuantiles(Ultimate_pred,p=0.975)
  lower_pred <- rowQuantiles(Ultimate_pred,p=0.025)
  
  Ultimate <- rowQuantiles(Ultimate_pred,p=0.5)
  Ultimate_5 <- rowQuantiles(Ultimate_pred,p=0.025)
  Ultimate_95 <- rowQuantiles(Ultimate_pred,p=0.975)
  Ultimate_25 <- rowQuantiles(Ultimate_pred,p=0.25)
  Ultimate_75 <- rowQuantiles(Ultimate_pred,p=0.75)
  
  date_list <- date_list[1:length(date_list)]
  counts <- as.data.frame(date_list) %>% mutate(n = Latest[1:(length(Latest)-1)])
  counts$gam_preds <- Ultimate[1:(length(Ultimate)-1)]
  counts$gam_ci_5 <- Ultimate_5[1:(length(Ultimate)-1)]
  counts$gam_ci_95 <- Ultimate_95[1:(length(Ultimate)-1)]
  counts$gam_ci_25 <- Ultimate_25[1:(length(Ultimate)-1)]
  counts$gam_ci_75 <- Ultimate_75[1:(length(Ultimate)-1)]
  counts$central_pred <- central_pred[1:(length(central_pred)-1)]
  counts$upper_pred <- upper_pred[1:(length(upper_pred)-1)]
  counts$lower_pred <- lower_pred[1:(length(lower_pred)-1)]
  
  counts <- cbind(counts, Ultimate_out[1:(nrow(Ultimate_pred)-1),])
  return(list(counts))
}


symptom_nowcast_region <- function(max_date,df_dates,ldaFit,ldaOut,lda,n_samples = 1000, date_list_full,denom = 7){
## this function is the nonpara regional nowcast by symptom onset date
  
  k = round(length(unique(ldaFit$origin))/denom)
  
  ### Model formula
  gamFit <- gam(value ~ 
                  s(as.numeric(origin),k= k, bs = "cr") +
                  s(as.numeric(origin),as.factor(region),k= k, bs = "fs") +
                  s(as.numeric(dev), k = 10, bs = "cr") +
                  s((dow_dev),bs = "re") +
                  s(as.factor(region),bs = "re") +
                  s((dow_origin),bs = "re"),
                data = ldaFit, family = "nb")
  
  
  mh_sample <- gam.mh(gamFit,ns=10000,burn=1000,t.df=40,rw.scale=.25,thin=1)$bs
  linv <- gamFit$family$linkinv
  
  
  glmFit <- gamFit

  counts_out <- data.frame()
  region_list <- c("London","Not-London")
  for (region_ in region_list){
  temp <- ldaOut %>% filter(region == region_)
  templda <- lda %>% filter(region == region_)
  Xp <- predict(gamFit, temp, type = "lpmatrix")
  
  resMeanAy_out <- matrix(0,length(unique(ldaOut$origin)),n_samples)
  Ultimate_out <- matrix(0,length(unique(ldaOut$origin))+1,n_samples)
  resMeanTot_out <- matrix(0,1,n_samples)
  Latest <- getLatestCumulative(incr2cum(as.triangle(templda)))
  Latest <- Latest[-(1:(length(Latest) - length(unique(ldaOut$origin))))]
  Latest <- c(Latest, total = sum(Latest))

  
  for (i in 1:n_samples)
  { 
    fits <- linv(Xp %*% mh_sample[i,])
    fits <- rnbinom(n = length(fits),mu = fits,size = gamFit$family$getTheta(TRUE))
    resMeanAy_out[,i] <- tapply(fits,factor(temp$origin),sum)
    resMeanTot_out[i] <- sum(resMeanAy_out[,i])
    IBNR <- round(c(resMeanAy_out[,i], resMeanTot_out[i]))
    Ultimate_out[,i] <- Latest + IBNR
  }
  
  
  
  Ultimate <- rowQuantiles(Ultimate_pred,p=0.5)
  Ultimate_5 <- rowQuantiles(Ultimate_pred,p=0.10)
  Ultimate_95 <- rowQuantiles(Ultimate_pred,p=0.90)
  Ultimate_25 <- rowQuantiles(Ultimate_pred,p=0.25)
  Ultimate_75 <- rowQuantiles(Ultimate_pred,p=0.75)
  
  date_list <- date_list_full[2:length(date_list_full)]
  counts <- as.data.frame(date_list) %>% mutate(n = Latest[1:(length(Latest)-1)])
  counts$gam_preds <- Ultimate[1:length(Ultimate)-1]
  counts$gam_ci_5 <- Ultimate_5[1:length(Ultimate)-1]
  counts$gam_ci_95 <- Ultimate_95[1:length(Ultimate)-1]
  counts$gam_ci_25 <- Ultimate_25[1:length(Ultimate)-1]
  counts$gam_ci_75 <- Ultimate_75[1:length(Ultimate)-1]
  
  counts <- cbind(counts, Ultimate_out[1:nrow(Ultimate_pred)-1,])
  counts$region <- region_
  counts_out <- rbind(counts,counts_out)
  }
  return(list(counts_out))
}


nowcast_specimen_region <- function(max_date,df_dates,ldaFit,ldaOut,lda,n_samples = 1000, date_list_full,denom = 7){
## this function is the nonpara regional nowcast by specimen date
  
  k = round(length(unique(ldaFit$origin))/denom)
  
  ### Model formula
  gamFit <- gam(value ~ 
                  s(as.numeric(origin),k= k, bs = "cr") +
                  s(as.numeric(origin),as.factor(region),k= k, bs = "fs") +
                  s(as.numeric(dev), by = as.factor(weekend_reporting), k = 10, bs = "cr") +
                  s(as.factor(region),bs = "re") +
                  s((dow_dev),bs = "re") +
                  s((dow_origin),bs = "re"),
                data = ldaFit, family = "nb")
  
  
  mh_sample <- gam.mh(gamFit,ns=10000,burn=1000,t.df=40,rw.scale=.25,thin=1)$bs
  linv <- gamFit$family$linkinv
  
  
  glmFit <- gamFit

  counts_out <- data.frame()
  region_list <- c("London","Not-London")
  for (region_ in region_list){
    temp <- ldaOut %>% filter(region == region_)
    templda <- lda %>% filter(region == region_)
    Xp <- predict(gamFit, temp, type = "lpmatrix")
    
    resMeanAy_out <- matrix(0,length(unique(ldaOut$origin)),n_samples)
    Ultimate_out <- matrix(0,length(unique(ldaOut$origin))+1,n_samples)
    resMeanTot_out <- matrix(0,1,n_samples)
    Latest <- getLatestCumulative(incr2cum(as.triangle(templda)))
    Latest <- Latest[-(1:(length(Latest) - length(unique(ldaOut$origin))))]
    Latest <- c(Latest, total = sum(Latest))

    
    for (i in 1:n_samples)
    { 
      fits <- linv(Xp %*% mh_sample[i,])
      fits <- rnbinom(n = length(fits),mu = fits,size = gamFit$family$getTheta(TRUE))
      resMeanAy_out[,i] <- tapply(fits,factor(temp$origin),sum)
      resMeanTot_out[i] <- sum(resMeanAy_out[,i])
      IBNR <- round(c(resMeanAy_out[,i], resMeanTot_out[i]))
      Ultimate_out[,i] <- Latest + IBNR
    }
    
    
    
    Ultimate <- rowQuantiles(Ultimate_pred,p=0.5)
    Ultimate_5 <- rowQuantiles(Ultimate_pred,p=0.10)
    Ultimate_95 <- rowQuantiles(Ultimate_pred,p=0.90)
    Ultimate_25 <- rowQuantiles(Ultimate_pred,p=0.25)
    Ultimate_75 <- rowQuantiles(Ultimate_pred,p=0.75)
    
    
    
    date_list <- date_list_full[2:length(date_list_full)]
    counts <- as.data.frame(date_list) %>% mutate(n = Latest[1:(length(Latest)-1)])
    counts$gam_preds <- Ultimate[1:length(Ultimate)-1]
    counts$gam_ci_5 <- Ultimate_5[1:length(Ultimate)-1]
    counts$gam_ci_95 <- Ultimate_95[1:length(Ultimate)-1]
    counts$gam_ci_25 <- Ultimate_25[1:length(Ultimate)-1]
    counts$gam_ci_75 <- Ultimate_75[1:length(Ultimate)-1]
    
    counts <- cbind(counts, Ultimate_out[1:nrow(Ultimate_pred)-1,])
    counts$region <- region_
    counts_out <- rbind(counts,counts_out)
  }
  return(list(counts_out))
}

