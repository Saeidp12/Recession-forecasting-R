library(tidyverse)
library(lubridate)
library(naivebayes)
library(WeightedROC)
library(R.matlab)

# This is the data set generated in MATLAB
df_matfile <- readMat('Matlab/freddata.mat')

df_matfile_data <- df_matfile[[1]]
df_matfile_dates <- df_matfile[[2]]
df_matfile_names <- c('Dates','NAPMPI', 'PAYEMS', 'PAYEMS_RT','S&P 500', 'T10YFFM', 'GS10', 'TB3MS', 'GS10_TB3MS')
df_matfile_tcode <- df_matfile$tcode

#--------------- prepare_missing.m
prepare_missing <- function(df, t_code){
  m <- ncol(df)
  if(length(t_code)!= m){
    stop('Number of dimensions in the raw data differs from the length of code vectors')
  }
  transform_func <- function(x){
    col <- which(df == x, arr.ind=T)[1,2]
    name <- names(df)[col]
    code <- t_code[col]
    n <- length(x)
    small <- 1e-06
    if(code == 1){
      return(x)
    }
    if(code == 2){
      return(c(NA, diff(x)))
    }
    if(code == 3){
      
      return(c(rep(NA, 2), diff(x, lag=1, differences=2)))
    }
    if(code == 4){
      if(min(x) < small){
        stop(cat("small values in column ", name))
      }
      else{
        return(log(x))
      }
    }
    if(code == 5){
      if(min(x) > small){
        return(c(NA, diff(log(x))))
      }
    }
    if(code == 6){
      if(min(x) > small){
        return(c(NA,NA, diff(log(x), lag=1, differences=2)))
      }
      else{stop(cat("small values in column ", name))}
    }
    if(code == 7){
      a = (x[2:n] - x[1:n-1])/x[1:n-1]
      return(c(NA,NA, diff(a)))
    }
  }
  df <- (apply(df, 2, transform_func))
  return(df[-c(1,2),])
}


#------------------- remove_outliers.m
remove_outliers <- function(df){
  medians <- apply(df, 2, median, na.rm=T)
  threshold <- function(x){
    col <- which(df == x, arr.ind=T)[1,2]
    out <- abs(x - medians[col])>10*IQR(df[,col], na.rm = T)
    indices <- which(out==T)
    x[indices] <- NA
    return(x)
  }
  return(apply(df, 2, threshold))
}

#------------------- generate data

generate <- function(df, t_code){
  # select only numeric variables
  df <- dplyr::select_if(as_tibble(df), is.numeric) %>% as.data.frame
  # prepare missing
  df_1 <- prepare_missing(df, t_code)
  # remove outliers
  df_2 <- remove_outliers(df_1)
  
  return(df_2)
}



#-------------------- Main Recession Function
#--------------- Lag Structure Settings
n_lag = 10
#--------------- Forecast Settings
forecast_horizons = c(0, 3, 6, 9, 12)
n_horz = length(forecast_horizons)

forecast_blind = 18
use_alt_labels = TRUE # use recession within horizons instead of recession at horizons
forecast_start_date = 1973
forecast_end_date = Inf # End of the data

#--------------- Model settings
use_model = 3
#1 : Logistic Regression
#2 : Naive Bayes
#3 : Modified Naive Bayes

#--------------- Performance Evaluation Settings
use_timedep_weights = TRUE # using time-dependent error weights instead of uniform weights

#-------------- macro data
df <- read.csv('Matlab/macro_data.csv', header=T)
names(df) <- df_matfile_names

transform_code <- df[1,-1]
df_raw <- df[-1,]
dates_raw <- df_raw$Dates
df_raw <- generate(df_raw[,-1], transform_code)

df_raw <- data.frame(Dates=dates_raw[-c(1,2)], df_raw)

# transforming the dates
df_raw$Dates <- as.Date(df_raw$Dates,"%m/%d/%Y")

# select a subset of macro_vars
selected_cols <- c('Dates','NAPMPI','PAYEMS_RT', 'S.P.500', 'T10YFFM')
macro_selected <- df_raw[,selected_cols]

#--------------- recession data
rec_df <- read.csv('Matlab/NBER_recessions.csv', header=T)
rec_df$Date <- as.Date(rec_df$Date)
rec_labels <- rec_df$Recession
rec_dates <- rec_df$Date

#--------------- Define and generate lags of macro variables

dates <- macro_selected$Dates

macro_lagged <- matrix(NA, nrow=nrow(macro_selected), ncol=(n_lag+1)*(ncol(macro_selected)-1))
for(i in (n_lag+1):nrow(macro_selected)){
  data_lagged = macro_selected[(i-n_lag):i,-1]
  macro_lagged[i,] = c(as.matrix(data_lagged))
}

macro_final <- macro_lagged[(n_lag+1):nrow(macro_lagged),]

macro_final <- data.frame(Dates=dates[-c(1:10)], macro_final)

rec_df_common <- rec_df %>% filter(Date %in% macro_final$Dates)
rec_labels_common <- rec_df_common$Recession
macro_common <- macro_final
n_obs <- nrow(macro_common)


#--------------- Create Data Structures and Indices for Forecast

forecast_start_date <- as.Date('1973-01-01')
forecast_end_date <- max(macro_common$Dates) # selected or common

forecast_start_idx <- which(macro_common$Dates == forecast_start_date)
forecast_end_idx <- which(macro_common$Dates == forecast_end_date)

forecast_probs <- matrix(NA, n_obs+max(forecast_horizons), n_horz)

#------------ Iterate through forecast dates and horizons, generating forecast probabilities along the way

for (forecast_i in forecast_start_idx:forecast_end_idx){
  pct_progress = 100*(forecast_i-forecast_start_idx)/(forecast_end_idx-forecast_start_idx)
  for(horizon_j in 1:n_horz){
    forecast_horizon <- forecast_horizons[horizon_j]
    
    train_idx <-  1:(forecast_i-1-forecast_blind-forecast_horizon)
    train_X <- macro_common[train_idx,-1] # remove dates
    train_y <- rec_labels_common[train_idx+forecast_horizon]
    
    test_idx <- forecast_i
    test_X <- macro_common[test_idx,-1]
    
    if(use_alt_labels){
      alt_train_y = train_y;
      for(d_i in (1+forecast_horizon):length(train_y)){
        if(train_y[d_i-forecast_horizon] == 0 && sum(train_y[(d_i-forecast_horizon):d_i])){
          alt_train_y[d_i] = 1
        }
      }
      train_y = alt_train_y;
    }
    alt_train_y <- train_y # needed the non-factorized version for Naive Bayes
    train_y <- factor(train_y)
    if(use_model == 1){
      # logistic regression
      df <- data.frame(train_X, y=train_y)
      glm_fit <- glm(y~., data=df, family=binomial(link='logit'))
      prob <- predict(glm_fit, test_X, type='response')
      probs <- c(1-prob, prob)
    }
    if(use_model == 2){
      # naive Bayes
      df <- data.frame(train_X, y=train_y) 
      nb_fit <- naive_bayes(y~., data=df)
      test_X <- as.data.frame(test_X)
      prob <- predict(nb_fit, test_X, type='prob')
      probs <- prob
    }
    if(use_model == 3){
      last_pos_prob <- forecast_probs[forecast_i+forecast_horizon-1,horizon_j]
      last_neg_prob <- 1-last_pos_prob
      
      total_neg <- sum(train_y[1:(length(train_y)-1)] == 0)
      total_pos <- sum(train_y[1:(length(train_y)-1)] == 1)
      
      pos_after_neg <-  sum(train_y[which(train_y[1:(length(train_y)-1)] == 0)+1]==1)
      pos_after_neg_prob = pos_after_neg/total_neg
      
      neg_after_neg <- sum(train_y[which(train_y[1:(length(train_y)-1)] == 0)+1]==0)
      neg_after_neg_prob <- neg_after_neg/total_neg
      
      n_y <- length(train_y)
      pos_after_pos_prob <- sum(train_y[which(train_y[1:(n_y - 1)]==1)+1]==1)/sum(train_y[1:(n_y-1)]==1);
      neg_after_pos_prob <- sum(train_y[which(train_y[1:(n_y - 1)]==1)+1]==0)/sum(train_y[1:(n_y-1)]==1);
      
      pos_prior <- pos_after_neg_prob*last_neg_prob+pos_after_pos_prob*last_pos_prob
      neg_prior <- neg_after_neg_prob*last_neg_prob+neg_after_pos_prob*last_pos_prob 
      
      modified_priors = c(neg_prior, pos_prior)
      df <- data.frame(train_X, y=train_y) 
      test_X <- as.data.frame(test_X)
      if(any(is.na(modified_priors))){
        model <- naive_bayes(y~., data=df)
        prob <- predict(model, test_X, type='prob')
        probs <- prob
      }
      if(!any(is.na(modified_priors))){
        model <- naive_bayes(y~., data=df, prior = modified_priors)
        prob <- predict(model,test_X, type='prob')
        probs <- prob
      }
    }
    forecast_probs[forecast_i+forecast_horizon,horizon_j] <- probs[2]
  }
}

#----------- Construct Forecast Weights for Performance Calculations
forecast_weights <- rep(1, nrow(forecast_probs))
if(use_timedep_weights){
  n_rec <- length(rec_labels_common)
  label_edges <-  ((rec_labels_common[1:(n_rec-1)]==0) & (rec_labels_common[2:n_rec]==1)) | 
    ((rec_labels_common[1:(n_rec-1)]==1) & (rec_labels_common[2:n_rec]==0))
  label_edges_idxs <- c(0, which(label_edges), length(label_edges)+1)
  for(i in 1:(length(label_edges_idxs)-1)){
    idx_range <- (label_edges_idxs[i]+1):label_edges_idxs[i+1]
    range_size <- length(idx_range)
    
    weights <- dnorm(seq(0,1, by=(1/(range_size-1))), 0.5, 1)
    scaled_weights <- (weights - min(weights))/(max(weights)-min(weights))
    
    if(i==1){
      scaled_weights[1:round(length(scaled_weights)/2)]=max(scaled_weights)
    }
    if(i==(length(label_edges_idxs)-1)){
      scaled_weights[round(length(scaled_weights)/2):length(scaled_weights)]=max(scaled_weights)
    }
    forecast_weights[idx_range] = scaled_weights*range_size/sum(scaled_weights)
  }
  forecast_weights[which(forecast_weights==0)] = 0.001
}

#------------ Assemble forecasts into a table for output and review
# This is the matlab file which is loaded for comparison
# forecast_table_m <- read.table('Matlab/forecast_table.txt', sep=',', header=TRUE)


forecast_dates <- seq.Date(macro_common$Dates[1], forecast_end_date %m+% years(1), by='month')

forecast_years <- lubridate::year(forecast_dates)
forecast_months <- lubridate::month(forecast_dates)

colnames(forecast_probs) <- c('H0','H3','H6','H9','H12')

forecast_table <- data.frame(Year = forecast_years, Month=forecast_months, forecast_probs)
forecast_table$recession <- c(rec_labels_common, rep(NA, 12))

#----------- Calculate and print summary statistics for forecast performance
performance_table <- data.frame(Statistic=c('MAE','MSE','WAUC','ACC','F1'), H0=rep(0, 5),
                                            H3 = rep(0, 5), H6=rep(0, 5),
                                            H9=rep(0, 5), H12=rep(0, 5))

for(horizon_i in 1:n_horz){
  horizon_probs <- forecast_probs[, horizon_i]
  horizon_idx <- which(!is.na(horizon_probs))
  
  rec_outcomes <- forecast_table$recession
  rec_idx <- which(!is.na(rec_outcomes))
  
  common_idx = base::intersect(horizon_idx,rec_idx)
  horizon_probs_common = horizon_probs[common_idx]
  rec_outcomes_common = rec_outcomes[common_idx]
  
  forecast_weights_common = forecast_weights[common_idx]
  
  MAE <- mean(abs(horizon_probs_common-rec_outcomes_common)*forecast_weights_common)
  MSE <- mean(((horizon_probs_common-rec_outcomes_common)^2)*forecast_weights_common)
  wroc <- WeightedROC(horizon_probs_common, rec_outcomes_common,
                      weight=forecast_weights_common)
  wauc <- WeightedAUC(wroc)
  
  TP <- sum((rec_outcomes_common==1 & horizon_probs_common>0.5)*forecast_weights_common)
  TN <- sum((rec_outcomes_common==0 & horizon_probs_common<0.5)*forecast_weights_common)
  FP <- sum((rec_outcomes_common==0 & horizon_probs_common>0.5)*forecast_weights_common)
  FN <- sum((rec_outcomes_common==1 & horizon_probs_common<0.5)*forecast_weights_common)

  precision <- TP/(TP+FP)
  recall<- TP/(TP+FN)
  
  F1 <- max(0,2*(precision*recall)/(precision+recall))
  ACC <- (TP+TN)/(TP+TN+FP+FN)
  
  performance_table[,horizon_i+1] <- c(MAE,MSE,wauc, ACC, F1)
}

performance_table$AVG <- apply(performance_table[,-1], 1, mean)

write.table(performance_table,'E:/performance_table.txt', row.names = F, sep=',')
write.table(forecast_table, 'E:/forecast_table.txt', row.names = F,sep=',')

pr <- read.table('Matlab/performance_table.txt', header=T, sep=',')
pr
performance_table
t1 <- read.table('Matlab/forecast_table.txt', header=T, sep=',')
t1 %>% tail
forecast_table %>% tail











