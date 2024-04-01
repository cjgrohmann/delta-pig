find_mortality_episodes <- function(global.dataset, X, Y, cohort.ID, padding, lambda, magnitude, duration, proportion) {
  
  ### Initialize and extract parameters
  X <- unname(unlist(global.dataset[, X]))
  Y <- unname(unlist(global.dataset[, Y]))
  cohort.ID <- unname(unlist(global.dataset[, cohort.ID]))
  padding <- padding
  lambda <- lambda
  magnitude <- magnitude
  proportion <- proportion
  output_list <- list()
  
  output_list[["parameters"]] <- data.frame(Parameter = c("padding", "lambda", "magnitude", "duration", "proportion"),
                                            Value = c(padding, lambda, magnitude, duration, proportion))
  
  ### Build new dataset
  df <- data.frame(COHORT_ID = cohort.ID,
                   X = X,
                   Y = Y)
  
  ### Summarize missing data by cohort
  missingness_report <- aggregate(cbind(X, Y) ~ COHORT_ID, data = df, function(x) {sum(is.na(x))}, na.action = NULL)
  missingness_report$X <- as.numeric(missingness_report$X)
  missingness_report$Y <- as.numeric(missingness_report$Y)
  colnames(missingness_report)[2:3] <- c("X_MISS", "Y_MISS")
  
  output_list[["missingness_report"]] <- missingness_report
  
  ### Merge missing data report with dataframe and filter
  df <- merge(df, missingness_report, by = "COHORT_ID")
  df <- df[df$X_MISS == 0 & df$Y_MISS == 0,]
  df <- df[, 1:3]
  
  output_list[["filtered_dataset"]] <- df
  rm(missingness_report)
  
  ### Split dataframe into subsets by cohort ID
  df_list <- split(df, f = df$COHORT_ID)
  output_list[["cohort_subsets"]] <- df_list
  rm(df)
  
  ### Add padding to cohort subsets
  padded_df_list <- list()
  cohort_names <- names(df_list)
  
  for (i in cohort_names) {
    
    subset_df <- df_list[[i]]
    min_days <- min(subset_df$X)
    max_days <- max(subset_df$X)
  
    leading_pad <- data.frame(COHORT_ID = rep("LEADING", times = padding),
                              X = (min_days + (padding * -1)):(min_days - 1),
                              Y = rep(0, times = padding))
    
    lagging_pad <- data.frame(COHORT_ID = rep("LAGGING", times = padding),
                              X = (max_days + 1):(max_days + (padding)),
                              Y = rep(0, times = padding))
    
    padded_df <- rbind.data.frame(leading_pad, subset_df, lagging_pad)
    
    padded_df_list[[i]] <- padded_df
    
  }
  
  rm(i, subset_df, min_days, max_days, leading_pad, lagging_pad, padded_df)
  
  ### Fit smoothing splines
  fit_list <- list()
  fitted_df_list <- list()
  
  for (i in cohort_names) {
    
    subset_df <- padded_df_list[[i]]
    SS_X <- subset_df$X
    SS_Y <- subset_df$Y
    
    #print(paste("Fitting smoothing spline regression equation for cohort", i, sep = " "))
    
    ## Fit spline regression model and extract fitted values
    fit <- smooth.spline(x = SS_X, y = SS_Y, lambda = lambda, cv = FALSE)
    EDF <- fit[["df"]]
    subset_df$EDF <- EDF
    subset_df$Y_HAT <- fitted(fit)
    
    ## Identify change points
    subset_df$Y_HAT_SLOPE <- c(NA, diff(subset_df$Y_HAT))
    subset_df$Y_HAT_SLOPE_LEAD1 <- dplyr::lead(subset_df$Y_HAT_SLOPE)
    subset_df$CHANGEPOINT <- ifelse(sign(subset_df$Y_HAT_SLOPE) != sign(subset_df$Y_HAT_SLOPE_LEAD1), 1, 0)
    subset_df$VALLEY <- ifelse(sign(subset_df$Y_HAT_SLOPE) == -1 & sign(subset_df$Y_HAT_SLOPE_LEAD1) == 1, 1, 0)
    subset_df$PEAK <- ifelse(sign(subset_df$Y_HAT_SLOPE) == 1 & sign(subset_df$Y_HAT_SLOPE_LEAD1) == -1, 1, 0)
    subset_df[subset_df$X == 1, "CHANGEPOINT"] <- 1
    subset_df[subset_df$X == 1, "VALLEY"] <- 1
    
    fit_list[[i]] <- fit
    fitted_df_list[[i]] <- subset_df
    
  }
  
  rm(i, subset_df, SS_X, SS_Y, fit, EDF)
  
  output_list[["fitted_models"]] <- fit_list
  
  ### Prepare final dataset for episode processing
  output_df <- do.call("rbind.data.frame", fitted_df_list)
  output_df <- output_df[output_df$COHORT_ID != "LEADING" & output_df$COHORT_ID != "LAGGING",]
  rownames(output_df) <- seq_len(nrow(output_df))
  
  ### Identify mortality sequences
  output_df <- output_df %>%
    dplyr::group_by(COHORT_ID) %>%
    dplyr::mutate(SEQUENCE_ID = cumsum(VALLEY),
                  SEQUENCE_ID = paste(COHORT_ID, SEQUENCE_ID, sep = "_")) %>%
    dplyr::ungroup(.)
  
  ### Identify start and end of mortality sequences
  output_df <- output_df %>%
    dplyr::group_by(SEQUENCE_ID) %>%
    dplyr::mutate(Y_HAT_SLOPE_MIN = rank(Y_HAT_SLOPE, ties.method = "first"),
                  Y_HAT_SLOPE_MAX = rank(-Y_HAT_SLOPE, ties.method = "last"),
                  START = ifelse(Y_HAT_SLOPE_MAX == 1, 1, 0),
                  END = ifelse(Y_HAT_SLOPE_MIN == 1, 1, 0)) %>%
    dplyr::ungroup(.) %>%
    dplyr::select(COHORT_ID, SEQUENCE_ID, X, Y, Y_HAT, Y_HAT_SLOPE, EDF, CHANGEPOINT, VALLEY, START, PEAK, END)
  
  ### Identify mortality subsequences
  output_df <- output_df %>%
    dplyr::mutate(START_NA = ifelse(START == 1, START, NA),
                  END_NA = ifelse(END == 1, END, NA),
                  SUBSEQUENCE_ID = coalesce(START_NA, END_NA),
                  SUBSEQUENCE_ID = ifelse(is.na(SUBSEQUENCE_ID) == TRUE, 0, SUBSEQUENCE_ID)) %>%
    dplyr::group_by(SEQUENCE_ID) %>%
    dplyr::mutate(SUBSEQUENCE_ID = cumsum(SUBSEQUENCE_ID),
                  SUBSEQUENCE_ID = ifelse(is.na(END_NA) == TRUE, SUBSEQUENCE_ID, 1)) %>%
    dplyr::ungroup(.) %>%
    dplyr::rename(SUBSEQUENCE_TYPE = SUBSEQUENCE_ID) %>%
    dplyr::mutate(SUBSEQUENCE_ID = paste(SEQUENCE_ID, SUBSEQUENCE_TYPE, sep = "_"),
                  EPISODE = ifelse(SUBSEQUENCE_TYPE == 1, 1, 0)) %>%
    dplyr::select(COHORT_ID:SEQUENCE_ID, SUBSEQUENCE_ID, X:END, SUBSEQUENCE_TYPE:EPISODE)
  
  ### Calculate subseqeunce statistics
  output_df <- output_df %>%
    dplyr::mutate(Y_GREATER_THAN_0 = ifelse(Y > 0, 1, 0)) %>%
    dplyr::group_by(SUBSEQUENCE_ID) %>%
    dplyr::mutate(DURATION = max(rank(SUBSEQUENCE_ID, ties.method = "first")),
                  Y_SUM = sum(Y),
                  Y_GREATER_THAN_0_SUM = sum(Y_GREATER_THAN_0)) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(MAGNITUDE = Y_SUM / DURATION,
                  PROPORTION = Y_GREATER_THAN_0_SUM / DURATION)
  
  ### Identify mortality episodes
  output_df <- output_df %>%
    dplyr::mutate(MAGNITUDE_FILTER = MAGNITUDE >= magnitude,
                  DURATION_FILTER = DURATION <= duration,
                  PROPORTION_FILTER = PROPORTION >= proportion) %>%
    dplyr::group_by(SUBSEQUENCE_ID) %>%
    dplyr::mutate(MAX_START = max(START),
                  MAX_PEAK = max(PEAK), 
                  MAX_END = max(END)) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(TYPE_FILTER = MAX_START == 1 & MAX_PEAK == 1 & MAX_END == 1,
                  GLOBAL_FILTER = MAGNITUDE_FILTER == TRUE & DURATION_FILTER == TRUE & PROPORTION_FILTER == TRUE & TYPE_FILTER == TRUE,
                  EPISODE = ifelse(GLOBAL_FILTER == TRUE & SUBSEQUENCE_TYPE == 1, 1, 0)) %>%
    dplyr::select(COHORT_ID:PROPORTION)
  
  ### Final columns
  output_df <- output_df %>%
    dplyr::mutate(PADDING_SETTING = padding,
                  LAMBDA_SETTING = lambda,
                  MAGNITUDE_SETTING = magnitude,
                  DURATION_SETTING = duration,
                  PROPORTION_SETTING = proportion)
    
  return(output_df)
  
}

