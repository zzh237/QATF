# QATF Real data New
# real_data.R

library(tidyverse)
library(foreach)
library(doParallel)


Combined_dat <- read.csv("happiness/Combined_dat.csv")
# par(mfrow = c(3, 3), mar = c(2, 2, 2, 2), oma = c(0, 0, 0, 0))
# for (i in 1:9) {
#   # Generate the plot with the title set to the column name
#   plot(sort(Combined_dat[, 1 + i]), 
#        main = colnames(Combined_dat)[1 + i], 
#        xlab = "", ylab = "", lwd = 1)
# }

min_vals <- Combined_dat %>%
  select(2:10) %>%
  summarise(across(everything(), min, na.rm = TRUE)) %>%
  unlist()

max_vals <- Combined_dat %>%
  select(2:10) %>%
  summarise(across(everything(), max, na.rm = TRUE)) %>%
  unlist()

Scaled_Centered <- Combined_dat %>% 
  mutate(across(2:10, ~ (. - min_vals[cur_column()]) / (max_vals[cur_column()] - min_vals[cur_column()]))) %>%
  mutate(Ladder.score = Ladder.score - mean(Ladder.score, na.rm = TRUE))

# Set up parallel backend
num_cores <- detectCores() - 1  # Use all cores except one
cl <- makeCluster(num_cores)
registerDoParallel(cl)

splits <- sample(rep(1:10, length.out = nrow(Scaled_Centered)))

# Parallel loop using foreach
results <- foreach(i = 1:10) %dopar% {
  source("algorithms.R")
  
  test_indices <- which(splits == i)
  test_data <- Scaled_Centered[test_indices, ]
  train_data <- Scaled_Centered[-test_indices, ]
  n = nrow(train_data)
  
  # Model fitting
  res5 <- qatf_cv(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.5, 0, prints = FALSE)
  lambda5 <- res5$LAMBDAS[which.min(rev(res5$MEANS))]
  fit5 <- fit_qatf(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.5, 0, lambda5)
  test_set_mean <- MSE(predict_fit(train_data[, 2:10], fit5$fit, test_data[, 2:10]), test_data[, 11])


  res9 <- qatf_cv(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.9, 0, prints = FALSE)
  lambda9 <- res9$LAMBDAS[which.min(rev(res9$MEANS))]
  fit9 <- fit_qatf(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.9, 0, lambda9)
  pred9 <- predict_fit(train_data[, 2:10], fit9$fit, test_data[, 2:10])

  res1 <- qatf_cv(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.1, 0, prints = FALSE)
  lambda1 <- res1$LAMBDAS[which.min(rev(res1$MEANS))]
  fit1 <- fit_qatf(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.1, 0, lambda1)
  pred1 <- predict_fit(train_data[, 2:10], fit1$fit, test_data[, 2:10])
  test_set_coverage8 <- mean(test_data[, 11] <= pred9 & test_data[, 11] >= pred1)

  list(
    test_set_mean = test_set_mean,
    test_set_coverage8 = test_set_coverage8,
    lambda5 = lambda5, 
    lambda1 = lambda1,
    lambda9 = lambda9
  )

}

# Stop cluster
stopCluster(cl)

# Collect and calculate the average MSE and coverage metrics
mean_mse <- mean(sapply(results, function(x) x$test_set_mean))
mean_coverage8 <- mean(sapply(results, function(x) x$test_set_coverage8))

# Convert the results into a data frame
average_metrics <- data.frame(
  Metric = c("Mean MSE", "Mean Coverage 0.8"),
  Value = c(mean_mse, mean_coverage8)
)

write.csv(average_metrics, "average_metrics.csv", row.names = FALSE)
write.csv(data.frame(Split = splits), file = "splits.csv", row.names = FALSE)

lambda_df <- data.frame(Lambda5 = sapply(results, function(x) x$lambda5), 
                        Lambda1 = sapply(results, function(x) x$lambda1),
                        Lambda9 = sapply(results, function(x) x$lambda9))
write.csv(lambda_df, file = "lambda.csv", row.names = FALSE)



# Their predictors were: 
# Satisfaction with freedom of choice
# Satisfaction with job
# Statisfaction with community
# Trust in national government
# % Rural population
# % Females with secondary education
# Mortality rate, under 5
# Log gross national income
# % Labor force unemployed
# Life expectancy at birth (years) 
# % internet users
# Lof sceintific journal articles published

# Our predictors
# Log GDP per-capita (WH)
# Social Support (WH)
# Healthy life expectancy at birth (WH)
# Freedom to make life choices (WH)
# Perceptions of corruption (WH)
# Years of Compulsory Education (WDI)
# % Rural population (WDI)
# Log Scientific journal articles published (WDI)
# % internet users (WDI)



# # Plots
# n = nrow(Scaled_Centered)
# lambda_qres <- qatf_cv(t(Scaled_Centered[, 2:10]), Scaled_Centered[, 11], n, 9, 0.5, 0)
# lambda <- lambda_qres$LAMBDAS[which.min(lambda_qres$MEANS)]
# outq <- fit_qatf(t(Scaled_Centered[, 2:10]), Scaled_Centered[, 11], n, 9, 0.5, 0, lambda)
# 
# for (i in 1:9) {
#   # Generate the plot
#   plot(outq[[2]][, i][outq[[3]][i, ]] + mean(Combined_dat[, 11]), 
#        xaxt = 'n', yaxt = 'n', main = paste("Plot", i), xlab = "", ylab = "", lwd = 1,
#        ylim = c(0, 10))
#   
#   # # If I need labels
#   # # Calculate the original x-axis labels for this plot
#   # x_original <- seq(0, 1, length.out = 5) * (max_vals[i] - min_vals[i]) + min_vals[i]
#   # # Add custom x-axis labels with the original values
#   # axis(1, at = seq(0, 1, length.out = 5), labels = round(x_original, 2))
# }
