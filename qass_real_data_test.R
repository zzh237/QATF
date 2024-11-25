# Test QASS without parallelization 

# QATF Real data
# real_data.R

library(tidyverse)
library(foreach)
library(doParallel)
source("algorithms.R")

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
  summarise(across(everything(), \(x) min(x, na.rm = TRUE))) %>%
  unlist()

max_vals <- Combined_dat %>%
  select(2:10) %>%
  summarise(across(everything(),  \(x) max(x, na.rm = TRUE))) %>%
  unlist()

Scaled_Centered <- Combined_dat %>% 
  mutate(across(2:10, ~ (. - min_vals[cur_column()]) / (max_vals[cur_column()] - min_vals[cur_column()]))) %>%
  mutate(Ladder.score = Ladder.score - mean(Ladder.score, na.rm = TRUE))

splits <- sample(rep(1:10, length.out = nrow(Scaled_Centered)))

num <- 2
test_set_mean <- numeric(num)
test_set_coverage9 <- numeric(num)
test_set_coverage8 <- numeric(num)
test_set_coverage6 <- numeric(num)

test_set_mean_sorted <- numeric(num)
test_set_coverage9_sorted <- numeric(num)
test_set_coverage8_sorted <- numeric(num)
test_set_coverage6_sorted <- numeric(num)

for (i in 1:num) {
  source("algorithms.R")
  
  test_indices <- which(splits == i)
  test_data <- Scaled_Centered[test_indices, ]
  train_data <- Scaled_Centered[-test_indices, ]
  n = nrow(train_data)
  
  # Model fitting
  res5 <- qass_cv(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.5, prints = FALSE)
  lambda5 <- res5$LAMBDAS[which.min(rev(res5$MEANS))]
  fit5 <- fit_qass(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.5, lambda5)
  pred5 <- predict_fit(train_data[, 2:10], fit5$fit, test_data[, 2:10])
  
  # Coverage calculations
  res95 <- qass_cv(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.95, prints = FALSE)
  lambda95 <- res95$LAMBDAS[which.min(res95$MEANS)]
  fit95 <- fit_qass(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.95, lambda95)
  pred95 <- predict_fit(train_data[, 2:10], fit95$fit, test_data[, 2:10])
  
  res9 <- qass_cv(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.9, prints = FALSE)
  lambda9 <- res9$LAMBDAS[which.min(res9$MEANS)]
  fit9 <- fit_qass(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.9, lambda9)
  pred9 <- predict_fit(train_data[, 2:10], fit9$fit, test_data[, 2:10])
  
  res8 <- qass_cv(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.8, prints = FALSE)
  lambda8 <- res8$LAMBDAS[which.min(res8$MEANS)]
  fit8 <- fit_qass(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.8, lambda8)
  pred8 <- predict_fit(train_data[, 2:10], fit8$fit, test_data[, 2:10])
  
  res2 <- qass_cv(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.2, prints = FALSE)
  lambda2 <- res2$LAMBDAS[which.min(res2$MEANS)]
  fit2 <- fit_qass(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.2, lambda2)
  pred2 <- predict_fit(train_data[, 2:10], fit2$fit, test_data[, 2:10])
  
  res1 <- qass_cv(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.1, prints = FALSE)
  lambda1 <- res1$LAMBDAS[which.min(res1$MEANS)]
  fit1 <- fit_qass(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.1, lambda1)
  pred1 <- predict_fit(train_data[, 2:10], fit1$fit, test_data[, 2:10])
  
  res05 <- qass_cv(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.05, prints = FALSE)
  lambda05 <- res05$LAMBDAS[which.min(res05$MEANS)]
  fit05 <- fit_qass(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.05, lambda05)
  pred05 <- predict_fit(train_data[, 2:10], fit05$fit, test_data[, 2:10])
  
  test_set_mean[i] <- MSE(pred5, test_data[, 11])
  test_set_coverage9[i] <- mean(test_data[, 11] <= pred95 & test_data[, 11] >= pred05)
  test_set_coverage8[i] <- mean(test_data[, 11] <= pred9 & test_data[, 11] >= pred1)
  test_set_coverage6[i] <- mean(test_data[, 11] <= pred8 & test_data[, 11] >= pred2)
  
  
  ### Sorting all 7 quantiles together ###
  # pred_matrix <- cbind(pred95, pred9, pred8, pred5, pred2, pred1, pred05)
  # pred_sorted <- t(apply(pred_matrix, 1, sort, decreasing = TRUE))
  
  # pred95_sorted <- pred_sorted[, 1]
  # pred9_sorted  <- pred_sorted[, 2]
  # pred8_sorted  <- pred_sorted[, 3]
  # pred5_sorted  <- pred_sorted[, 4]
  # pred2_sorted  <- pred_sorted[, 5]
  # pred1_sorted  <- pred_sorted[, 6]
  # pred05_sorted <- pred_sorted[, 7]
  
  # test_set_mean_sorted <- MSE(pred5_sorted, test_data[, 11])
  # test_set_coverage9_sorted <- mean(test_data[, 11] <= pred95_sorted & test_data[, 11] >= pred05_sorted)
  # test_set_coverage8_sorted <- mean(test_data[, 11] <= pred9_sorted & test_data[, 11] >= pred1_sorted)
  # test_set_coverage6_sorted <- mean(test_data[, 11] <= pred8_sorted & test_data[, 11] >= pred2_sorted)
  
  ### Sorting quantile pairs only ###
  test_set_mean_sorted[i] <- test_set_mean
  test_set_coverage9_sorted[i] <- mean(test_data[, 11] <= max(pred05, pred95) & test_data[, 11] >= min(pred05, pred95))
  test_set_coverage8_sorted[i] <- mean(test_data[, 11] <= max(pred1, pred9) & test_data[, 11] >= min(pred1, pred9))
  test_set_coverage6_sorted[i] <- mean(test_data[, 11] <= max(pred2, pred8) & test_data[, 11] >= min(pred2, pred8))
}

# Convert the results into a data frame
average_metrics <- data.frame(
  Metric = c("Mean MSE", "Mean Coverage 0.9", "Mean Coverage 0.8", "Mean Coverage 0.6", 
             "Mean MSE (sorted)", "Mean Coverage 0.9 (sorted)", 
             "Mean Coverage 0.8 (sorted)", "Mean Coverage 0.6 (sorted)"),
  Value = c(mean(test_set_mean), mean(test_set_coverage9), mean(test_set_coverage8), mean(test_set_coverage6), 
            mean(test_set_mean_sorted), mean(test_set_coverage9_sorted), mean(test_set_coverage8_sorted), mean(test_set_coverage6_sorted))
)
write.csv(average_metrics, "average_metrics_s_test.csv", row.names = FALSE)

cat("Average MSE: ", mean_mse, "\n",
    "Average Coverage (95% CI): ", mean_coverage9, "\n",
    "Average Coverage (90% CI): ", mean_coverage8, "\n",
    "Average Coverage (80% CI): ", mean_coverage6, "\n")



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



