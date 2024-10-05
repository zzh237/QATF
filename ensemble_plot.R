library(tidyverse)
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
  summarise(across(everything(), min, na.rm = TRUE)) %>%
  unlist()

max_vals <- Combined_dat %>%
  select(2:10) %>%
  summarise(across(everything(), max, na.rm = TRUE)) %>%
  unlist()

Scaled_Centered <- Combined_dat %>% 
  mutate(across(2:10, ~ (. - min_vals[cur_column()]) / (max_vals[cur_column()] - min_vals[cur_column()]))) %>%
  mutate(Ladder.score = Ladder.score - mean(Ladder.score, na.rm = TRUE))


lambdas <- read.csv("happiness/lambda.csv")
lambda5 <- lambdas$Lambda5
lambda1 <- lambdas$Lambda1
lambda9 <- lambdas$Lambda9
splits <- read.csv("happiness/splits.csv")


png(filename = "ensemble_plot.png", width = 1000, height = 1000)
par(mfrow = c(3, 3))
first_plot <- TRUE

for (i in 1:10) {
  
  test_indices <- which(splits == i)
  test_data <- Scaled_Centered[test_indices, ]
  train_data <- Scaled_Centered[-test_indices, ]
  n = nrow(train_data)
  
  outq <- fit_qatf(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.5, 0, lambda5[i] / 10)
  outq1 <- fit_qatf(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.1, 0, lambda1[i] / 10)
  outq9 <- fit_qatf(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.9, 0, lambda9[i] / 10)
  for (j in 1:9) {
    if (first_plot) {
      plot(outq[[2]][, j][outq[[3]][j, ]] + mean(Combined_dat[, 11]),
           xaxt = 'n', yaxt = 'n', main = colnames(Combined_dat)[j + 1],
           xlab = "", ylab = "", lwd = 1, type = "l", ylim = c(4, 7))

      x_range <- seq(1, length(outq[[2]][, j][outq[[3]][j, ]]), length.out = 5)
      x_original <- seq(0, 1, length.out = 5) * (max_vals[j] - min_vals[j]) + min_vals[j]
      axis(1, at = x_range, labels = round(x_original, 2))
      
      y_range <- c(4.5, 5, 5.5, 6, 6.5)
      axis(2, at = y_range, labels = y_range) 
      
      # Add upper and lower quantiles
      lines(outq1[[2]][, j][outq1[[3]][j, ]] + mean(Combined_dat[, 11]), lwd = 1, col = "purple")
      lines(outq9[[2]][, j][outq9[[3]][j, ]] + mean(Combined_dat[, 11]), lwd = 1, col = "orange")
      
    } else {
      par(mfg = c((j - 1) %/% 3 + 1, (j - 1) %% 3 + 1))
      
      lines(outq[[2]][, j][outq[[3]][j, ]] + mean(Combined_dat[, 11]), lwd = 1)
      
      # Add upper and lower quantiles
      lines(outq1[[2]][, j][outq1[[3]][j, ]] + mean(Combined_dat[, 11]), lwd = 1, col = "purple")
      lines(outq9[[2]][, j][outq9[[3]][j, ]] + mean(Combined_dat[, 11]), lwd = 1, col = "orange")
    }
  }
  first_plot <- FALSE
}

# Close the PNG device and save the plot
dev.off()




png(filename = "shade_plot.png", width = 1000, height = 1000)
par(mfrow = c(3, 3))
first_plot <- TRUE

for (i in 1:10) {
  
  test_indices <- which(splits == i)
  test_data <- Scaled_Centered[test_indices, ]
  train_data <- Scaled_Centered[-test_indices, ]
  n = nrow(train_data)
  
  # Fit the quantile models
  outq <- fit_qatf(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.5, 0, lambda5[i] / 100)
  outq1 <- fit_qatf(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.1, 0, lambda1[i] / 100)
  outq9 <- fit_qatf(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.9, 0, lambda9[i] / 100)
  
  for (j in 1:9) {
    y_median <- outq[[2]][, j][outq[[3]][j, ]] + mean(Combined_dat[, 11])
    y_lower <- outq1[[2]][, j][outq1[[3]][j, ]] + mean(Combined_dat[, 11])
    y_upper <- outq9[[2]][, j][outq9[[3]][j, ]] + mean(Combined_dat[, 11])
    
    if (any(y_lower > y_upper)) {
      indices_to_adjust <- which(y_lower > y_upper)
      y_lower[indices_to_adjust] <- y_median[indices_to_adjust]
      y_upper[indices_to_adjust] <- y_median[indices_to_adjust]
    }
    
    if (first_plot) {
      plot(NULL, xaxt = 'n', yaxt = 'n', main = colnames(Combined_dat)[j + 1],
           xlab = "", ylab = "", type = "n", ylim = c(4, 7), xlim = c(1, length(y_lower)))
      
      x_range <- seq(1, length(y_lower), length.out = 5)
      x_original <- seq(0, 1, length.out = 5) * (max_vals[j] - min_vals[j]) + min_vals[j]
      axis(1, at = x_range, labels = round(x_original, 2))
      
      y_range <- c(4.5, 5, 5.5, 6, 6.5)
      axis(2, at = y_range, labels = y_range)
      
    
      polygon(c(1:length(y_lower), rev(1:length(y_upper))),
              c(y_lower, rev(y_upper)),
              col = rgb(0.7, 0.7, 0.7, 0.5), border = NA)
      
    } else {
      par(mfg = c((j - 1) %/% 3 + 1, (j - 1) %% 3 + 1))
      
      polygon(c(1:length(y_lower), rev(1:length(y_upper))),
              c(y_lower, rev(y_upper)),
              col = rgb(0.7, 0.7, 0.7, 0.5), border = NA)
    }
  }
  first_plot <- FALSE
}
dev.off()
