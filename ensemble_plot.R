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


plotnames <- c("Log-GDP per Capita", "Social Support", "Healthy Life Expectancy at Birth", 
               "Freedom to Make Life Choices", "Perceptions of Corruption", "Years of Compulsory Education", 
               "Internet Access Percentage", "Rural Population Percentage", "Log Journal Articles Published")

lambdas <- read.csv("happiness/lambdas.csv")
lambda5 <- lambdas$Lambda5
lambda1 <- lambdas$Lambda1
lambda9 <- lambdas$Lambda9
splits <- read.csv("happiness/splits2.csv")



#################
# Ensemble Plot # 
#################


png(filename = "ensemble_plot.png", width = 1000, height = 1000)
par(mfrow = c(3, 3))
first_plot <- TRUE

for (i in 1:10) {
  
  test_indices <- which(splits == i)
  test_data <- Scaled_Centered[test_indices, ]
  train_data <- Scaled_Centered[-test_indices, ]
  n = nrow(train_data)
  
  outq <- fit_qatf(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.5, 0, lambda5[i])
  outq1 <- fit_qatf(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.1, 0, lambda1[i])
  outq9 <- fit_qatf(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.9, 0, lambda9[i])
  for (j in 1:9) {
    if (first_plot) {
      plot(outq[[2]][, j][outq[[3]][j, ]] + mean(Combined_dat[, 11]),
           xaxt = 'n', yaxt = 'n', main = plotnames[j],
           xlab = "", ylab = "", lwd = 1, type = "l", ylim = c(4, 7))

      x_range <- seq(1, length(outq[[2]][, j][outq[[3]][j, ]]), length.out = 5)
      x_original <- seq(0, 1, length.out = 5) * (max_vals[j] - min_vals[j]) + min_vals[j]
      axis(1, at = x_range, labels = round(x_original, 2))
      
      y_range <- c(4.5, 5, 5.5, 6, 6.5)
      axis(2, at = y_range, labels = y_range) 
      
      # Add upper and lower quantiles
      # lines(outq1[[2]][, j][outq1[[3]][j, ]] + mean(Combined_dat[, 11]), lwd = 1, col = "purple")
      # lines(outq9[[2]][, j][outq9[[3]][j, ]] + mean(Combined_dat[, 11]), lwd = 1, col = "orange")
      
    } else {
      par(mfg = c((j - 1) %/% 3 + 1, (j - 1) %% 3 + 1))
      
      lines(outq[[2]][, j][outq[[3]][j, ]] + mean(Combined_dat[, 11]), lwd = 1)
      
      # Add upper and lower quantiles
      # lines(outq1[[2]][, j][outq1[[3]][j, ]] + mean(Combined_dat[, 11]), lwd = 1, col = "purple")
      # lines(outq9[[2]][, j][outq9[[3]][j, ]] + mean(Combined_dat[, 11]), lwd = 1, col = "orange")
    }
  }
  first_plot <- FALSE
}

dev.off()


################
# Shading Plot # 
################

png(filename = "shade_plot.png", width = 1000, height = 1000)
par(mfrow = c(3, 3))
first_plot <- TRUE

for (i in 1:10) {
  
  test_indices <- which(splits == i)
  test_data <- Scaled_Centered[test_indices, ]
  train_data <- Scaled_Centered[-test_indices, ]
  n = nrow(train_data)
  
  # Fit the quantile models
  outq <- fit_qatf(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.5, 0, lambda5[i])
  outq1 <- fit_qatf(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.1, 0, lambda1[i])
  outq9 <- fit_qatf(t(train_data[, 2:10]), train_data[, 11], n, 9, 0.9, 0, lambda9[i])
  
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
      plot(NULL, xaxt = 'n', yaxt = 'n', main = plotnames[j],
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




#################
#  Single Plot  # 
#################

png(filename = "overall_plot.png", width = 1100, height = 1000)

par(mfrow = c(3, 3), mar = c(5, 6, 2, 2))
# Plots
n = nrow(Scaled_Centered)
# res5 <- qatf_cv(t(Scaled_Centered[, 2:10]), Scaled_Centered[, 11], n, 9, 0.5, 0)
# lambda5 <- res5$LAMBDAS[which.min(rev(res5$MEANS))]
# outq <- fit_qatf(t(Scaled_Centered[, 2:10]), Scaled_Centered[, 11], n, 9, 0.5, 0, lambda5)
# 
# res1 <- qatf_cv(t(Scaled_Centered[, 2:10]), Scaled_Centered[, 11], n, 9, 0.1, 0)
# lambda1 <- res1$LAMBDAS[which.min(rev(res1$MEANS))]
# outq1 <- fit_qatf(t(Scaled_Centered[, 2:10]), Scaled_Centered[, 11], n, 9, 0.1, 0, lambda1)
# 
# res9 <- qatf_cv(t(Scaled_Centered[, 2:10]), Scaled_Centered[, 11], n, 9, 0.9, 0)
# lambda9 <- res9$LAMBDAS[which.min(rev(res9$MEANS))]
# outq9 <- fit_qatf(t(Scaled_Centered[, 2:10]), Scaled_Centered[, 11], n, 9, 0.9, 0, lambda9)
  
for (i in 1:9) {
  # Generate the plot
  plot(outq[[2]][, i][outq[[3]][i, ]] + mean(Combined_dat[, 11]),
       xaxt = 'n', yaxt = 'n', main = plotnames[i],
       xlab = "", lwd = 2, type = "l", ylim = c(3, 8), 
       ylab = (if (i %in% c(1, 4, 7)) "Happiness Index" else ""), 
       cex.lab = 1.5, cex.main = 1.5, font.lab = 2, font.main = 2)
    
  x_range <- seq(1, length(outq[[2]][, i][outq[[3]][i, ]]), length.out = 5)
  x_original <- seq(0, 1, length.out = 5) * (max_vals[i] - min_vals[i]) + min_vals[i]
  axis(1, at = x_range, labels = round(x_original, 2))
    
  y_range <- c(4, 5, 6, 7)
  axis(2, at = y_range, labels = y_range) 
    
  lines(outq1[[2]][, i][outq1[[3]][i, ]] + mean(Combined_dat[, 11]), lwd = 2, col = "purple")
  lines(outq9[[2]][, i][outq9[[3]][i, ]] + mean(Combined_dat[, 11]), lwd = 2, col = "orange")
}

dev.off()

