# QATF algorithms.R

# install.packages("devtools")
# library(devtools)
# install_github("zzh237/detrendr")
# install_github("glmgen/genlasso")
# not sure if this next one is necessary - Zhi can you check? 
# install_github("statsmaths/glmgen", subdir="R_pkg/glmgen")
# install.packages("fields")
# install.packages("plotly")

library(detrendr)
# trace(get_model, edit=TRUE)
library(glmgen)
library(fields)
library(tidyverse)


MSE <- function(a, b){
  len = length(a)
  return(norm(a - b, type="2")**2/len) 
}

# These functions require all inputs to be passed in
# Return an object with the best_mse, best_lambda, and best_fit
# optional parameters  allow for better control
get_mse <- function(x, y, y_star, n, d, k, alpha = 10**-4, max_t = 50, prints = TRUE, plots = FALSE){
  lambda_list <- 10**seq(3, -7, length.out=50)
  if (plots) {
    lambda_res <- numeric(50)
    lambda_i <- 1
  }
  
  y_mean <- mean(y)
  
  best_mse <- Inf 
  best_lambda <- Inf
  best_trend_hat <- rep(0, times = n)
  
  ord <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  # Calculate order for each row of x, and sort each row of x
  # Doing ahead of time saves cost
  for (j in 1:nrow(x)) { 
    ord[j, ] <- order(x[j, ])
    x[j, ] <- x[j, ][ord[j, ]]
  }
  
  for (lambda in lambda_list) {
    zero_matrix <- as.data.frame(matrix(0, nrow = n, ncol = d)) 
    trend_list <- as.data.frame(zero_matrix) 
    # we initialize all component functions to be 0
    
    trend_hat_prev <- rep(0, times = n)
    t <- 1
    repeat {
      for (j in 1:d) {
        # calculate jth partial residual using components not equal to j
        if (d == 2) {resp <- y - y_mean - as.numeric(trend_list[, -j])}
        else {resp <- y - y_mean - as.numeric(t(rowSums(trend_list[, -j])))}
        
        fit <- as.numeric(trendfilter(x[j, ], resp[ord[j, ]], k = k, lambda = lambda, thinning = FALSE)$beta)
        
        # then unorder the fit 
        trend_list[, j] <- fit[order(ord[j, ])] 
      }
      trend_hat <- rowSums(trend_list)
      
      if (all((abs(trend_hat - trend_hat_prev) <= alpha))) {break}
      else if (t >= max_t) {break}
      
      trend_hat_prev <- trend_hat
      t <- t + 1
    }
    
    current_mse <- MSE(y_star, trend_hat + y_mean) 
    if (prints) {cat("lambda of ", lambda, " achieved true MSE of ", current_mse, "\n")}
    if (plots) {
      lambda_res[lambda_i] <- current_mse
      lambda_i <- lambda_i + 1
    }
    
    if (current_mse < best_mse) {
      best_mse <- current_mse
      best_lambda <- lambda
      best_trend_hat <- trend_hat
      best_components <- trend_list
    }
  }
  if(plots) {
    # Basic plot with log scale on x-axis
    plot(lambda_list, lambda_res, log="x", 
         xlab="Lambda (log scale)", ylab="MSE", main="Ablation Plot ATF")
    
  }
  return(list("MSE" = best_mse, "LAMBDA" = best_lambda, "FIT" = best_trend_hat, "COMP" = best_components))
}
get_mse_s <- function(x, y, y_star, n, d, tau, alpha = 10**-4, max_t = 50, prints = TRUE, plots = FALSE){
  lambda_list <- 10**seq(0, -16, length.out=50)
  if (plots) {
    lambda_res <- numeric(50)
    lambda_i <- 1
  }
  # cubic splines of this form seem to prefer very small lambda values
  
  best_mse <- Inf 
  best_lambda <- Inf
  best_s_trend_hat <- rep(0, times = n)
  
  ord <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  # Calculate order for each row of x, and sort each row of x
  # Doing ahead of time saves cost
  for (j in 1:nrow(x)) { 
    ord[j, ] <- order(x[j, ])
    x[j, ] <- x[j, ][ord[j, ]]
  }
  
  
  for (lambda in lambda_list) {
    s_trend_list <- as.data.frame(matrix(0, nrow = n, ncol = d)) 
    # we initialize all component functions to be 0
    
    s_trend_hat_prev <- rep(0, times = n)
    t <- 1
    repeat {
      for (j in 1:d){
        # calculate jth partial residual using components not equal to j
        if (d == 2) {resp <- y - as.numeric(s_trend_list[, -j])}
        else {resp <- y - as.numeric(t(rowSums(s_trend_list[, -j])))}
        
        # order inputs
        fit <- qsreg(x[j, ], resp[ord[j, ]], lam = lambda, alpha = tau)$fitted.values
        # unorder fit
        s_trend_list[, j] <- fit[order(ord[j, ])] 
      }
      s_trend_hat <- rowSums(s_trend_list)
      
      if (all((abs(s_trend_hat - s_trend_hat_prev) <= alpha))) {break}
      else if (t >= max_t) {break}
      
      s_trend_hat_prev <- s_trend_hat
      t <- t + 1
    }
    
    current_mse <- MSE(y_star, s_trend_hat) 
    if (prints) {cat("lambda of ", lambda, " achieved true MSE of ", current_mse, "\n")}
    if (plots) {
      lambda_res[lambda_i] <- current_mse
      lambda_i <- lambda_i + 1
    }
    
    if (current_mse < best_mse) {
      best_mse <- current_mse
      best_lambda <- lambda
      best_s_trend_hat <- s_trend_hat
      best_components <- s_trend_list
    }
  }
  if(plots) {
    # Basic plot with log scale on x-axis
    plot(lambda_list, lambda_res, log="x", 
         xlab="Lambda (log scale)", ylab="MSE", main="Ablation Plot QSS")
    
  }
  return(list("MSE" = best_mse, "LAMBDA" = best_lambda, "FIT" = best_s_trend_hat, "COMP" = best_components))
}
get_mse_q <- function(x, y, y_star, n, d, tau, k, alpha = 10**-4, max_t = 50, prints = TRUE, plots = FALSE){
  lambda_list <- 10**seq(3 + k, -3 + k, length.out=50)
  if (plots) {
    lambda_res <- numeric(50)
    lambda_i <- 1
  }
  
  best_mse <- Inf 
  best_lambda <- Inf
  best_q_trend_hat <- rep(0, times = n)
  
  ord <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  # Calculate order for each row of x, and sort each row of x
  # Doing ahead of time saves cost
  for (j in 1:nrow(x)) { 
    ord[j, ] <- order(x[j, ])
    x[j, ] <- x[j, ][ord[j, ]]
  }
  
  for (lambda in lambda_list) {
    q_trend_list <- as.data.frame(matrix(0, nrow = n, ncol = d)) 
    # we initialize all component functions to be 0
    
    q_trend_hat_prev <- rep(0, times = n)
    t <- 1
    repeat {
      for (j in 1:d){
        # calculate jth partial residual using components not equal to j
        if (d == 2) {resp <- y - as.numeric(q_trend_list[, -j])}
        else {resp <- y - as.numeric(t(rowSums(q_trend_list[, -j])))}
        
        # for some reason, get_trend's k is one above expected (e.g. 2 is linear fit)
        # fit on ordered response
        fit <- get_trend(resp[ord[j, ]], tau, lambda, k+1)
        # get_trend requires equally spaced points
        
        # unorder fit
        q_trend_list[, j] <- fit[order(ord[j, ])]
      }
      q_trend_hat <- rowSums(q_trend_list)
      
      if (all((abs(q_trend_hat - q_trend_hat_prev) <= alpha))) {break}
      else if (t >= max_t) {break}
      
      q_trend_hat_prev <- q_trend_hat
      t <- t + 1
    }
    
    current_mse <- MSE(y_star, q_trend_hat) 
    if (prints) {cat("lambda of ", lambda, " achieved true MSE of ", current_mse, "\n")}
    if (plots) {
      lambda_res[lambda_i] <- current_mse
      lambda_i <- lambda_i + 1
    }
    
    if (current_mse < best_mse) {
      best_mse <- current_mse
      best_lambda <- lambda
      best_q_trend_hat <- q_trend_hat
      best_components <- q_trend_list
    }
  }
  if(plots) {
    # Basic plot with log scale on x-axis
    plot(lambda_list, lambda_res, log="x", 
         xlab="Lambda (log scale)", ylab="MSE", main="Ablation Plot QATF")
    
  }
  return(list("MSE" = best_mse, "LAMBDA" = best_lambda, "FIT" = best_q_trend_hat, "COMP" = best_components))
}

# get_mse uses 50 lambdas in 10^5 to 10^-9
# get_mse uses 50 lambdas in 10^1 to 10^-14
# get_mse uses 50 lambdas in 10^5 to 10^-9

# These functions employ the algorithm for a single lambda, 
# and output the fit, coefficients for each component, permutation matrix of inputs.
# Designed for use in plotting after best lambda is acquired using get_mse functions. 
fit_atf <- function(x, y, n, d, k, lambda, alpha = 10**-4, max_t = 50) {
  ord <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  # Calculate order for each row of x, and sort each row of x
  # Doing ahead of time saves cost
  for (j in 1:nrow(x)) { 
    ord[j, ] <- order(x[j, ])
    x[j, ] <- x[j, ][ord[j, ]]
  }
  
  y_mean <- mean(y)
  zero_matrix <- as.data.frame(matrix(0, nrow = n, ncol = d)) 
  trend_list <- as.data.frame(zero_matrix) 
  # we initialize all component functions to be 0
  
  trend_hat_prev <- rep(0, times = n)
  t <- 1
  repeat {
    for (j in 1:d) {
      # calculate jth partial residual using components not equal to j
      if (d == 2) {resp <- y - y_mean - as.numeric(trend_list[, -j])}
      else {resp <- y - y_mean - as.numeric(t(rowSums(trend_list[, -j])))}
      
      fit <- as.numeric(trendfilter(x[j, ], resp[ord[j, ]], k = k, lambda = lambda, thinning = FALSE)$beta)
      
      # then unorder the fit 
      trend_list[, j] <- fit[order(ord[j, ])] 
    }
    trend_hat <- rowSums(trend_list)
    
    if (all((abs(trend_hat - trend_hat_prev) <= alpha))) {break}
    else if (t >= max_t) {break}
    
    trend_hat_prev <- trend_hat
    t <- t + 1
  }
  
  return(list("fit" = trend_hat, "components" = trend_list, "order" = ord))
}
fit_qass <- function(x, y, n, d, tau, lambda, alpha = 10**-4, max_t = 50) {
  ord <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  # Calculate order for each row of x, and sort each row of x
  # Doing ahead of time saves cost
  for (j in 1:nrow(x)) { 
    ord[j, ] <- order(x[j, ])
    x[j, ] <- x[j, ][ord[j, ]]
  }
  
  s_trend_list <- as.data.frame(matrix(0, nrow = n, ncol = d)) 
  # we initialize all component functions to be 0
  
  s_trend_hat_prev <- rep(0, times = n)
  t <- 1
  repeat {
    for (j in 1:d){
      # calculate jth partial residual using components not equal to j
      if (d == 2) {resp <- y - as.numeric(s_trend_list[, -j])}
      else {resp <- y - as.numeric(t(rowSums(s_trend_list[, -j])))}
      
      # order inputs
      fit <- qsreg(x[j, ], resp[ord[j, ]], lam = lambda, alpha = tau)$fitted.values
      # unorder fit
      s_trend_list[, j] <- fit[order(ord[j, ])] 
    }
    s_trend_hat <- rowSums(s_trend_list)
    
    if (all((abs(s_trend_hat - s_trend_hat_prev) <= alpha))) {break}
    else if (t >= max_t) {break}
    
    s_trend_hat_prev <- s_trend_hat
    t <- t + 1
  }
  
  return(list("fit" = s_trend_hat, "components" = s_trend_list, "order" = ord))
}
fit_qatf <- function(x, y, n, d, tau, k, lambda, alpha = 10**-4, max_t = 50) {
  ord <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  # Calculate order for each row of x, and sort each row of x
  # Doing ahead of time saves cost
  for (j in 1:nrow(x)) { 
    ord[j, ] <- order(x[j, ])
    x[j, ] <- x[j, ][ord[j, ]]
  }
  
  q_trend_list <- as.data.frame(matrix(0, nrow = n, ncol = d)) 
  # we initialize all component functions to be 0
  
  q_trend_hat_prev <- rep(0, times = n)
  t <- 1
  repeat {
    for (j in 1:d){
      # calculate jth partial residual using components not equal to j
      if (d == 2) {resp <- y - as.numeric(q_trend_list[, -j])}
      else {resp <- y - as.numeric(t(rowSums(q_trend_list[, -j])))}
      
      # for some reason, get_trend's k is one above expected (e.g. 2 is linear fit)
      # fit on ordered response
      fit <- get_trend(resp[ord[j, ]], tau, lambda, k+1)
      # get_trend requires equally spaced points
      
      # unorder fit
      q_trend_list[, j] <- fit[order(ord[j, ])]
    }
    q_trend_hat <- rowSums(q_trend_list)
    
    if (all((abs(q_trend_hat - q_trend_hat_prev) <= alpha))) {break}
    else if (t >= max_t) {break}
    
    q_trend_hat_prev <- q_trend_hat
    t <- t + 1
  }
  
  return(list("fit" = q_trend_hat, "components" = q_trend_list, "order" = ord))
}

# TODO:
# Add fit function / combine fit function to do CV 
# So we can do feature selection for real data 




