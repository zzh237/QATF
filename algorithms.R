# QATF algorithms.R

library(detrendr)
# trace(get_model, edit=TRUE)
# needs gurobi license: https://portal.gurobi.com/iam/licenses/list/
library(glmgen)
library(fields)
library(tidyverse)
library(FNN)
library(e1071)


MSE <- function(a, b){
  len = length(a)
  return(norm(a - b, type="2")**2/len) 
}

check <- function(u, tau) {u*(tau - 1*ifelse(u < 0, 1, 0))}

# These functions require all inputs to be passed in
# Return a range of MSEs and range of Lambdas
# optional parameters  allow for better control
get_mse <- function(x, y, y_star, n, d, k, alpha = 10**-4, max_t = 50, prints = TRUE, plots = FALSE){
  lambda_list <- 10**seq(3, -7, length.out=50)
  lambda_res <- numeric(50)
  
  y_mean <- mean(y)
  
  ord <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  # Calculate order for each row of x, and sort each row of x
  # Doing ahead of time saves cost
  for (j in 1:nrow(x)) { 
    ord[j, ] <- order(x[j, ])
    x[j, ] <- x[j, ][ord[j, ]]
  }
  
  lambda_i <- 1
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
    
    lambda_res[lambda_i] <- MSE(y_star, trend_hat + y_mean) 
    
    if (prints) {cat("lambda of ", lambda, " achieved true MSE of ", lambda_res[lambda_i], "\n")}
    lambda_i <- lambda_i + 1
  }
  if(plots) {
    # Basic plot with log scale on x-axis
    plot(lambda_list, lambda_res, log="x", 
         xlab="Lambda (log scale)", ylab="MSE", main="Ablation Plot ATF")
    
  }
  return(list("MSES" = lambda_res, "LAMBDAS" = lambda_list))
}
get_mse_s <- function(x, y, y_star, n, d, tau, alpha = 10**-4, max_t = 50, prints = TRUE, plots = FALSE){
  lambda_list <- 10**seq(0, -16, length.out=50)
  lambda_res <- numeric(50)
  # cubic splines of this form seem to prefer very small lambda values
  
  ord <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  # Calculate order for each row of x, and sort each row of x
  # Doing ahead of time saves cost
  for (j in 1:nrow(x)) { 
    ord[j, ] <- order(x[j, ])
    x[j, ] <- x[j, ][ord[j, ]]
  }
  
  lambda_i <- 1
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
    
    lambda_res[lambda_i] <- MSE(y_star, s_trend_hat) 
    if (prints) {cat("lambda of ", lambda, " achieved true MSE of ", lambda_res[lambda_i], "\n")}
    
    lambda_i <- lambda_i + 1
  }
  if(plots) {
    # Basic plot with log scale on x-axis
    plot(lambda_list, lambda_res, log="x", 
         xlab="Lambda (log scale)", ylab="MSE", main="Ablation Plot QSS")
    
  }
  return(list("MSES" = lambda_res, "LAMBDAS" = lambda_list))
}
get_mse_q <- function(x, y, y_star, n, d, tau, k, alpha = 10**-4, max_t = 50, prints = TRUE, plots = FALSE){
  lambda_list <- 10**seq(4 + k, -3 + k, length.out=50)
  lambda_res <- numeric(50)
  
  ord <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  # Calculate order for each row of x, and sort each row of x
  # Doing ahead of time saves cost
  for (j in 1:nrow(x)) { 
    ord[j, ] <- order(x[j, ])
    x[j, ] <- x[j, ][ord[j, ]]
  }
  
  lambda_i <- 1
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
    
    lambda_res[lambda_i] <- MSE(y_star, q_trend_hat) 
    if (prints) {cat("lambda of ", lambda, " achieved true MSE of ", lambda_res[lambda_i], "\n")}
    lambda_i <- lambda_i + 1

  }
  if(plots) {
    # Basic plot with log scale on x-axis
    plot(lambda_list, lambda_res, log="x", 
         xlab="Lambda (log scale)", ylab="MSE", main="Ablation Plot QATF")
  }
  return(list("MSES" = lambda_res, "LAMBDAS" = lambda_list))
}

# get_mse uses 50 lambdas in 10^3 to 10^-7
# get_mse_s uses 50 lambdas in 10^0 to 10^-16
# get_mse uses 50 lambdas in 10^(4+k) to 10^(-3+k)

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


predict_fit <- function(data, fit, loc, method = "1-NN") {
  if (grepl("-NN$", method)) {
    k <- as.numeric(sub("(\\d+)-NN$", "\\1", method))
    return(knn.reg(train = data, test = loc, y = fit, k = k)$pred)
  } 
  else if (method == "SVM") {
    svm_model <- svm(data, fit)
    return(predict(svm_model, loc))
  }
  else if (method == "GBM") {stop("Not Implemented")}
  else if (method == "GPR") {stop("Not Implemented")}
  else {stop("Not Implemented")}
}
qatf_cv <- function(x, y, n, d, tau, k, alpha = 10**-4, nfolds = 5, max_t = 50, prints = FALSE){
  # Evaluates QATF for a range of lambdas via CV with the Check Function
  # Predictions are constructed via 1-NN
  # May choose to pass this in as an argument 
  folds <- sample(rep(1:nfolds, length.out = n))
  
  x_sep <- vector("list", nfolds)
  y_sep <- vector("list", nfolds)
  ord_sep <- vector("list", nfolds)

  for (fold in 1:nfolds) {
    train_idx <- which(folds != fold)
    x_sep[[fold]] <- x[, train_idx]
    y_sep[[fold]] <- y[train_idx]
    
    ord_sep[[fold]] <- matrix(NA, nrow = nrow(x_sep[[fold]]), ncol = ncol(x_sep[[fold]]))
    for (j in 1:nrow(x)) { 
      # breakdown of indexing here is easier to follow in other algorithms
      ord_sep[[fold]][j, ] <- order(x_sep[[fold]] [j, ])
      x_sep[[fold]][j, ] <- x_sep[[fold]][j, ][ord_sep[[fold]][j, ]]
    }
  }
  
  
  lambda_list <- 10**seq(3 + k, -2 + k, length.out=50)
  lambda_res <- numeric(50)
  
  
  lambda_i <- 1
  for (lambda in lambda_list) {
    cv_mean <- numeric(nfolds)
    
    for (fold in 1:nfolds) {
      fold_n <- ncol(x_sep[[fold]])
      q_trend_list <- as.data.frame(matrix(0, nrow = fold_n, ncol = d)) 
      
      q_trend_hat_prev <- rep(0, times = fold_n)
      t <- 1
      repeat {
        for (j in 1:d){
          # calculate jth partial residual using components not equal to j
          if (d == 2) {resp <- y_sep[[fold]] - as.numeric(q_trend_list[, -j])}
          else {resp <- y_sep[[fold]] - as.numeric(t(rowSums(q_trend_list[, -j])))}
          
          # for some reason, get_trend's k is one above expected (e.g. 2 is linear fit)
          # fit on ordered response
          fit <- get_trend(resp[ord_sep[[fold]][j, ]], tau, lambda, k+1)
          # get_trend requires equally spaced points
          
          # unorder fit
          q_trend_list[, j] <- fit[order(ord_sep[[fold]][j, ])]
        }
        q_trend_hat <- rowSums(q_trend_list)
        
        if (all((abs(q_trend_hat - q_trend_hat_prev) <= alpha))) {break}
        else if (t >= max_t) {break}
        
        q_trend_hat_prev <- q_trend_hat
        t <- t + 1
      }
      
      test_idx <- which(folds == fold)
      cv_mean[fold] <- mean(check(y[test_idx] - predict_fit(t(x_sep[[fold]]), q_trend_hat, t(x[, test_idx]), method = "1-NN"), tau))
      if (prints) {cat("\tPrediction mean check of fold ", fold, " is ", cv_mean[fold], " \n", sep = "")}
    }
    
    lambda_res[lambda_i] <- mean(cv_mean) 
    if (prints) {cat("lambda of ", lambda, " achieved ", nfolds, "-fold CV mean check of ",
                     lambda_res[lambda_i], " \n", sep = "")}
    lambda_i <- lambda_i + 1
    
  }
  return(list("MEANS" = lambda_res, "LAMBDAS" = lambda_list))
  
} 

qass_cv <- function(x, y, n, d, tau, alpha = 10**-4, nfolds = 5, max_t = 50, prints = FALSE){
  # Evaluates QASS for a range of lambdas via CV with the Check Function
  # Predictions are constructed via 1-NN
  # May choose to pass this in as an argument 
  folds <- sample(rep(1:nfolds, length.out = n))
  
  x_sep <- vector("list", nfolds)
  y_sep <- vector("list", nfolds)
  ord_sep <- vector("list", nfolds)
  
  for (fold in 1:nfolds) {
    train_idx <- which(folds != fold)
    x_sep[[fold]] <- x[, train_idx]
    y_sep[[fold]] <- y[train_idx]
    
    ord_sep[[fold]] <- matrix(NA, nrow = nrow(x_sep[[fold]]), ncol = ncol(x_sep[[fold]]))
    for (j in 1:nrow(x)) { 
      # breakdown of indexing here is easier to follow in other algorithms
      ord_sep[[fold]][j, ] <- order(x_sep[[fold]] [j, ])
      x_sep[[fold]][j, ] <- x_sep[[fold]][j, ][ord_sep[[fold]][j, ]]
    }
  }
  
  
  lambda_list <- 10**seq(0, -16, length.out=50)
  lambda_res <- numeric(50)
  
  
  lambda_i <- 1
  for (lambda in lambda_list) {
    cv_mean <- numeric(nfolds)
    
    for (fold in 1:nfolds) {
      fold_n <- ncol(x_sep[[fold]])
      s_trend_list <- as.data.frame(matrix(0, nrow = fold_n, ncol = d)) 
      
      s_trend_hat_prev <- rep(0, times = fold_n)
      t <- 1
      repeat {
        for (j in 1:d){
          # calculate jth partial residual using components not equal to j
          if (d == 2) {resp <-  y_sep[[fold]] - as.numeric(s_trend_list[, -j])}
          else {resp <-  y_sep[[fold]] - as.numeric(t(rowSums(s_trend_list[, -j])))}
          
          # order inputs
          fit <- qsreg(x_sep[[fold]][j, ], resp[ord_sep[[fold]][j, ]], lam = lambda, alpha = tau)$fitted.values
          # unorder fit
          s_trend_list[, j] <- fit[order(ord_sep[[fold]][j, ])]
        }
        
        s_trend_hat <- rowSums(s_trend_list)
        
        if (all((abs(s_trend_hat - s_trend_hat_prev) <= alpha))) {break}
        else if (t >= max_t) {break}
        
        s_trend_hat_prev <- s_trend_hat
        t <- t + 1
      }
      
      test_idx <- which(folds == fold)
      cv_mean[fold] <- mean(check(y[test_idx] - predict_fit(t(x_sep[[fold]]), s_trend_hat, t(x[, test_idx]), method = "1-NN"), tau))
      if (prints) {cat("\tPrediction mean check of fold ", fold, " is ", cv_mean[fold], " \n", sep = "")}
    }
    
    lambda_res[lambda_i] <- mean(cv_mean) 
    if (prints) {cat("lambda of ", lambda, " achieved ", nfolds, "-fold CV mean check of ",
                     lambda_res[lambda_i], " \n", sep = "")}
    lambda_i <- lambda_i + 1
    
  }
  return(list("MEANS" = lambda_res, "LAMBDAS" = lambda_list))
  
} 




  


