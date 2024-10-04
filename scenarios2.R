# QATF scenarios.R

library(plotly)
library(tidyverse)
library(extraDistr)

# Scenario functions wrap our scenarios
# construct plots withing
# and return the true quantile fit, y_star_q, and the data, y
# for tau = 0.5, y_star_q = y_star
scenario0 <- function(n, d, tau, plots = FALSE) {
  # Scenario 0
  # g_0 is piecewise constant function alternating between 1 and -1, 
  # at j + 2 breakpoints
  # x drawn randomly from uniform distribution for each component
  # f_0 <- a_j * piecewise_function - b_j with b_j such that mean(f_0) = 0 
  # and a_j such that norm(f_0) = 1
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  x_list <- matrix(NA, nrow = d, ncol = n)
  y_list <- matrix(NA, nrow = d, ncol = n)
  
  piecewise_constant <- function(x, breakpoints, values) {
    intervals <- findInterval(x, breakpoints)
    return(values[intervals])
  }
  
  for (j in 1:d) {
    x_list[j, ] <- runif(n, 0, 1)
    
    # Define breakpoints for piecewise constant function
    breakpoints <- seq(0, 1, length.out = j + 2)^2
    values <- rep(c(1, -1), length.out = j + 1)  # alternating between 1 and -1
    
    # Piecewise constant function
    g_0 <- piecewise_constant(x_list[j, ], breakpoints, values)
    b_j <- mean(g_0)
    a_j <- 1 / (norm(g_0 - b_j, type = "2") / sqrt(n))
    y_list[j, ] <- a_j * g_0 - a_j * b_j
    
    if (plots) {
      filename <- paste(getwd(), '/sce0_plot_', j, '.png', sep = '')
      png(filename)
      plot(y_list[j, ]~x_list[j, ], main = paste('Scenario 0 j =', j), xlab = 'x', ylab = 'y')
      dev.off()
    }
  }
  
  # Sum of each column of y_list
  y_star <- colSums(y_list)
  
  # T(3) errors
  y <- y_star + rt(n, 3)
  y_star_q <- y_star + qt(tau, 3)
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(x_list, y, y_star_q))
}
scenario1 <- function(n, d, tau, plots = FALSE) {
  # Scenario 1
  # g_0(x) <- sin(2*pi/(x + 0.1)**(j/10))
  # x <- 1/n equally spaced
  # f_0 <- a_j*g_0 - b_j w/ b_j s.t. mean(f_0) = 0 and a_j s.t. norm(f_0) = 1
  # y = f_0(x) + epsilon_i
  # epsilon_i normal errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  x_list <- matrix(NA, nrow = d, ncol = n)
  y_list <- matrix(NA, nrow = d, ncol = n)
  
  for (j in 1:d) {
    x_list[j, ] <- sample(seq(0, 1, length.out = n), replace = FALSE)
    
    # Doppler-like
    g_0 <- sin(2 * pi / (x_list[j, ] + 0.1)^(j / 10))
    b_j <- mean(g_0)
    a_j <- 1 / (norm(g_0 - b_j, type = "2") / sqrt(n))
    y_list[j, ] <- a_j * g_0 - a_j * b_j
    
    if (plots) {
      filename <- paste(getwd(), '/sce1_plot_', j, '.png', sep = '')
      png(filename)
      plot(y_list[j, ]~x_list[j, ], main = paste('Scenario 1 j =', j), xlab = 'x', ylab = 'y')
      dev.off()
    }
  }
  
  # Sum of each column of y_list
  y_star <- colSums(y_list)
  
  # Normal Errors
  y <- y_star + rnorm(n, 0, 1)
  y_star_q <- y_star + qnorm(tau, 0, 1)
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(x_list, y, y_star_q))
}
scenario2 <- function(n, d, tau, plots = FALSE) {
  # Scenario 2
  # g_0(x) <- sin(2*pi/(x + 0.1)**(j/10))
  # x drawn randomly from uniform distribution for each component
  # f_0 <- a_j*g_0 - b_j w/ b_j s.t. mean(f_0) = 0 and a_j s.t. norm(f_0) = 1
  # y = f_0(x) + epsilon_i
  # epsilon_i cauchy errors
  
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  x_list <- matrix(NA, nrow = d, ncol = n)
  y_list <- matrix(NA, nrow = d, ncol = n)
  
  for (j in 1:d) {
    x_list[j, ] <- runif(n, 0, 1)
    
    # Doppler-like
    g_0 <- sin(2 * pi / (x_list[j, ] + 0.1)^(j / 10))
    b_j <- mean(g_0)
    a_j <- 1 / (norm(g_0 - b_j, type = "2") / sqrt(n))
    y_list[j, ] <- a_j * g_0 - a_j * b_j
    
    if (plots) {
      filename <- paste(getwd(), '/sce2_plot_', j, '.png', sep = '')
      png(filename)
      plot(y_list[j, ]~x_list[j, ], main = paste('Scenario 2 j =', j), xlab = 'x', ylab = 'y')
      dev.off()
    }
  }
  
  # Sum of each column of y_list
  y_star <- colSums(y_list)
  
  # Cauchy Errors
  y <- y_star + rcauchy(n, 0, 1)
  y_star_q <- y_star + qcauchy(tau, 0, 1)
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(x_list, y, y_star_q))
}
scenario3 <- function(n, d, tau, plots = FALSE) {
  # Scenario 3
  # g_0(x) <- sin(2*pi/(x + 0.1)**(j/10))
  # x drawn randomly from uniform distribution for each component
  # f_0 <- a_j*g_0 - b_j w/ b_j s.t. mean(f_0) = 0 and a_j s.t. norm(f_0) = 1
  # y = f_0(x) + epsilon_i
  # epsilon_i lognormal
  
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  x_list <- matrix(NA, nrow = d, ncol = n)
  y_list <- matrix(NA, nrow = d, ncol = n)
  
  for (j in 1:d) {
    x_list[j, ] <- runif(n, 0, 1)
    
    # Doppler-like
    g_0 <- sin(2 * pi / (x_list[j, ] + 0.1)^(j / 10))
    b_j <- mean(g_0)
    a_j <- 1 / (norm(g_0 - b_j, type = "2") / sqrt(n))
    y_list[j, ] <- a_j * g_0 - a_j * b_j
    
    if (plots) {
      filename <- paste(getwd(), '/sce3_plot_', j, '.png', sep = '')
      png(filename)
      plot(y_list[j, ]~x_list[j, ], main = paste('Scenario 3 j =', j), xlab = 'x', ylab = 'y')
      dev.off()
    }
  }
  
  # Sum of each column of y_list
  y_star <- colSums(y_list)
  
  # Cauchy Errors
  y <- y_star + rlnorm(n, 0, 1) - 1
  y_star_q <- y_star + qlnorm(tau, 0, 1) - 1
  # y_avg <- y_star + rlnorm(n, 0, 1) - exp(1/2)
  # y_star_avg <- y_star - exp(1/2)
  
  # if (tau != 0.5) { warning("Skewed Distribution, extra outputs.")}
  return(list(x_list, y, y_star_q))
}
scenario4 <- function(n, d, tau, plots = FALSE) {
  # Scenario 4
  # g_0(x) <- sin(2*pi/(x + 0.1)**(j/10))
  # x drawn randomly from uniform distribution for each component
  # f_0 <- a_j*g_0 - b_j w/ b_j s.t. mean(f_0) = 0 and a_j s.t. norm(f_0) = 1
  # y = f_0(x) + epsilon_i
  # epsilon_i heteroskedastic t(2) errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  x_list <- matrix(NA, nrow = d, ncol = n)
  y_list <- matrix(NA, nrow = d, ncol = n)
  
  for (j in 1:d) {
    x_list[j, ] <- runif(n, 0, 1)
    
    # Doppler-like
    g_0 <- sin(2 * pi / (x_list[j, ] + 0.1)^(j / 10))
    b_j <- mean(g_0)
    a_j <- 1 / (norm(g_0 - b_j, type = "2") / sqrt(n))
    y_list[j, ] <- a_j * g_0 - a_j * b_j
    
    if (plots) {
      filename <- paste(getwd(), '/sce4_plot_', j, '.png', sep = '')
      png(filename)
      plot(y_list[j, ]~x_list[j, ], main = paste('Scenario 4 j =', j), xlab = 'x', ylab = 'y')
      dev.off()
    }
  }
  
  # Sum of each column of y_list
  y_star <- colSums(y_list)
  
  # Heteroskedastic T errors
  # function of index
  e <- (seq(1, n, 1)**.5)/(n**.5)  # i**.5 / n**.5
  y <- y_star + e * rt(n, 2)
  y_star_q <- y_star + e * qt(tau, 2)
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(x_list, y, y_star_q))
}
scenario5 <- function(n, d=3, tau, plots = FALSE) {
  # Scenario 5
  # i <- 1:n
  # g_1(x) <- (cos(6*pi*x) + 0.1)
  # g_2(x) <- piecwise constant
  # g_3(x) <- exp(3*x)*sin(4*pi*x)
  # # x drawn randomly from uniform distribution for each component
  # f_0 <- a_j*g_0 - b_j w/ b_j s.t. mean(f_0) = 0 and a_j s.t. norm(f_0) = 1
  # y = f_0(x) + epsilon_i
  # epsilon_i t(3) errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  if (d != 3) {
    stop("Scenario 5 is specifically designed for d = 3")
  }
  
  x_list <- matrix(NA, nrow = d, ncol = n)
  y_list <- matrix(NA, nrow = d, ncol = n)
  
  if (!plots) {par(mfrow = c(1, 3))}
  for (j in 1:d) {
    x_list[j, ] <- sample(seq(0, 1, length.out = n), replace = FALSE)
    
    g_0 <- switch(j,
                  cos(6*pi*x_list[j, ]) + 0.1,
                  ifelse(x_list[j, ] < 0.5, 3*x_list[j, ], 3*(1 - x_list[j, ])) + 0.1,
                  exp(3*x_list[j, ])*sin(4*pi*x_list[j, ]),
                  0*x_list[j, ])
    b_j <- mean(g_0)
    a_j <- 1 / (norm(g_0 - b_j, type="2")/sqrt(n))
    y_list[j, ] <- a_j*g_0 - a_j*b_j
    
    ord <- order(x_list[j, ])
    
    if (plots) {
      filename <- paste(getwd(), '/sce5_plot_', j, '.png', sep = '')
      png(filename)
      plot(y_list[j, ]~x_list[j, ], main = paste('Scenario 5 j =', j), xlab = 'x', ylab = 'y')
      dev.off()
    }
    else {
      plot(y_list[j, ][ord], col = "black", ylab = paste("Component ", j))
    }
  }
  y_star <- colSums(y_list)
  
  # t(3) errors
  y <- y_star + rt(n, 3)
  y_star_q <- y_star + qt(tau, 3)
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(x_list, y, y_star_q, y_list))
}
scenario6 <- function(n, d=4, tau, plots = FALSE) {
  # Scenario 6
  # i <- 1:n
  # g_1(x) <- −(t−1/2)**2
  # g_2(x) <- heterogenous sin wave
  # g_3(x) <- useless dimension
  # g_4(x) <- exp(3*x)*sin(4*pi*x)
  # # x drawn randomly from uniform distribution for each component
  # f_0 <- a_j*g_0 - b_j w/ b_j s.t. mean(f_0) = 0 and a_j s.t. norm(f_0) = 1
  # y = f_0(x) + epsilon_i
  # epsilon_i t(3) errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  if (d != 4) {
    stop("Scenario 6 is specifically designed for d = 4")
  }
  
  x_list <- matrix(NA, nrow = d, ncol = n)
  y_list <- matrix(NA, nrow = d, ncol = n)
  
  if (!plots) {par(mfrow = c(1, 4))}
  for (j in 1:d) {
    x_list[j, ] <- sample(seq(0, 1, length.out = n), replace = FALSE)
    
    g_0 <- switch(j,
                  -(x_list[j, ] - 1/2)**2,
                  1.5*sin(4*pi*x_list[j, ]) + ifelse(x_list[j, ] > 0.5, 0, sin(16*pi*x_list[j, ])),
                  sample(c(rep(1, ceiling(0.001*n)), rep(0.1, n-ceiling(0.001*n))), n),
                  exp(3*x_list[j, ])*sin(4*pi*x_list[j, ]),
                  0*x_list[j, ])
    b_j <- mean(g_0)
    a_j <- 1 / (norm(g_0 - b_j, type="2")/sqrt(n))
    y_list[j, ] <- a_j*g_0 - a_j*b_j
    
    ord <- order(x_list[j, ])
    
    if (plots) {
      filename <- paste(getwd(), '/sce6_plot_', j, '.png', sep = '')
      png(filename)
      plot(y_list[j, ]~x_list[j, ], main = paste('Scenario 6 j =', j), xlab = 'x', ylab = 'y')
      dev.off()
    } else {
      plot(y_list[j, ][ord], col = "black", ylab = paste("Component ", j))
    }
  }
  y_star <- colSums(y_list)
  
  # t(3) errors
  y <- y_star + rt(n, 3)
  y_star_q <- y_star + qt(tau, 3)
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(x_list, y, y_star_q, y_list))
}
scenario7 <- function(n, d=2, tau, plots = FALSE) {
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario. Please ensure inputs are each scalar values.")
  }
  if (d != 2) {
    stop("Scenario 7 is specifically designed for d = 2")
  }
  
  x_list <- matrix(NA, nrow = d, ncol = n)
  y_list <- matrix(NA, nrow = d, ncol = n)
  ord <- matrix(NA, nrow = d, ncol = n)
  
  for (j in 1:d) {
    x_list[j, ] <- sample(seq(0, 1, length.out = n), n, replace = FALSE)
    # x_list[j, ] <- seq(0, 1, length.out = n)
    
    g_0 <- switch(j,
                  -(x_list[j, ] - 1/2)**2,
                  1/2 * cos(6 * pi * x_list[j, ]) + 0.1)
    b_j <- mean(g_0)
    a_j <- 1 / (norm(g_0 - b_j, type = "2") / sqrt(n))
    y_list[j, ] <- a_j * (g_0 - b_j)
    
    ord[j, ] <- order(x_list[j, ])
    
    if (plots) {
      filename <- paste(getwd(), '/sce7_plot_', j, '.png', sep = '')
      png(filename)
      plot(y_list[j, ]~x_list[j, ], main = paste('Scenario 7 j =', j), xlab = 'x', ylab = 'y')
      dev.off()
    }
  }
  
  y_star <- colSums(y_list)
  
  plot <- plot_ly(x = ~x_list[1, ], y = ~x_list[2, ], z = ~y_star, 
                  type = "scatter3d", mode = "markers", 
                  marker = list(size = 3, color = ~y_star, colorscale = 'Viridis')) %>%
    layout(scene = list(
      xaxis = list(title = "x1"),
      yaxis = list(title = "x2"),
      zaxis = list(title = "y")
    ))
  
  print(plot)
  
  # t(3) errors
  y <- y_star + rt(n, 3)
  y_star_q <- y_star + qt(tau, 3)
  
  if (tau != 0.5) {
    warning("Tau != 0.5. Only use output for QATF!")
  }
  
  return(list(x_list, y, y_star_q, y_list))
}

scenario8 <- function(n, d, tau, plots = FALSE) {
  # Scenario 8
  # Piecewise constant function based on distance
  # x drawn randomly from uniform distribution for each component
  # epsilon_i cauchy errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario. Please ensure inputs are each scalar values.")
  }
  
  x_list <- matrix(NA, nrow = d, ncol = n)
  y_list <- matrix(NA, nrow = d, ncol = n)
  
  for (j in 1:d) {
    x_list[j, ] <- runif(n, 0, 1)
    
    # Piecewise constant based on distance
    g_0 <- ifelse(x_list[j, ] < 0.5, 1, -1)
    b_j <- mean(g_0)
    a_j <- 1 / (norm(g_0 - b_j, type = "2") / sqrt(n))
    y_list[j, ] <- a_j * g_0 - a_j * b_j
    
    if (plots) {
      filename <- paste(getwd(), '/sce8_plot_', j, '.png', sep = '')
      png(filename)
      plot(y_list[j, ]~x_list[j, ], main = paste('Scenario 8 j =', j), xlab = 'x', ylab = 'y')
      dev.off()
    }
  }
  
  y_star <- colSums(y_list)
  y <- y_star + rcauchy(n, 0, 1)
  y_star_q <- y_star + qcauchy(tau, 0, 1)
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!") }
  return(list(x_list, y, y_star_q))
}

scenario9 <- function(n, d, tau, plots = FALSE) {
  # Scenario 9
  # piecewise constant
  # x drawn randomly from uniform distribution for each component
  # f_0 <- a_j*piecewise_constant - b_j with b_j such that mean(f_0) = 0 
  # and a_j such that norm(f_0) = 1
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  x_list <- matrix(NA, nrow = d, ncol = n)
  y_list <- matrix(NA, nrow = d, ncol = n)
  
  piecewise_constant <- function(x) {
    return(ifelse(x < 0.5, 1, -1))
  }
  
  for (j in 1:d) {
    x_list[j, ] <- runif(n, 0, 1)
    g_0 <- piecewise_constant(x_list[j, ])
    b_j <- mean(g_0)
    a_j <- 1 / (norm(g_0 - b_j, type = "2") / sqrt(n))
    y_list[j, ] <- a_j * g_0 - a_j * b_j
    
    if (plots) {
      filename <- paste(getwd(), '/sce9_plot_', j, '.png', sep = '')
      png(filename)
      plot(y_list[j, ]~x_list[j, ], main = paste('Scenario 9 j =', j), xlab = 'x', ylab = 'y')
      dev.off()
    }
  }
  
  # Sum of each column of y_list
  y_star <- colSums(y_list)
  
  # Cauchy Errors
  y <- y_star + rcauchy(n, 0, 1)
  y_star_q <- y_star + qcauchy(tau, 0, 1)
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(x_list, y, y_star_q))
}

scenario10 <- function(n, d, tau, plots = FALSE) {
  # Scenario 10
  # g_0(x) is piecewise constant based on the distance to fixed points
  # x drawn randomly from uniform distribution for each component
  # f_0 <- a_j*g_0 - b_j with b_j such that mean(f_0) = 0 
  # and a_j such that norm(f_0) = 1
  # y = f_0(x) + epsilon_i
  # epsilon_i Cauchy errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  x_list <- matrix(NA, nrow = d, ncol = n)
  y_list <- matrix(NA, nrow = d, ncol = n)
  
  piecewise_constant <- function(x) {
    # Fixed points
    q1 <- rep(1/4, length(x))
    q2 <- rep(1/2, length(x))
    q3 <- rep(3/4, length(x))
    
    dist_q1 <- sqrt((x - q1)^2)
    dist_q2 <- sqrt((x - q2)^2)
    dist_q3 <- sqrt((x - q3)^2)
    
    min_dist <- pmin(dist_q1, dist_q2, dist_q3)
    values <- ifelse(min_dist == dist_q1, 1, ifelse(min_dist == dist_q2, -1, 0))
    return(values)
  }
  
  for (j in 1:d) {
    x_list[j, ] <- runif(n, 0, 1)
    g_0 <- piecewise_constant(x_list[j, ])
    b_j <- mean(g_0)
    a_j <- 1 / (norm(g_0 - b_j, type = "2") / sqrt(n))
    y_list[j, ] <- a_j * g_0 - a_j * b_j
    
    if (plots) {
      filename <- paste(getwd(), '/sce10_plot_', j, '.png', sep = '')
      png(filename)
      plot(y_list[j, ]~x_list[j, ], main = paste('Scenario 10 j =', j), xlab = 'x', ylab = 'y')
      dev.off()
    }
  }
  
  y_star <- colSums(y_list)
  y <- y_star + rcauchy(n, 0, 1)
  y_star_q <- y_star + qcauchy(tau, 0, 1)
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!")}
  return(list(x_list, y, y_star_q))
}


scenario11 <- function(n, d, tau, plots = FALSE) {
  # Scenario 11
  # g_0(x) is piecewise constant based on a complex pattern
  # x drawn randomly from uniform distribution for each component
  # f_0 <- a_j*g_0 - b_j with b_j such that mean(f_0) = 0 
  # and a_j such that norm(f_0) = 1
  # y = f_0(x) + epsilon_i
  # epsilon_i Cauchy errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  x_list <- matrix(NA, nrow = d, ncol = n)
  y_list <- matrix(NA, nrow = d, ncol = n)
  
  piecewise_constant <- function(x) {
    # Define a complex piecewise constant pattern
    q1 <- rep(1/3, length(x))
    q2 <- rep(2/3, length(x))
    q3 <- rep(1/2, length(x))
    
    dist_q1 <- sqrt((x - q1)^2)
    dist_q2 <- sqrt((x - q2)^2)
    dist_q3 <- sqrt((x - q3)^2)
    
    min_dist <- pmin(dist_q1, dist_q2, dist_q3)
    values <- ifelse(min_dist == dist_q1, 2, ifelse(min_dist == dist_q2, -2, 1))
    return(values)
  }
  
  for (j in 1:d) {
    x_list[j, ] <- runif(n, 0, 1)
    g_0 <- piecewise_constant(x_list[j, ])
    b_j <- mean(g_0)
    a_j <- 1 / (norm(g_0 - b_j, type = "2") / sqrt(n))
    y_list[j, ] <- a_j * g_0 - a_j * b_j
    
    if (plots) {
      filename <- paste(getwd(), '/sce11_plot_', j, '.png', sep = '')
      png(filename)
      plot(y_list[j, ] ~ x_list[j, ], main = paste('Scenario 11 j =', j), xlab = 'x', ylab = 'y')
      dev.off()
    }
  }
  
  y_star <- colSums(y_list)
  y <- y_star + rcauchy(n, 0, 1)
  y_star_q <- y_star + qcauchy(tau, 0, 1)
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!") }
  return(list(x_list, y, y_star_q))
}

scenario12 <- function(n, d, tau, plots = FALSE) {
  # Scenario 12
  # f_aux(x) = sqrt(sum(x)) + I(x1 < 0.5)
  # x drawn randomly from uniform distribution for each component
  # y = f_aux(x) + sqrt(x * beta) * t(2) errors
  # epsilon_i Cauchy errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario. Please ensure inputs are each scalar values.")
  }
  
  x_list <- matrix(runif(n * d, 0, 1), nrow = d, ncol = n)
  beta <- rep(0.5, d)
  beta[1] <- 1
  
  f_aux <- function(x) {
    sqrt(rowSums(x)) + as.numeric(x[, 1] < 0.5)
  }
  
  y_list <- apply(x_list, 2, function(x) {
    f_aux(matrix(x, nrow = 1)) + sqrt(sum(x * beta)) * rt(1, df = 2)
  })
  
  if (plots) {
    for (j in 1:d) {
      filename <- paste(getwd(), '/sce12_plot_', j, '.png', sep = '')
      png(filename)
      plot(y_list[j, ] ~ x_list[j, ], main = paste('Scenario 12 j =', j), xlab = 'x', ylab = 'y')
      dev.off()
    }
  }
  
  y_star <- y_list
  y <- y_star + rcauchy(n, 0, 1)
  y_star_q <- y_star + qcauchy(tau, 0, 1)
  
  if (tau != 0.5) {
    warning("Tau != 0.5. Only use output for QATF!")
  }
  
  return(list(x_list, y, y_star_q))
}

scenario13 <- function(n, d, tau, plots = FALSE) {
  # Scenario 13
  # g_0(x) is a complex nested function
  # x drawn randomly from uniform distribution for each component
  # f_0 <- g3(g2(g1(x)))
  # y = f_0(x) + epsilon_i
  # epsilon_i t(3) errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario. Please ensure inputs are each scalar values.")
  }
  
  g1 <- function(X) {
    return(cbind(sqrt(X[, 1]^2 + rowSums(X[, 2:d])), rowSums(X)^3))
  }
  
  g2 <- function(X) {
    return(cbind(abs(X[, 1]), apply(X, 1, prod)))
  }
  
  g3 <- function(X) {
    return(X[, 1] + sqrt(rowSums(X)))
  }
  
  x_list <- matrix(runif(n * d), nrow = d, ncol = n)
  y_list <- g3(g2(g1(t(x_list))))
  
  if (plots) {
    filename <- paste(getwd(), '/sce13_plot.png', sep = '')
    png(filename)
    plot(y_list, main = 'Scenario 13', xlab = 'x', ylab = 'y')
    dev.off()
  }
  
  y_star <- y_list
  y <- y_star + rcauchy(n, 0, 1)
  y_star_q <- y_star + qcauchy(tau, 0, 1)
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!") }
  return(list(x_list, y, y_star_q))
}

scenario14 <- function(n, d, tau, plots = FALSE) {
  # Scenario 14
  # g_0(x) is a sinusoidal function
  # x drawn randomly from uniform distribution for each component
  # f_0 <- sin(sum(cos(pi * x)))
  # y = f_0(x) + epsilon_i
  # epsilon_i normal errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario. Please ensure inputs are each scalar values.")
  }
  
  x_list <- matrix(runif(n * d), nrow = d, ncol = n)
  
  f_aux <- function(x) {
    return(sin(rowSums(cos(pi * t(x)))))
  }
  
  y_list <- f_aux(x_list)
  
  if (plots) {
    filename <- paste(getwd(), '/sce14_plot.png', sep = '')
    png(filename)
    plot(y_list, main = 'Scenario 14', xlab = 'x', ylab = 'y')
    dev.off()
  }
  
  y_star <- y_list
  y <- y_star + rcauchy(n, 0, 1)
  y_star_q <- y_star + qcauchy(tau, 0, 1)
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!") }
  return(list(x_list, y, y_star_q))
}

scenario15 <- function(n, d, tau, plots = FALSE) {
  # Scenario 15
  # g_0(x) is an exponential function based on distance from a fixed point
  # x drawn randomly from uniform distribution for each component
  # f_0 <- exp(-||x - c||^2)
  # y = f_0(x) + epsilon_i
  # epsilon_i normal errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario. Please ensure inputs are each scalar values.")
  }
  
  x_list <- matrix(runif(n * d), nrow = d, ncol = n)
  c <- rep(0.5, d)
  
  f_aux <- function(x) {
    return(exp(-rowSums((t(x) - c)^2)))
  }
  
  y_list <- f_aux(x_list)
  
  if (plots) {
    filename <- paste(getwd(), '/sce15_plot.png', sep = '')
    png(filename)
    plot(y_list, main = 'Scenario 15', xlab = 'x', ylab = 'y')
    dev.off()
  }
  
  y_star <- y_list
  y <- y_star + rcauchy(n, 0, 1)
  y_star_q <- y_star + qcauchy(tau, 0, 1)
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!") }
  return(list(x_list, y, y_star_q))
}


scenario16 <- function(n, d, tau, plots = FALSE) {
  # Scenario 16
  # Piece-wise constant function with different constants based on j
  # x drawn randomly from uniform distribution for each component
  # epsilon_i Cauchy errors
  
  if (length(n) != 1 || length(d) != 1 || length(tau) != 1) {
    stop("Scenario function is only suitable for a single scenario.\n
          Please ensure inputs are each scalar values.")
  }
  
  x_list <- matrix(NA, nrow = d, ncol = n)
  y_list <- matrix(NA, nrow = d, ncol = n)
  
  piecewise_function <- function(x_val, j) {
    if (x_val <= 0.2) {
      return(5 * j)
    } else if (x_val > 0.2 && x_val <= 0.4) {
      return(2 * j)
    } else if (x_val > 0.4 && x_val <= 0.6) {
      return(8 * j)
    } else {
      return(1 * j)
    }
  }
  
  for (j in 1:d) {
    x_list[j, ] <- runif(n, 0, 1)
    y_list[j, ] <- sapply(x_list[j, ], piecewise_function, j = j)
    
    if (plots) {
      filename <- paste(getwd(), '/sce16_plot_', j, '.png', sep = '')
      png(filename)
      plot(y_list[j, ] ~ x_list[j, ], main = paste('Scenario 16 j =', j), xlab = 'x', ylab = 'y')
      dev.off()
    }
  }
  
  # Sum of each column of y_list
  y_star <- colSums(y_list)
  
  # Cauchy Errors
  y <- y_star + rcauchy(n, 0, 1)
  y_star_q <- y_star + qcauchy(tau, 0, 1)
  
  if (tau != 0.5) { warning("Tau != 0.5. Only use output for QATF!") }
  return(list(x_list, y, y_star_q))
}

