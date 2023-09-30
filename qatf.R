# install.packages("devtools")
# library(devtools)
# install_github("zzh237/detrendr")
# install_github("glmgen/genlasso")
# library(detrendr)
# library(genlasso)
trace(get_model, edit=TRUE)
MSE <- function(a, b){
  len = length(a)
  return(norm(a - b, type="2")**2/len) 
}

# Professor suggested possibly increasing this to 6
lambda_list <- seq(1, 4.5, length.out=300)
10**(lambda_list)


lambda <- 5
<<<<<<< HEAD
tau <- c(0.1)
n <- 5000
=======
tau <- c(0.5)
n <- 1000
>>>>>>> 0ca003d1a4e61eb17fbb889c9db01a997580594e
d <- 10



## prepare the inputs X, and the true function
x <- seq(1, n, 1)

<<<<<<< HEAD
# get sceneria 1, 2, 3，6
=======
# get sceneria 1, 2, 3
>>>>>>> 0ca003d1a4e61eb17fbb889c9db01a997580594e
x <- x/n


# get scenerio 4
# half <- n/2
# x[1:half] <- 3*x[1:half]
# x <- replace(x, seq(half+1, n, 1), 3*n - 3*tail(x,half))
# x <- x/n

# get scenerio 5
# x<- cos(6*pi*x/n)

# get scenerio 6
half <- n/2
x[1:half] <- (0.25*(x[1:half]/n)**0.5 + 1.375)/3
x <- replace(x, seq(half+1, n, 1), (7*(tail(x,half)/n)**0.5-2)/3)



# get the simulated heterogeneously-smooth data
x_j <- x # we could add permutation for different j.


# scenerio 6
# get_g_j <- function(x_j, j){
#   g_0 <- (x_j + 0.1)*(j/10)
#   g_0_n <- norm(g_0, type="2")/sqrt(n)
#   a_j <- 1/g_0_n
#   g_0 <- a_j**(x_j + 0.1)*(j/10)
#   return(g_0)
# }

# scenerio 5
# get_g_j <- function(x_j, j){
#   g_0 <- (x_j + 0.1)*(j/10)
#   g_0_n <- norm(g_0, type="2")/sqrt(n)
#   a_j <- 1/g_0_n
#   g_0 <- a_j*(x_j + 0.1)*(j/10)
#   return(g_0)
# }


# scenerio 4
# get_g_j <- function(x_j, j){
#   g_0 <- (x_j + 0.1)**(j/10)
#   g_0_n <- norm(g_0, type="2")/sqrt(n)
#   a_j <- 1/g_0_n
#   g_0 <- a_j*(x_j + 0.1)**(j/10)
#   return(g_0)
# }

<<<<<<< HEAD
# # scenerio 1, 2, 3, Doppler-like
# get_g_j <- function(x_j, j){
#   g_0 <- sin(2*pi/(x_j + 0.1)**(j/10))
#   g_0_n <- norm(g_0, type="2")/sqrt(n)
#   a_j <- 1/g_0_n
#   g_0 <- a_j*sin(2*pi/(x_j + 0.1)**(j/10))
#   return(g_0)
# }
=======
# scenerio 1, 2, 3, Doppler-like
get_g_j <- function(x_j, j){
  g_0 <- sin(2*pi/(x_j + 0.1)**(j/10))
  g_0_n <- norm(g_0, type="2")/sqrt(n)
  a_j <- 1/g_0_n
  g_0 <- a_j*sin(2*pi/(x_j + 0.1)**(j/10))
  return(g_0)
}
>>>>>>> 0ca003d1a4e61eb17fbb889c9db01a997580594e





# create y star and y 
y_list <- data.frame()
for (j in 1:d){
  y_j <- get_g_j(x_j, j) 
  y_list <- rbind(y_list, y_j)
}
#get the y_star
y_list <- unname(y_list)
y_star <- colSums(y_list)
#get the y_i
sce1 <- rnorm(n, 0, 1)
# sce2 <- rcauchy(n, 0, 1)
# e <- x**0.5/n**0.5
# te <- rt(n, 2)
# sce3 <- e*te
# sce4 <- rt(n,3)
# sce5 <- rcauchy(n, 0, 1)
y <- y_star + sce1

#scenario 6
<<<<<<< HEAD
sce6 <-rt(n,2)
sce6_9q <- qt(0.9, 2)
y<-y_star+sce6
y_9q <- y_star+sce6_9q
plot(y_9q, type="l", col="black", ylab='0.9 quantile')
=======
# sce6 <-rt(n,2)
# sce6_9q <- qt(0.9, 2)
# y<-y_star*sce6
# y_9q <- y_star*sce6_9q
# plot(y_9q, type="l", col="black", ylab='0.9 quantile')
>>>>>>> 0ca003d1a4e61eb17fbb889c9db01a997580594e

y_mean <- mean(y)

get_mse <- function(k){
  trend_list <- data.frame() # local
  cum_trend_sofar <-0
  for (j in 1:d){
    resp <- y - cum_trend_sofar
    # since we are using trendfilter, rather than QTF, this is ATF
    fit_matrix <- trendfilter(resp, ord=k)$fit # I'm imagine here we are using all lambdas, pulled globally
    mses<-apply(fit_matrix,2,MSE,b=resp) # this and the following line take the best model
    trend_j <- fit_matrix[,which.min(mses)]
    trend_list<-  rbind(trend_list,trend_j)
    cum_trend_sofar <- cum_trend_sofar + trend_j # im not sure exactly how this accomplishes line (i)
  }
  
  ## get MSE
  trend_hat <- colSums(trend_list)
  MSE_trend_hat <- MSE(y_star, trend_hat)
  
  # using the quantile additive trend filtering
  q_trend_list <- data.frame()
  cum_trend_sofar <-0
  for (j in 1:d){
    resp <- y - cum_trend_sofar
    # same concepts as above, but now with QTF
    q_trend_j <- c(get_trend(resp, tau, lambda, k))
    q_trend_list<-  rbind(q_trend_list,q_trend_j)
    cum_trend_sofar <- cum_trend_sofar + q_trend_j
  }
  
  ## get MSE
  q_trend_hat <- colSums(q_trend_list)
  MSE_q_trend_hat <- MSE(y_star, q_trend_hat)
  return(list(MSE_trend_hat, MSE_q_trend_hat))
  
  # where does the MC sampling come in? 
}

# This looks like the same thing? 

# # using the additive trend filtering
# 
# trend_list <- data.frame()
# cum_trend_sofar <-0
# for (j in 1:d){
#   resp <- y - y_mean - cum_trend_sofar
#   fit_matrix <- trendfilter(resp, ord=k)$fit
#   mses<-apply(fit_matrix,2,MSE,b=resp)
#   trend_j <- fit_matrix[,which.min(mses)]
#   trend_list<-  rbind(trend_list,trend_j)
#   cum_trend_sofar <- cum_trend_sofar + trend_j
# }
# 
# ## get MSE
# trend_hat <- colSums(trend_list)
# MSE_trend_hat <- MSE(y_star, trend_hat)
# 
# # using the quantile additive trend filtering
# q_trend_list <- data.frame()
# cum_trend_sofar <-0
# for (j in 1:d){
#   resp <- y - y_mean - cum_trend_sofar
#   q_trend_j <- c(get_trend(resp, tau, lambda, k))
#   q_trend_list<-  rbind(q_trend_list,q_trend_j)
#   cum_trend_sofar <- cum_trend_sofar + q_trend_j
# }
# 
# ## get MSE
# q_trend_hat <- colSums(q_trend_list)
# MSE_q_trend_hat <- MSE(y_star, q_trend_hat)




mse1 = get_mse(k=1)
ATF1 = mse1[[1]]
QATF1 = mse1[[2]]
mse2 = get_mse(k=2)
ATF2 = mse2[[1]]
QATF2 = mse2[[2]]

data = data.frame(n = n, Scenario = 6,
                  tau=tau, QATF1 = QATF1, QATF2 = QATF2,
                  ATF1 =ATF1, ATF2 =ATF2)

# write data to a sample.csv file
write.table(data, file = "~/QATF/sample2.csv", append = TRUE, quote = FALSE,
            col.names = FALSE, row.names = FALSE) 



#> Using same lambda for all quantiles
plot(y, type="l", col="black") # this is the y with errors
plot(y_star, type="l", col="black") # this is the y_star, without errors


# lines(trend_hat~x, col="red")
# lines(q_trend_hat~x, col="blue")
# lines(trend[,1]~x, col="red")
# lines(trend[,2]~x, col="blue")


