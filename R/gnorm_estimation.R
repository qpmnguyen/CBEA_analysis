# Fitting generalized normal using ML approach 
fit_gen_gauss <- function(data, impute=1e-5,...){
  data <- as.vector(scale(data))
  start <- mean(data)/sqrt(var(data))
  L <- length(data)
  mu <- median(data)
  beta <- newtonRaphson(fun = g, x0 = start, data = data, impute=impute, ...)$root
  alpha <- ((beta/L)*sum(abs(data - mu)^beta))^(1/beta)
  return(list(mu = mu, beta = beta, alpha = alpha))
}


param <- fit_gen_gauss(data = test)
test2 <- rgnorm(1000, mu = param$mu, alpha = param$alpha, beta = param$beta)
qqplot(x = test, y = test2)


test <- as.vector(parameters$scores[[1]])
library(pracma)
g <- function(x, data, impute = 1e-5){
  L <- length(data)
  mu <- median(data)
  data[data == 0] <- 1e-5 # small-value imputation
  term1 <- 1 + (digamma(1/x)/x)
  term2 <- (sum(abs(data - mu)^x * log(abs(data - mu)))/sum(abs(data-mu)^x)) 
  term3 <- log((x/L)* sum(abs(data - mu)^x))/x
  term1 - term2 + term3
}
start <- mean(test)/sqrt(var(test))
pracma::newtonRaphson(fun = g, x0 = start, data = test)