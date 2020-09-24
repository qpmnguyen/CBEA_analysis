library(tidyverse)
library(VGAM)
library(SpiecEasi)

data <- readRDS(file = "data/stool_16S.rds")
abundant <- colSums(data)
idx <- order(abundant, decreasing = T)
data <- data[,idx]
data_eval <- data[,c(1:300)]

# fitting the zero inflated negative binomial distribution
result <- apply(data, 2, function(x){
  x <- as.numeric(x)
  res <- SpiecEasi::fitdistr(x, densfun = "zinegbin")$par
  return(res)
})

result <- do.call(rbind, result)
result <- as.data.frame(result)
result$par
qplot(result$size, geom = "histogram", xlim = c(0,1))
qplot(result$munb, geom = "histogram", xlim = c(0,10))
qplot(result$pstr0, geom = "histogram")

na.omit(result$pstr0)

min(result$mu, na.rm = T)
median(result$pstr0, na.rm = T)
median(result$size, na.rm = T)
median(result$munb, na.rm = T)

