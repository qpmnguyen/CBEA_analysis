library(tidyverse)
library(VGAM)
library(SpiecEasi)
library(patchwork)

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

result <- t(result)
result <- as.data.frame(result)


p1 <- qplot(result$size, geom = "histogram", xlim = c(0,1), fill = I("steelblue"), col = I("black")) + 
  geom_vline(aes(xintercept = median(result$size)), size = 1.5, col = "red") + theme_bw() + 
  labs(x = "Size", y = "Frequency") + 
  annotate("text", x = 0.3, y = 1000, label = "Median = 0.25", hjust = 0)

p2 <- qplot(result$munb, geom = "histogram", xlim = c(0,5), fill = I("steelblue"), col = I("black")) + 
  geom_vline(aes(xintercept = median(result$munb)), size = 1.5, col = "red") + theme_bw() + 
  labs(x = "Mean", y = "Frequency") + 
  annotate("text", x = 0.72, y = 1500, label = "Median = 0.66", hjust = 0)

p3 <- qplot(result$pstr0, geom = "histogram", xlim = c(0,1), fill = I("steelblue"), col = I("black")) + 
  geom_vline(aes(xintercept = median(result$pstr0)), size = 1.5, col = "red") + theme_bw() + 
  labs(x = "pstr0", y = "Frequency") + 
  annotate("text", x = 0.35, y = 1000, label = "Median = 0.31", hjust = 0)

fitting_plot <- (p1 / p2 / p3) + plot_annotation(tag_levels = "A")
ggsave(fitting_plot, filename = "docs/manuscript/figures/HMP_fit.png", dpi = 300, width = 10, height = 10)
