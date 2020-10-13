library(tidyverse)
library(VGAM)
library(edgeR)
library(patchwork)
library(HMP16SData)
library(phyloseq)

data <- V35() %>% subset(select = HMP_BODY_SUBSITE == "Stool" & VISITNO == 1) %>% as_phyloseq()
data <- subset_samples(data,!duplicated(RSID)) %>% 
  filter_taxa(function(x) (sum(x == 0)/length(x)) < 0.9, TRUE)
data <- prune_samples(sample_sums(data) >= 1000, data)


# Fitting negative binomial distribution using robust measures in edgeR
counts <- data@otu_table@.Data
design <- model.matrix(~1, data = data.frame(data@sam_data))
normFacts <- edgeR::calcNormFactors(counts, method = "TMM")
dge <- DGEList(counts = counts)
dge$samples$norm.factors <- normFacts
disp <- estimateDisp(y = dge, design = design, tagwise = TRUE)
fit <- glmFit(dge$counts,dispersion = disp$tagwise.dispersion, design = design)
fit$fitted.values

mu <- rowMeans(fit$fitted.values)
size <- 1/disp$tagwise.dispersion
pstr0 <- rowMeans((1+fit$fitted.values*disp$tagwise.dispersion)^(-1/disp$tagwise.dispersion))
result <- as.data.frame(cbind(mu,size,pstr0))
result$log_mu <- log1p(mu)


data_param <- list(log_mu = log1p(rowMeans(data@otu_table@.Data)), 
                   pstr0 = rowMeans(data@otu_table@.Data == 0))


p1 <- qplot(y = result$log_mu, x = result$log_mu)


cor(data_param$log_mu, result$log_mu)^2
cor(data_param$pstr0, result$pstr0)^2


# Plotting 
p1 <- qplot(result$size, geom = "histogram", fill = I("steelblue"), col = I("black")) + 
  geom_vline(aes(xintercept = median(result$size)), size = 1.5, col = "red") + theme_bw() + 
  labs(x = "Size", y = "Frequency") + 
  annotate("text", x = 0.2, y = 1000, label = glue("Median = {med}", med = round(median(result$size),2)), hjust = 0)

p2 <- qplot(result$mu, geom = "histogram", xlim = c(0,5), fill = I("steelblue"), col = I("black")) + 
  geom_vline(aes(xintercept = median(result$mu)), size = 1.5, col = "red") + theme_bw() + 
  labs(x = "Mean", y = "Frequency") + 
  annotate("text", x = 0.72, y = 1500, label = glue("Median = {med}", med = round(median(result$mu),2)), hjust = 0)

p3 <- qplot(result$pstr0, geom = "histogram", xlim = c(0,1), fill = I("steelblue"), col = I("black")) + 
  geom_vline(aes(xintercept = median(result$pstr0)), size = 1.5, col = "red") + theme_bw() + 
  labs(x = "pstr0", y = "Frequency") + 
  annotate("text", x = 0.83, y = 1000, label = glue("Median = {med}", med = round(median(result$pstr0),2)), hjust = 0)

fitting_plot <- (p1 / p2 / p3) + plot_annotation(tag_levels = "A")
fitting_plot
ggsave(fitting_plot, filename = "docs/manuscript/figures/HMP_fit.png", dpi = 300, width = 10, height = 10)
saveRDS(result, file = "R/fitted_values.rds")


