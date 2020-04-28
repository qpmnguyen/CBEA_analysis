library(compositions)
library(phyloseq)
library(vegan)
library(factoextra)
library(fpc)
library(NbClust)
library(patchwork)

data(GlobalPatterns)

org <- GlobalPatterns %>% otu_table(taxa_are_rows = T) %>% t() %>% as("matrix") %>% 
  as.data.frame() %>% slice(1)

nvar <- ncol(org)

noise <- unclass(alrInv(rnorm(nvar-1)))

perturbed_data <- unclass(acomp(org))

for (i in 1:50){
  noise <- unclass(alrInv(rnorm(nvar-1)))
  noise <- perturbe(x = acomp(org), y = acomp(noise))
  perturbed_data <- rbind(perturbed_data, unclass(noise))
}
rownames(perturbed_data) <- c("org", paste0("rep",1:50))
covariates <- data.frame(annotation = c("org", rep("replicates",50)))
rownames(covariates) <- rownames(perturbed_data)


perturbed_phylo <- phyloseq(otu_table(perturbed_data, taxa_are_rows = F), 
                            tax_table(tax_table(GlobalPatterns)), 
                            sample_data(covariates))

genus <- tax_glom(perturbed_phylo, taxrank = "Genus")
family <- tax_glom(perturbed_phylo, taxrank = "Family")

generate_plots <- function(physeq, dist = "euclidean"){
  if (dist != 'bray'){
    physeq <- physeq %>% transform_sample_counts(function(x) unclass(clr(x)))
  }
  d <- distance(physeq, method = dist)
  ord <- ordinate(physeq, distance = d, method = "MDS")
  plt <- ggplot(cbind(as.data.frame(ord$vectors), covariates), aes(x = Axis.1, y = Axis.2, col = annotation)) + 
    geom_point() 
  return(list(plot = plt, ord = ord,dist = d))
}

genus_results <- generate_plots(genus, dist = "manhattan")
family_results <- generate_plots(family, dist = "manhattan")
otu_results <- generate_plots(perturbed_phylo, dist = "manhattan")


p1 <- otu_results$plot + theme_bw()+ geom_point(size = 3) + theme(legend.position = "none") + labs(title = "OTU level")
p2 <- genus_results$plot + theme_bw() + geom_point(size = 3) + theme(legend.position = "none") + labs(title = "Genus level")
p3 <- family_results$plot + theme_bw() + geom_point(size = 3) + labs(title = "Family level")
patch <- p1 + p2 + p3

proc_ana <- protest(genus_results$ord$vectors, family_results$ord$vectors)

ggplot() + geom_point(data = as.data.frame(proc_ana$X), aes(x = Axis.1, y = Axis.2, col = "Genus"), size = 3) + 
  geom_point(data = as.data.frame(proc_ana$Yrot), aes(x = V1, y = V2, col = "Family"), size = 3) + 
  scale_color_discrete(name = "Ordination") + theme_bw()
plot(proc_ana)

ggsave(patch, file = "ordination_aggregation_noise.png", dpi = 300, width = 8, height = 8)


png(filename = "docs/figures/procrustes_error_family_genus_simulation.png", 
    res = 300, units = "in", height = 8, width = 8)
plot(proc_ana, kind = 0, xlim = c(-0.1, 0.07), ylim = c(-0.05, 0.1))
points(proc_ana, display = "target", pch = 19, col = "red")
points(proc_ana, display = "rotated", pch = 19, col = "blue")
lines(proc_ana, type  = "segment")
legend(x = -0.1, y = 0.1, legend = c("Genus", "Family"), pch = c(19,19), col = c("red", "blue"))
dev.off()





