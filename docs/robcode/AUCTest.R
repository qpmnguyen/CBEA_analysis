source("./SimUtils.R")

args = commandArgs(T)
if (length(args) >= 1) {
  param.to.vary = args[1]
} else {
  param.to.vary = "sparsity"
}

# Simulate multiple data sets and test each one for a range of either sparsity values or SDs
num.values = 9
sparsity.values = rep(0.8, num.values)
sd.values = rep(.34, num.values)
mean.values = rep(1, num.values)    
size.values = rep(30, num.values)    
if (param.to.vary == "sparsity") {
  sparsity.values = seq(from=0.1, to=0.9, by=0.1)
} else if (param.to.vary == "sd") {
  sd.values = seq(from=0.15, to=0.55, by=0.05)
} else if (param.to.vary == "mean") {
  mean.values = seq(from=0.6, to=1.4, by=0.1)
} else if (param.to.vary == "size") {
  size.values = seq(from=10, to=50, by=5)
} 

save.results=T
all.auc.values = list()
all.exec.times = list()
num.sims = 50
null.mean = 0.642 # based on PBMC scRNA-seq data 
null.sd = 0.34 # based on PBMC scRNA-seq data 
method.names = c("VAM","GSVA", "ssGSEA", "PCA", "z-score")
num.methods = length(method.names)
exec.times = rep(0, num.methods)
names(exec.times) = method.names

for (i in 1:num.sims) {
  message("Processing simulation ", i)  

  auc.values = matrix(0, nrow=num.methods, ncol=num.values)
  exec.times = matrix(0, nrow=num.methods, ncol=num.values)  
  rownames(auc.values) = method.names
  rownames(exec.times) = method.names
    
  for (j in 1:num.values) {    
    sparsity = sparsity.values[j]
    null.sd = sd.values[j]
    inflate.mean = mean.values[j]
    set.size = size.values[j]
    message("Testing for sparsity ", sparsity, 
        ", null sd ", null.sd, 
        ", inflated mean ", inflate.mean,         
        ", set size ", set.size,             
        " for simulation ", i)
    sim.results = simulateData(
      #n = 1000, p = 100,
      #n.inflate = 50, p.inflate = 10,        
      n = 2000, p = 500,        
      n.inflate = 50, p.inflate = set.size,
      null.mean = null.mean, null.sd = null.sd, 
      inflate.mean = inflate.mean, inflate.sd=null.sd,
      sparsity=sparsity)
    aucs = computeAUCs(data=sim.results$data, inflate.data=sim.results$inflate.data,
        p.set=set.size, n.inflate=50)
    auc.values[1,j] = aucs$vam
    auc.values[2,j] = aucs$gsva
    auc.values[3,j] = aucs$ssgsea 
    auc.values[4,j] = aucs$plage
    auc.values[5,j] = aucs$zscore     
    
    exec.times[1,j] = aucs$vam.time    
    exec.times[2,j] = aucs$gsva.time
    exec.times[3,j] = aucs$ssgsea.time 
    exec.times[4,j] = aucs$plage.time
    exec.times[5,j] = aucs$zscore.time
  }
  
  all.auc.values[[i]] = auc.values
  all.exec.times[[i]] = exec.times
}

if (save.results) {
if (param.to.vary == "sparsity") {
  saveRDS(all.auc.values, "./sim_results/AUCTestSparsity.RDS")
  saveRDS(all.exec.times, "./sim_results/ExecTimesSparsity.RDS")  
} else if (param.to.vary == "sd") {
  saveRDS(all.auc.values, "./sim_results/AUCTestNoise.RDS")
  saveRDS(all.exec.times, "./sim_results/ExecTimesNoise.RDS")    
} else if (param.to.vary == "mean") {
  saveRDS(all.auc.values, "./sim_results/AUCTestEffectSize.RDS")
  saveRDS(all.exec.times, "./sim_results/ExecTimesEffectSize.RDS")  
} else if (param.to.vary == "size") {
  saveRDS(all.auc.values, "./sim_results/AUCTestSetSize.RDS")
  saveRDS(all.exec.times, "./sim_results/ExecTimesSetSize.RDS")    
}
}

#print(all.auc.values)
