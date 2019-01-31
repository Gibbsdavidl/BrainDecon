
# getting importance for genes #

library(ggplot2)
library(readr)
library(stringr)
require(randomForest)
library(doParallel)

registerDoParallel(4)


load("data/big_rec_dat.rda")
load(file='results/brain_lasso_results_v2.rda')


#cluster
#1    Inh L4-6 SST B3GAT2
#2    Exc L5-6 RORB TTC12
#3    Exc L5-6 FEZF2 ABO

# have cellcluster as label for each cell type #

# have cellcluster, which is cell type and location level

# can also have cell with marker gene
clusterlabels <- strsplit(x = mtgsamp$cluster, ' ')
cset13 <- unlist( lapply(clusterlabels, function(a) paste(a[1], a[3], sep='_') ) )
table(cset13)

rfList <- list()

# for a single cell type and for a number of iterations 
dim(dat)
#[1] 15928 14340

for (ci in names(resList)[c(14,15,17,18)]) {
  
  # get the cell labels
  idcode <- ifelse(cset13 == ci, 1, 0)
  
  # get the genes of interest here
  genes <- resList[[ci]]$names
  
  # subset those genes
  datSub <- dat[,genes]
  
  if (ncol(datSub) > 1) {
    # fit random forest
    rffit <- randomForest(y = as.factor(idcode), x=datSub, ntree = length(genes) * 10, importance = T)
    
    #    put top variables into a list
    rfList[[ci]] <- rffit
  }
  else {
    rfList[[ci]] <- "no genes!"
  }
  
}

save(rfList, file='brain_rf_results_v2.rda')


