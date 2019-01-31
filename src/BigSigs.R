
library(ggplot2)
library(readr)
library(stringr)
require(randomForest)
library(doParallel)

registerDoParallel(4)


load("data/big_rec_dat.rda")
load('results/brain_rf_results_v2.rda')
load('results/brain_lasso_results_v2.rda')
load('data/cellids.rda')

# what we need #

# gene level for each cell type.

# subset the main dat matrix to the predictive sigs

# 
cells <- c('Astro_FGFR3',
           'OPC_PDGFRA',
           'Oligo_OPALIN',
           'Inh_PVALB',
           'Inh_VIP',
           'Inh_LAMP5',
           'Inh_SST',
           'Exc_FEZF2',
           'Exc_RORB',
           'Exc_THEMIS',
           'Micro_TYROBP',
           'Exc_LINC00507')

# Getting a matrix of gene expression by cell type.

idx <- cset13 %in% cells
celltype <- cset13[cset13 %in% cells]

datPrime <- dat[cset13 %in% cells,]

celPrime <- cset13[cset13 %in% cells]

datPrime <- cbind(data.frame(Cell=as.character(celPrime), stringsAsFactors = F), as.data.frame(datPrime))

datSum <- datPrime %>% group_by(Cell) %>% summarise_all(funs(sum))


# a cell signature for each cell type.
#   that remains predictive

impGenes <- list()
for (ci in datSum$Cell) {
  impTable <- rfList[[ci]]$importance
  impTable <- impTable[order(impTable[,4], decreasing = T),]
  theseGenes <- as.numeric(impTable[,4])
  names(theseGenes) <- rownames(impTable)
  impGenes[[ci]] <- theseGenes
}

# smallest number of genes for a cell type 
m <- min(unlist(lapply(impGenes, function(a) length(a))))
# 26

cellSigMat <- data.frame()
kappas <- c()
for (i in 1:m) {
  print(i)
  # for each cell type, get a gene and add the row  
  for (ci in datSum$Cell) {
    gi <- names(impGenes[[ci]])[i]
    rowi <- data.frame(t(datSum[,gi]))
    cellSigMat <- rbind(cellSigMat, (rowi/sum(rowi)))  
    if (nrow(cellSigMat) > 10) {
      kappas <- c(kappas, kappa(cellSigMat))
    }
  }
}

plot(kappas, type='l')

save(kappas, cellSigMat, impGenes, file='BigSig_Results.rda')

#

# then deconvolve the expression data.

# then compute gene levels for each cell type, 
# as long as we have gene expression measured for them.



