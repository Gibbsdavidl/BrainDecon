
# run cybersort

library(ggplot2)
library(readr)
library(dplyr)

load('results/BigSig_Results.rda')
load("data/mayo.rda")

protGenes <- read_tsv("data/protein-coding_gene.txt")

source('src/pseudosort.R')

# from BigSigs.R
cellNames <- c("Astro_FGFR3", "Exc_FEZF2", "Exc_LINC00507", "Exc_RORB", "Exc_THEMIS", "Inh_LAMP5", "Inh_PVALB", "Inh_SST", "Inh_VIP", "Micro_TYROBP", "Oligo_OPALIN", "OPC_PDGFRA")

mayot <- t(mayo[,-1])
colnames(mayot) <- mayo$X1

res0 <- list()
for (i in 8:300) {
  
  csigs <- cellSigMat[1:i,]
  colnames(csigs) <- cellNames
  
  sigGeneNames <- protGenes[match(x = rownames(csigs), table = protGenes$entrez_id),] %>% select('symbol')
  
  idx <-  which(!is.na(sigGeneNames$symbol))
  csigs <- csigs[idx,]
  
  rownames(csigs) <- sigGeneNames$symbol[idx]
  
  # got ID translation from BioMart
  genetable <- read.table('data/mart_export.txt', sep='\t', stringsAsFactors=F, header=T)
  
  # then mapped the table to the mayo data
  mayosymbols <-unique(genetable[match(rownames(mayot), genetable$Gene.stable.ID),])
  
  # pulled out the intersection of selected genes from single cell with mayo genes
  intgenes <- intersect(rownames(csigs), mayosymbols$HGNC.symbol)
  intgenes <- intgenes[!is.na(intgenes)]
  
  # took a subset on the mapping table.
  mayoSymbolSub <- mayosymbols[mayosymbols$HGNC.symbol %in% intgenes,]
  
  # missing these: "IL1RAPL2" "ZNF804B"  "PRODH"   
  
  # getting the tables in the right order and with proper names
  B <- csigs[intgenes,]
  X <- mayot[mayoSymbolSub$Gene.stable.ID, ]
  rownames(X) <- mayoSymbolSub$HGNC.symbol
  X <- X[intgenes,]
  
  B[is.na(B)] <- 0
  
  # for each sample (rows in the mayo data)
  #   deconvolve using signatures.
  mayoDeconv <- lapply(1:ncol(X), function(k) {
    PSEUDOSORT(matrix(X[,k], ncol=1), as.matrix(B))
  })
  
  mayoDeconvMat <- do.call('rbind', mayoDeconv)
  rownames(mayoDeconvMat) <- as.character(mayo$X1)
  
  
  vals <- 1:264
  for (j in 1:264) {
    # each row of mayoDeconvMat is a sample
    # and B is the cell signatures
    # this produces a single vector of gene expression.
    check1 <- as.matrix(B, ncol = 12) %*% matrix(mayoDeconvMat[j,], ncol=1)
    #qplot(x=as.numeric(X[,i]), y=as.numeric(check1)) + geom_smooth(method='lm')
    vals[j] <- (cor(as.numeric(X[,j]),as.numeric(check1), method = 'spearman'))
  } 

  print(mean(vals, na.rm=T))
  res0[[i]] <- vals
  
}

save(res0, file='cibersort_spearman_correlation_reconstructed_expr.rda')
  
#library(CellMix)
#library(DeconRNASeq)
#res0 <- DeconRNASeq(datasets=as.data.frame(X), signatures=as.data.frame(B))
#deconOut <- as.data.frame(res0$out.all)

#qplot(res0$out.all - mayoDeconvMat)
#qplot(as.numeric(mayoDeconvMat[,'BOligo_OPALIN']), deconOut$Oligo_OPALIN)



