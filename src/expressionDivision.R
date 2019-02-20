

# Working with the summarized deconvolution #

load('results/cibersort_bootstrap_summary.rda')

load('data/full_cell_expression.rda')

load('data/mayo_hgnc.rda')

cellExpr <- tdatSum
rownames(cellExpr) <- symbols
dim(cellExpr)
#[1] 14340    12



colnames(resMean) <- c("Astro_FGFR3", "Exc_FEZF2","Exc_LINC00507","Exc_RORB","Exc_THEMIS","Inh_LAMP5", "Inh_PVALB", "Inh_SST", "Inh_VIP", "Micro_TYROBP", "Oligo_OPALIN", "OPC_PDGFRA")

# have obs. express for patient pi
# have means and sds for each sample and data type. #
# and have expression for cell types.

# want fraction of expression for each gene in each sample for each cell

patientList <- list()

# For patient pi
for (pi in colnames(mayot3)[-c(1,2)]) {
  # For gene gk in expression for pi
  
  patientMat <- matrix(0, nrow=nrow(cellExpr), ncol=ncol(cellExpr))
  rownames(patientMat) <- rownames(cellExpr)
  colnames(patientMat) <- colnames(cellExpr)
  
  for (gk in rownames(cellExpr)) {
    if (gk %in% mayot3$HGNC.symbol) {
      exprVal <- exp(mayot3[mayot3$HGNC.symbol == gk,pi])
      cellVals <- cellExpr[gk,] / sum(cellExpr[gk,])
      deconVals <- resMean[pi,]
      
      portioningVec <- cellVals * deconVals
      portioningVec <- portioningVec / sum(portioningVec)
      
      patientMat[gk,] <- exprVal * portioningVec * 1000
    }
  }

    patientList[[pi]] <- patientMat
}


save(patientList, file='patientList_expr_portions.rda')




