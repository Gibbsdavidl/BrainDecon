
## split up the cell types ##

load('../patientList_expr_portions.rda')
library(data.table)

genesColumn <- data.table(data.frame(GeneSymbols=rownames(patientList[[1]])))

patientDataTables <- lapply(patientList, function(a) data.table(a))
rm(patientList); gc()

length(patientList)
# 264, one for each patient

dim(patientList[[1]])
#[1] 14340    12 ... we have 14K genes for 12 cells.

# a list of matrices by cell type
cellList <- list()

# for each cell type
for (ci in colnames(patientDataTables[[1]])) {
  x <- genesColumn
  #    for each patient
  for (pi in 1:264) {
    #       add that cell type column to the big matrix
    x <- cbind(x, patientDataTables[[pi]][,..ci])  
  }  
  cellList[[ci]] <- x  
}


save(cellList, file='mayoCellMats.rda')




