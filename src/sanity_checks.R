
# let's check the new cell matrices #

library(ggplot2)
library(readr)

load('../mayoCellMats.rda')

load('results/BigSig_Results.rda')

protGenes <- read_tsv("data/protein-coding_gene.txt")


impSymbols <- lapply(impGenes, function(a) {
  sapply(names(a), function(b) as.character(protGenes$symbol[protGenes$entrez_id == b])[1])  
})

# gene of interest
goi <- "MOG"  # found in Oligo_OPALIN, entrez 4340

#[1] "Astro_FGFR3"   "Exc_FEZF2"     "Exc_LINC00507" "Exc_RORB"      "Exc_THEMIS"    "Inh_LAMP5"     "Inh_PVALB"     "Inh_SST"      
#[9] "Inh_VIP"       "Micro_TYROBP"  "Oligo_OPALIN"  "OPC_PDGFRA"  

# let's build a matrix to plot 
a <- as.numeric(cellList[[10]][GeneSymbols == goi, -1])
b <- as.numeric(cellList[[11]][GeneSymbols == goi, -1])

df <- data.frame(CellType=c(rep.int(1,264), rep.int(2,264)), ExprVals=c(a,b))

df$CellNames <- ifelse(df$CellType == 1, yes = "Micro_TYROBP", no = "Oligo_OPALIN")

qplot(data=df, x=log(ExprVals), col=CellNames, geom='density',
      main = 'Cell Type Deconvolved Expression, gene: MOG, impt for Oligos', )


