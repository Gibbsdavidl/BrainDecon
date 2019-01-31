# getting the cell-type specific data from the Barres lab

library(org.Mm.eg.db)
library(readxl)
library(tidyr)
library(dplyr)
library(biomaRt)
library(RUnit)

counts <- read_excel("data/barreslab_rnaseq.xlsx")
keys <- counts$`Gene symbol`
lookup <- AnnotationDbi::select(org.Mm.eg.db, keys=keys, keytype="ALIAS", columns=c("ENTREZID", "ENSEMBL"))
ensmbl.counts <- left_join(counts, lookup, by = c("Gene symbol" = "ALIAS"))

#-----------------------------------------------------------------------------------
# make a map of all human and mouse homologs from biomart
ensembl = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mouse_mapping <- getBM(attributes = c("ensembl_gene_id", "entrezgene"), mart = ensembl)
mouse_mapping <- rename(mouse_mapping, "mouse_gene_id" = "ensembl_gene_id")
mouse_mapping <- rename(mouse_mapping, "mouse_entrez" = "entrezgene")

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
human_mapping <- getBM(attributes = c("ensembl_gene_id","entrezgene"), mart = ensembl)
human_mapping <- rename(human_mapping, "human_gene_id" = "ensembl_gene_id")
human_mapping <- rename(human_mapping, "human_entrez" = "entrezgene")

human_mouse <- getBM(attributes = c("ensembl_gene_id","mmusculus_homolog_ensembl_gene"), mart = ensembl)
# get rid of ids with no homolog
#human_mouse <- human_mouse[which(human_mouse$mmusculus_homolog_ensembl_gene != ""),] # optional

human_mouse <- full_join(human_mouse, human_mapping, by = c("ensembl_gene_id" = "human_gene_id"))
human_mouse <- full_join(human_mouse, mouse_mapping, by = c("mmusculus_homolog_ensembl_gene" = "mouse_gene_id"))
#------------------------------------------------------------------------------------

human_mouse$mouse_entrez_str <- as.character(human_mouse$mouse_entrez)
human_mouse_naomit <- human_mouse[!is.na(human_mouse$mouse_entrez_str),]

ensmbl.counts2 <- inner_join(ensmbl.counts, human_mouse_naomit, by = c("ENTREZID" = "mouse_entrez_str"))

ensmbl.counts3 <- ensmbl.counts2 %>% 
  dplyr::select("Gene symbol", "ENSEMBL", "mouse_entrez", "human_entrez", "ensembl_gene_id", "Astrocytes", "Neuron", 
                "Oligodendrocyte Precursor Cell", "Newly Formed Oligodendrocyte", "Myelinating Oligodendrocytes", 
                "Microglia", "Endothelial Cells") %>% 
  rename("human_ensembl" = "ensembl_gene_id") %>%
  rename("mouse_ensemble" = "ENSEMBL")

barres.df <- ensmbl.counts3

barres.percent <- round(barres.df[5:11]/rowSums(barres.df[5:11]), 3)
barres.percent <- cbind(barres.df[1:4], barres.percent)

#------------------------------------------------------------------------------------
# function to give cell-type percentages using a gene as input


getCelltype <- function(gene.name, barres.df, cutoff){
  results <- NA

  if(gene.name %in% barres.df$mgi_symbol)
  {
    results <- subset(barres.df, mgi_symbol == gene.name)[1,-c(1:4)]
  }
  else if(gene.name %in% barres.df$mouse_ensembl)
  {
    results <- subset(barres.df, mouse_ensembl == gene.name)[1, -c(1:4)]
  }
  else if(gene.name %in% barres.df$hgnc_symbol)
  {
    results <- subset(barres.df, hgnc_symbol == gene.name)[1,-c(1:4)]
  }
  else if(gene.name %in% barres.df$human_ensembl)
  {
    results <- subset(barres.df, human_ensembl == gene.name)[1, -c(1:4)]
  }
  else{
    browser()
    xyz = 99
  }
  results <- as.data.frame(t(results))
  colnames(results)[1] <- "percent"
  results$celltype <- rownames(results)

  #results <- results %>% arrange(desc(percent)) %>% filter(percent >= cutoff)
  return(results)
}

topCelltype <- function(gene.name, df, cutoff){
  names <- names(which(as.list(tbl.sub[gene.name,]) > cutoff))
  names.collapsed <- paste(names, collapse=",")
  return(names.collapsed)
  }
test_topCelltype <- function(){
  checkEquals(topCelltype("DYNC1I1", w.targets.celltype, 0.1), "Neuron")
  checkEquals(topCelltype("DYNC1I1", w.targets.celltype, 0.99), "")
  checkEquals(topCelltype("DYNC1I1", w.targets.celltype, 0.0),
              "Astrocytes,Neuron,Oligodendrocyte Precursor Cell,Newly Formed Oligodendrocyte,Myelinating Oligodendrocytes,Microglia,Endothelial Cells")

  checkEquals(x, "Neuron")
}
runTests <- function()
{

  test_topCelltype()
}

w.targets.celltype$celltype <- unlist(lapply(rownames(w.targets.celltype), function(gene.name) topCelltype(gene.name, w.targets.celltype, 0.2)))
df <- w.targets.celltype[,"celltype", drop = FALSE]
#------------------------------------------------------------------------------------

w.targets <- scan("~/Alzheimers/sage/wall_of_targets", what="", sep="\n")
# swap out one name for another
w.targets <- gsub("FAM63A", "MINDY1", w.targets)
w.targets <- gsub("U1-70k", "SNRNP70", w.targets)
w.targets <- gsub("U1-C", "SNRPC", w.targets)
w.targets <- gsub("SmN", "SMN1", w.targets)
w.targets <- gsub("SmB", "SNRPB", w.targets)
w.targets <- gsub("NGFRAP1", "BEX3", w.targets)

all.wall <- lapply(w.targets, getCelltype, barres.df = barres.percent, cutoff = 0.2)
names(all.wall) <- w.targets

# without a cutoff
all.wall.full <- lapply(w.targets, getCelltype, barres.df = barres.percent, cutoff = 0.0)
names(all.wall.full) <- w.targets



#---------------------------------------------------------------------------------------
# Subset the data for the Wall of Targets

temp <- barres.percent[which(barres.percent$hgnc_symbol %in% w.targets),]
temp$mgi_symbol <- NULL
temp$mouse_ensembl <- NULL
temp$human_ensembl <- NULL
dups <- which(duplicated(temp$hgnc_symbol))
temp <- temp[-dups,]
