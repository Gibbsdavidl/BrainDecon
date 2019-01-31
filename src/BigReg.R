
#library(ggplot2)
#library(readr)
#library(stringr)
#require(glmnet)

# use regularized regression to select genes #

#load("~/Work/Brain_Deconvo/rdata/mtg_exon_intron.rda")
# 15928 cells x 50281 genes

#geneNames <- read_csv('data/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_genes-rows.csv')
# 50281 for measure mac rows and number genes in geneNames

#load("~/Work/Brain_Deconvo/rdata/samp_info.rda")
#idname <- "Exc L2-3 LINC00507 FREM3"
#idcode <- ifelse(mtgsamp$cluster == idname, 1, 0)
#protGenes <- read_tsv("~/Data/protein-coding_gene.txt")
#dat <- dat[rownames(dat) %in% protGenes$entrez_id,]
#dat <- t(dat)
#geneVar <- apply(dat, 2, var, na.rm=T)
#library(ggplot2)
#qplot(geneVar)
#summary(geneVar)
#sum(geneVar < 0.7178)
#qplot(geneVar[geneVar > 0.7178])
#dat2 <- dat[,geneVar > 0.7178]
# saved the above #

library(ggplot2)
library(readr)
library(stringr)
require(glmnet)
library(doParallel)

registerDoParallel(4)

load("~/Work/Brain_Deconvo/big_rec_dat.rda")

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

resList <- list()

# for a single cell type and for a number of iterations 
dim(dat)
#[1] 15928 14340

for (ci in cset13) {
  
  idcode <- ifelse(cset13 == ci, 1, 0)
  
  # first do t-test for each variable ... remove genes that have T < 10
  tlist <- apply(dat, 2, function(a) abs(t.test(a[idcode == 0], a[idcode == 1])$statistic))
  
  # get the top 25%
  n <- 0.25 * length(tlist)
  idx <- order(tlist, decreasing = T)[1:n]
  
  #    fit lasso
  fitdat <- dat[,idx]
  fitdatvar <- apply(fitdat, 2, var, na.rm=T)
  fitdat <- fitdat[,fitdatvar > 2]
  
  # needed to remove some low variance genes in case of micro_tyrobp
  cvfit <- cv.glmnet(y = idcode, x=fitdat, family = "binomial", type.measure = "class", parallel = T)
  
  # put top variables into a list
  cvcof <- coef(cvfit, s = "lambda.min")
  nonZeds <- cvcof[cvcof[,1] != 0,]
  
  theseRes <- list()
  theseRes[['n']] <- n
  theseRes[['names']] <- names(nonZeds)[-1]
  theseRes[['lambda.min']] <- cvfit$lambda.min
  theseRes[['cverror']] <- cvfit$cvm[which(cvfit$lambda == cvfit$lambda.min)]
  theseRes[['numNonZed']] <- cvfit$nzero[which(cvfit$lambda == cvfit$lambda.min)]
  theseRes[['cvfit']] <- cvfit
  
  resList[[ci]] <- theseRes
}

save(resList, file='brain_lasso_results_v2.rda')

  