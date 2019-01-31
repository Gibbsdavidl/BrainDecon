

# data prep
library(readr)

# use regularized regression to select genes #

#load('rdata/samp_info.rda')
mtgsamp <- read_csv('data/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_samples-columns.csv')
dat1 <- read_csv('data/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_exon-matrix.csv.gz')
dat2 <- read_csv('data/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_intron-matrix.csv.gz')
all(dat1$X1 == dat2$X1)

dat1 <- dat1[,-1]
dat2 <- dat2[,-1]
dat <- dat1+dat2

rm(dat1, dat2); gc();

save(dat, "~/Work/Brain_Deconvo/rdata/mtg_exon_intron.rda")

load("~/Work/Brain_Deconvo/rdata/mtg_exon_intron.rda")
# 15928 cells x 50281 genes

geneNames <- read_csv('data/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_genes-rows.csv')
# 50281 for measure mac rows and number genes in geneNames

mtgsamp <- read_csv('data/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_samples-columns.csv')

protGenes <- read_tsv("data/protein-coding_gene.txt")
dat <- dat[rownames(dat) %in% protGenes$entrez_id,]
dat <- t(dat)
geneVar <- apply(dat, 2, var, na.rm=T)
library(ggplot2)
#qplot(geneVar)
summary(geneVar)
sum(geneVar < 0.7178)
qplot(geneVar[geneVar > 0.7178])
dat2 <- dat[,geneVar > 0.7178]

dat <- dat2
save(dat, 'big_rec_dat.rda')
