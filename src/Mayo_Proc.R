

# proc mayo

load("data/mayo.rda")

# got ID translation from BioMart
genetable <- read.table('data/mart_export.txt', sep='\t', stringsAsFactors=F, header=T)

mayot2 <- cbind(data.frame(Gene.stable.ID=rownames(mayot)), mayot)

mayot3 <- inner_join(genetable, mayot2)

save(mayot3, file='data/mayo_hgnc.rda')