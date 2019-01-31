
# processing the cibersort bootstrap

library(ggplot2)

load("~/Work/Brain_Deconvo/BrainRepo/cibersort_bootstrap.rda")

resMean <- res0[[1]]

for (i in 2:1000) {resMean <- resMean + res0[[i]]}

resMean <- resMean/1000

resSD <- res0[[1]]
for (i in 1:264) {
  for (j in 1:12) {
    vals <- c()
    for (k in 1:1000) {
      vals <- c(vals, res0[[k]][i,j])
    }
    resSD[i,j] <- sd(vals, na.rm = T)
  }
}



qplot(data=as.data.frame(res1), x=BAstro_FGFR3, geom='density')

qplot(data=as.data.frame(res1), x=BOligo_OPALIN, geom='density')

