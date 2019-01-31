
# PSEUDOSORT
# A reimplementation of CIBERSORT
# David L Gibbs, dgibbs@systemsbiology.org
# June 9, 2017

#A variety of GEP deconvolution methods have been proposed,
#most of which model an mRNA mixture m by a system of linear equations,
#corresponding to a weighted sum of cell type-specific GEPs3–9,12,13,17,19.

#Let B denote a GEP signature matrix and let f denote a vector consisting
#of the unknown fractions of each cell type in the mixture. Then the problem
#of GEP deconvolution can be represented by m = f x B, provided that B contains
#more marker genes than cell types (i.e., the system is overdetermined4).
#The preponderance of genes in whole transcriptome studies renders this
#requirement trivial in practice. If the linearity argument is biologically
#plausible, as previous studies imply4,12,13,20, then genes with expression
#profiles enriched in each cell type can be leveraged to impute unknown cell
#fractions from mixture profiles5.

#CIBERSORT, ν-SVR is applied with a linear kernel to solve for f,
#and the best result from three values of ν = {0.25, 0.5, 0.75}
#is saved, where ‘best’ is defined as the lowest root mean squared error
#between m and the deconvolution result, f x B.

#Our current implementation of CIBERSORT executes ν-SVR using the
#‘svm’ function in the R package, ‘e1071’.
#Regression coefficients are extracted with the following R command:

#coef <- t(model$coefs) %*% model$SV

#Negative SVR regression coefficients are subsequently set to zero
#(as done for LLSR), and the remaining regression coefficients are
#normalized to sum to 1, yielding a final vector of estimated cell type
#fractions, f (notably, f denotes relative, not absolute fractions of each
#cell type from B in m). To decrease running time and promote better
#overall performance, both B and m are each normalized to zero mean and
#unit variance prior to running CIBERSORT. As previously suggested for
#other linear deconvolution methods, CIBERSORT works best on expression
#values in non-log linear space[20].

# m is the mixture matrix.
# B is the signature matrix
# will return f, the unknown fractions.

library(e1071)

#' Use weightNorm to normalize the SVM weights
#'
#' @param w  The weight vector from fitting an SVM
#' @export
#' @usage w1 <- weightNorm(t(fit1$coefs) %*% fit1$SV), where fit comes from  <- svm(m~B, nu=0.25, kernel="linear")
weightNorm <- function(w) {
  w[w<0] <- 0
  return(w/sum(w))
}


#' Use PSEUDOSORT to estimate the cell count percentage
#'   To Do: confirm correct model & p-value estimation via monty carlo methods & make parallel & split funcs.
#'
#' @param m  a matrix represenging the mixture (genes X 1 sample)
#' @param B  a matrix representing the references (genes X cells)
#' @export
#' @usage w2 <- PSEUDOSORT(y1, X)
PSEUDOSORT <- function(m,B) {

  # Here, the m should be subset to match B #

  # split into functions:
  # (1) main -- calls (2) and (3)
  # (2) cibersort model selection function, returns weights
  # (3) p-value function, calls (2) with random m's

  # three models are fit with different values of nu
  fit1 <- svm(m~B, nu=0.25, kernel="linear", scale=T, type="nu-regression")
  fit2 <- svm(m~B, nu=0.50, kernel="linear", scale=T, type="nu-regression")
  fit3 <- svm(m~B, nu=0.75, kernel="linear", scale=T, type="nu-regression")
  # these w's are the cell fractions
  w1 <- weightNorm(t(fit1$coefs) %*% fit1$SV)
  w2 <- weightNorm(t(fit2$coefs) %*% fit2$SV)
  w3 <- weightNorm(t(fit3$coefs) %*% fit3$SV)
  # return the model with the smallest mean sq error
  err1 <- sqrt( sum( (m - B %*% t(w1))^2 )/nrow(m) )
  err2 <- sqrt( sum( (m - B %*% t(w2))^2 )/nrow(m) )
  err3 <- sqrt( sum( (m - B %*% t(w3))^2 )/nrow(m) )
  resIdx <- which(c(err1,err2,err3) == min(c(err1,err2,err3)))[1]
  return(list(w1,w2,w3)[[resIdx]])
}
