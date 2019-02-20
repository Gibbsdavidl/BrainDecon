
# David L Gibbs
# ISB
# 2018-1019


# goal
# 1.  cell specific gene expression signatures from single cell data
#
# 2.  deconvolved bulk brain tissue for cell type quantities
#
# 3.  using 1. & 2., estimate bulk expression present in specific cell types.

# workflow description #

# Data source:
# Allen Brain Institute
# http://celltypes.brain-map.org/rnaseq

#Middle Temporal Gyrus (MTG): This RNA-Seq data set is created from intact nuclei derived from frozen human brain specimens, to survey cell type diversity in the human middle temporal gyrus (MTG). In total, 15,928 nuclei from 8 human tissue donors ranging in age from 24-66 years were analyzed. Analysis of these transcriptional profiles reveals approximately 75 transcriptionally distinct cell types, subdivided into 45 inhibitory neuron types, 24 excitatory neuron types, and 6 non-neuronal types.

#----------------------------------------------------------------------
# 0. Preprocessing  src/dataPrep/proc_mtg_data.R
#
# 0a. Summed the exon and intron mapped single cell RNA-seq data:  
# 0b. Subset genes to be only protein coding genes. 
# 0c. Worked on coarse grained cell clusters... cell type + main marker gene (17 cell types).
# 

#-------------------------------------------------------------------

# 1. Used lasso to select genes for each cell type.

#######  BigReg.R  ################# 

#    1a.  Each cell type (in turn) was labeled as 1, all other types labeled as 0
#    1b.  Lasso used to select predictive genes
#    1c.  Top 25% of genes taken by T statistic magnitude.
#    1d.  Usually about 3500 genes
#    1e.  Cross validation on Lasso model... usually selected about x genes
#    1f.  Those genes taken, expression matrix subset, and random forest used to check predictive power.
#    1g.  Get importance of predictive genes.  Most important and T-test values can be used for selection.

#------------------------------------------------------------------

# 2. Cell signatures and gene expression deconvolution.

####### BigRF.R ##########

# Uses the draft cell signatures to predict cell types.
# End up with good predictive signatures for:
# decent cell signatures
  #Astro_FGFR3
  #OPC_PDGFRA
  #Oligo_OPALIN
  #Inh_PVALB
  #Inh_VIP
  #Inh_LAMP5
  #Inh_SST
  #Exc_FEZF2
  #Exc_RORB
  #Exc_THEMIS
  #Micro_TYROBP  # fixed
  #Exc_LINC00507 # 11% error on 2400 cells  # fixed removed 252 cells #

#------------------------------------------------------------------

####### BigSigs.R ##########

# 3. Then we take the big sets of genes, make a matrix of cell signatures,
#    using the original gene expression for each cell type, summed, avg.

#    Also the gene signature, for each gene across cells, sums to 1.

#    The signature matrix is created in order of important genes, for 
#    each cell type, found by the random forest predictions.

#    And try to find the smallest matrix by minimizing the condition number
#    which is estimated by kappa()

#    Found kappa is minimum at 39 genes


#------------------------------------------------------------------

################  MayoCibersortBootStrap.R #########################

#  Now to do deconvolution.

# Noticed that correlations very high for small number of genes ... but also
# easier to correlate ~20 genes.

# Rather then just take 39 genes...
# set up a probability for taking each of the signature genes.
# Since they're ranked top to bottom interms of how predictive they are.

# Prob of taking gene == (place in rank higher is better) / sum of ranks.

# top rank = 312 / sum(312 to 1)

#top probs
#  0.006389776 0.006369296 0.006348816 0.006328336 0.006307856 0.006287376

#mid probs
#  0.004362251 0.004341771 0.00432129

#bottom probs
#  2.662407e-04 2.457606e-04 2.252806e-04

# So, top gene is 0.00639 / 2.6624e-4
# [1] 24.000  times more likely to be picked!


# Did 1000 iterations where a gene signature was sampled
# and cibersort used for deconvolution

# the 1000 data frame results were summed up.

# found that many of the 'zeros' we saw when using a static signature, were just small numbers.

# Can also get a SD out too.


#------------------------------------------------------------------

################  expressionDivision.R #########################

# Then we divided the expression from the full cell matrix...
# which is average expression over cells for a cell type.

#  Divided_Gene_Expression_i_patent_j = Gene_ij * Cell_Sig_jk * Cell_Quantity_jk  

# so, taking the gene expression value for gene i for patient j,
# we multiply in the normalized expression expected for this cell type (Cell_Sig_jk)
# and multiply in the quantity of cell for patient j.









