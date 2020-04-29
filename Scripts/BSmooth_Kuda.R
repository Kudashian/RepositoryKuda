##Script for BSmooth by Kudakwashe Nyamupangedengu - 04 July 2019
#####################################################################
#Set to the correct working directory
setwd("C:/Users/kuda1/Documents/Postgraduate/Masters/Experimental Work/Bioinformatics/Scripts/BSmooth")
##Install BSSeq package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("bsseq") ##BSSeq is the class of which BSmooth is an object thereof
############################Load BSseq and all its tools################################################
BiocManager::install(version='devel')

BiocManager::install("bsseqData")
BiocManager::install("rtracklayer")
BiocManager::install("Rsamtools")
library(bsseq)
library(bsseqData)
source("C:/Users/kuda1/Documents/Postgraduate/Masters/Experimental Work/Bioinformatics/Scripts/BSmooth/MakeBSseqData.R")
###################Load Data###############################################
R1_CTRL_24 = read.table()
R2_CTRL_24 = read.table()
R3_CTRL_24 = read.table()
R1_D_24 = read.table()
R2_D_24 = read.table()
R3_D_24 = read.table()
R1_CTRL_72 = read.table()
R2_CTRL_72 = read.table()
R3_CTRL_72 = read.table()
R1_D_72 = read.table()
R2_D_72 = read.table()
R3_D_72 = read.table()
R1_CTRL_HEK = read.table()
R2_CTRL_HEK = read.table()
R1_D_HEK = read.table()
R2_D_HEK = read.table()
ListSampleNames <- c("R1_CTRL-24","R2_CTRL_24","R3_CTRL_24","R1_D_24","R2_D_24","R3_D_24",
                     "R1_CTRL-72","R2_CTRL_72","R3_CTRL_72","R1_D_72","R2_D_72","R3_D_72",
                     "R1_CTRL_HEK","R2_CTRL_HEK","R1_D_HEK","R2_D_HEK")
BSObj <- makebsseqdata(list(), sampleNames(ListSampleNames))

###################First look at the data you have inputted #####################################################
## Granges methods also work on BSseq objects
head(granges(BS.cancer.ex, seqnames(BS.cancer.ex)))
dim()
ncol() ##Columns -> samples
nrow() ##Rows -> methylation loci
head(granges(BS.cancer.ex, n = 5)) ##BSseq also has granges object which has general genomic regions.
##Command shows the location of first 5 identified methylation loci
head(getCoverage(, type = "M"), n = 5) ##BSseq has M matrix which has the number of reads covering methylation on a single loci.
##Columns -> num. of samples. Row -> methylation loci
head(getCoverage(, type = "Cov"), n = 5) ##BSseq has a Cov matrix which has the total number of reads covering a single loci.
####Columns -> num. of samples. Row -> methylation loci

###########Apply BSmooth algorithm to data - estimates CpG site methylation by taking methylation info from neighbouring CpG sites#####################

BSObj_smoothed <- BSmooth(BSseq = BSObj,
                  ns = 1, ##Minimum number of loci in a smoothing window. 1 because I want to include any region a CpG.
                    h = 14, ##Minimum number of bases in a smoothing window. 14 because the shortest gene is 14 bases (testis determining factor).
                      maxGap = 10^8, ##Maximum gap between two methylation loci
                        keep.se = FALSE, ##Estimates standard errors to ne kept
                          ## chunkdim = NULL,
                            ## level = NULL, For storage of data either in-memory (null) or on-disk (HDF5Array file type)
                              verbose = FALSE) #Progress reports to be kept

############## Visualize the data before analyzing#####################################

plotManyRegions(BSseq = BSObj_smoothed)               ##Regions - data frame with GRanges columns. Default is Null: Entire BSSeq columns plotted
                                                      ##Main - Plot title
                                                      ##col - Color of methylation estimates
                                                      ##Addticks - Ticks showing location of methylation loci
                                                      ##Annotrack - An object of GRanges. Each component is a track, each track plotted as solid bars. E.g. CpG islands, exons, promoters
                                                      ##BSseqStat - Object of class BSseqStat. Adds a panel to show t-statistics.
                                                      ##Stat - Follows up to BSseqStat
                                                      ##stat.col - color of the stat plot
                                                      ##regionCol - Color used for highlighting the region

BinomFit_BSObj_smoothed <- binomialGoodnessOfFit(BSseq = BSObj_smoothed)  ##Tests whether the number of reads supporting methylation are independent and identically distributed across samples.
PoisGoF_BSObj_smoothed <- poissonGoodnessOfFit(BSseq = BSObj_smoothed)    ##nQuantiles - number of evenly-spaced quantiles stored in the return object.
                                                                          ##Method - how the parameter is estimated
                                                                          ##Type - are the chisq, p-values or both being plotted

plot(BinomFit_BSObj_smoothed)
plot(PoisGoF_BSObj_smoothed)

round(colMeans(getCoverage(BSObj_smoothed)), 1)##The average coverage of CpGs on the two chromosomes
length(BSObj_smoothed) ##Number of CpGs identified
sum(rowSums(getCoverage(BSObj_smoothed) >= 1) == 6) ##Number of CpGs which are covered by at least 1 read in all 6 samples
sum(rowSums(getCoverage(BSObj_smoothed)) == 0) ## Number of CpGs with 0 coverage

############Decide on quality control parameters for the data###################
BS.cov <- getCoverage(BSObj_smoothed)
keepLoci.ex <- which(rowSums(BS.cov[, BSObj_smoothed$Type == ""] >= 2) >= 2 &
                       rowSums(BS.cov[, BSObj_smoothed$Type == ""] >= 2) >= 2) ## $ refers to a specific column in a specific data frame.In thia instance - summarizes data where rows are features of interest and columns are samples
length(keepLoci.ex)
## This is the CpGs where at least 2 samples have at least 2X coverage
BSObj_Smoothed_fit <- BSObj_Smoothed_fit[keepLoci.ex]
###################Perform t-statistics on smoothed data####################################################
##Perform t-statistics on smoothed data

BSObj_Smoothed_tstat <- BSmooth.tstat(BSObj_Smoothed_fit, ##The BSseq object that you are working on
                                    group1 = c("C1", "C2", "C3"), ##The samples you want to compare, in this case - Cancer
                                    group2 = c("N1", "N2", "N3"), ## Against the other sample - Normal
                                    estimate.var = "group2", ## The sample which the statistic will be based on - Control
                                    local.correct = TRUE,
                                    verbose = TRUE)

plot(BSObj_Smoothed_tstat)
