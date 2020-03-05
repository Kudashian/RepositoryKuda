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
###################Load Data###############################################
data("BS.cancer.ex")
data("BS.chr22")
###################First look at the data you have inputted #####################################################
## Granges methods also work on BSseq objects
head(granges(BS.cancer.ex, seqnames(BS.cancer.ex)))
dim(BS.cancer.ex)
ncol(BS.cancer.ex) ##Columns -> samples
nrow(BS.cancer.ex) ##Rows -> methylation loci
head(granges(BS.cancer.ex, n = 5)) ##BSseq also has granges object which has general genomic regions.
##Command shows the location of first 5 identified methylation loci
head(getCoverage(BS.cancer.ex, type = "M"), n = 5) ##BSseq has M matrix which has the number of reads covering methylation on a single loci.
##Columns -> num. of samples. Row -> methylation loci
head(getCoverage(BS.cancer.ex), n = 5) ##BSseq has a Cov matrix which has the total number of reads covering a single loci
####Columns -> num. of samples. Row -> methylation loci

###########Apply BSmooth algorithm to data - estimates CpG site methylation by taking methylation info from neighbouring CpG sites#####################

BS.chr22_smoothed <- BSmooth(BSseq = BS.chr22,
                  ns = 70, ##Minimum number of loci in a smoothing window
                    h = 1000, ##Minimum number of bases in a smoothing window
                      maxGap = 10^8, ##Maximum gap between two methylation loci
                        keep.se = FALSE, ##Estimates standard errors to ne kept
                          ## chunkdim = NULL,
                            ## level = NULL, For storage of data either in-memory (null) or on-disk (HDF5Array file type)
                              verbose = FALSE) #Progress reports to be kept
BS.cancer.ex.Smoothed <- BSmooth(BSseq = BS.cancer.ex,
                                 ns = 70,
                                  h = 1000,
                                    maxGap = 10^8,
                                      keep.se = FALSE,
                                        BPPARAM = MulticoreParam(workers = 1), 
                                          verbose = TRUE)
############## Visualize the data before analyzing#####################################

plotManyRegions(BSseq = BS.cancer.ex.fit)
##plotRegion(BSseq = BS.chr22_smoothed)
BinomFit_BS.Chr22 <- binomialGoodnessOfFit(BSseq = BS.chr22_smoothed) ##Tests whether the number of reads supporting
                                                                          ##methylation are independent and identically distributed across samples
plot(BinomFit_BS.Chr22)
PoisGoF_BS.Chr22 <- poissonGoodnessOfFit(BSseq = BS.chr22_smoothed)
plot(PoisGoF_BS.Chr22)

round(colMeans(getCoverage(BS.cancer.ex)), 1)##The average coverage of CpGs on the two chromosomes
length(BS.cancer.ex) ##Number of CpGs identified
sum(rowSums(getCoverage(BS.cancer.ex) >= 1) == 6) ##Number of CpGs which are covered by at least 1 read in all 6 samples
sum(rowSums(getCoverage(BS.cancer.ex)) == 0) ## Number of CpGs with 0 coverage

############Decide on quality control parameters for the data###################
BS.cov <- getCoverage(BS.cancer.ex.fit)
keepLoci.ex <- which(rowSums(BS.cov[, BS.cancer.ex$Type == "cancer"] >= 2) >= 2 &
                       rowSums(BS.cov[, BS.cancer.ex$Type == "normal"] >= 2) >= 2) ## $ Summarizes data where rows are features of interest and columns are samples
length(keepLoci.ex)
## This is the CpGs where at least 2 cancer and 2 normal samples have at least 2X coverage
BS.cancer.ex.fit <- BS.cancer.ex.fit[keepLoci.ex]
###################Perform t-statistics on smoothed data####################################################
##Perform t-statistics on smoothed data

BS.cancer.ex.tstat <- BSmooth.tstat(BS.cancer.ex.fit, ##The BSseq object that you are working on
                                    group1 = c("C1", "C2", "C3"), ##The samples you want to compare, in this case - Cancer
                                    group2 = c("N1", "N2", "N3"), ## Against the other sample - Normal
                                    estimate.var = "group2", ## The sample which the statistic will be based on - Control 
                                    local.correct = TRUE,
                                    verbose = TRUE)

plot(BS.cancer.ex.tstat)
