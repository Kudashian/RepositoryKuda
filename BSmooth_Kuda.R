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
head(granges(BSObj, seqnames(BSObj)))
ncol(BSObj) ##Columns -> samples
nrow(BSObj) ##Rows -> methylation loci
head(granges(BSObj, n = 5)) ##BSseq also has granges object which has general genomic regions.
##Command shows the location of first 5 identified methylation loci
head(getCoverage(BSObj, type = "M"), n = 5) ##BSseq has M matrix which has the number of reads covering methylation on a single loci.
##Columns -> num. of samples. Row -> methylation loci
head(getCoverage(BSObj, type = "Cov"), n = 5) ##BSseq has a Cov matrix which has the total number of reads covering a single loci.
####Columns -> num. of samples. Row -> methylation loci
getMeth()   ##Obtain methylation estimates for BSSeq objects. Smoothed and raw.

###########Apply BSmooth algorithm to data - estimates CpG site methylation by taking methylation info from neighbouring CpG sites#####################

BSObj_smoothed <- BSmooth(BSseq = BSObj,
                  ns = 1, ##Minimum number of loci in a smoothing window. 1 because I want to include any region a CpG.
                    h = 14, ##Minimum number of bases in a smoothing window. 14 because the shortest gene is 14 bases (testis determining factor).
                      maxGap = 10^8, ##Maximum gap between two methylation loci
                        keep.se = FALSE, ##Estimates standard errors to ne kept
                          ## chunkdim = NULL,
                            ## level = NULL, For storage of data either in-memory (null) or on-disk (HDF5Array file type)
                              verbose = FALSE) #Progress reports to be kept

###################Perform t-statistics on smoothed data####################################################

BSObj_Smoothed_tstat <- BSmooth.tstat(BSObj_Smoothed, ##The BSseq object that you are working on
                                        group1 = c(), ##The samples you want to compare - Treatment
                                          group2 = c(), ## Against the other sample - Control
                                            estimate.var = "group2", ## The sample which the statistic will be based on - Control
                                              local.correct = TRUE,
                                                verbose = TRUE)

Plot(BSObj_Smoothed_tstat)
plotManyRegions(BSObj_Smoothed_tstat)
##################Determine the differentially methylated regions in each sample######################
BSObj_DMR <- dmrFinder(BSObj_Smoothed&Controlled,
                          cutoff = ,          ##The cutoff of the t-statistics. This should be a vector of length two giving the (low, high)
                            qcutoff = ,       ##In case cutoff is NULL, compute the cutoff using these quantiles of the t-statistic
                              maxGap = ,      ##If two methylation loci are separated by this distance, break a possible DMR. This guarantees that the return DMRs have CpGs that are this distance from each other.
                                stat = ,      ##Which statistic should be used?
                                  verbose = ) ##Should the function be verbose?
############## Visualize the data before analyzing after performing t-statistics#####################################

plotManyRegions(BSseq = BSObj_DMR,               ##Regions - data frame with GRanges columns. Default is Null: Entire BSSeq columns plotted
                            Main = "",                                ##Main - Plot title
                              col = ,                              ##col - Color of methylation estimates
                                Addticks = ,                            ##Addticks - Ticks showing location of methylation loci
                                  Annotrack = ,                          ##Annotrack - An object of GRanges. Each component is a track, each track plotted as solid bars. E.g. CpG islands, exons, promoters
                                    BSseqStat = BSObj_Smoothed_tstat                        ##BSseqStat - Object of class BSseqStat. Adds a panel to show t-statistics.
                                      stat = ,                      ##Stat - Follows up to BSseqStat
                                        stat.col = ,                    ##stat.col - color of the stat plot
                                          regionCol = )                  ##regionCol - Color used for highlighting the region

BinomFit_BSObj_smoothed <- binomialGoodnessOfFit(BSseq = BSObj_DMR)  ##Tests whether the number of reads supporting methylation are independent and identically distributed across samples.
PoisGoF_BSObj_smoothed <- poissonGoodnessOfFit(BSseq = BSObj_DMR)    ##nQuantiles - number of evenly-spaced quantiles stored in the return object.
                                                                                                            ##Method - how the parameter is estimated
                                                                                                            ##Type - are the chisq, p-values or both being plotted
plot(BinomFit_BSObj_smoothed)
plot(PoisGoF_BSObj_smoothed)

round(colMeans(getCoverage(BSObj_smoothed)), 1)##The average coverage of CpGs on the two chromosomes
length(BSObj_smoothed) ##Number of CpGs identified
sum(rowSums(getCoverage(BSObj_smoothed) >= 1) == 6) ##Number of CpGs which are covered by at least 1 read in all 6 samples
sum(rowSums(getCoverage(BSObj_smoothed)) == 0) ## Number of CpGs with 0 coverage

############Decide on quality control parameters for the data###################
BS.cov <- getCoverage(BSObj_Smoothed_tstat)
keepLoci.ex <- which(rowSums(BS.cov[, BSObj_Smoothed_tstat$Type == ""] >= 2) >= 2 &
                      rowSums(BS.cov[, BSObj_Smoothed_tstat$Type == ""] >= 2) >= 2) ## $ refers to a specific column in a specific data frame.In thia instance - summarizes data where rows are features of interest and columns are samples
length(keepLoci.ex)
## This is the CpGs where at least 2 samples have at least 2X coverage
BSObj_Smoothed&Controlled <- BSObj_Smoothed_tstat[keepLoci.ex]
