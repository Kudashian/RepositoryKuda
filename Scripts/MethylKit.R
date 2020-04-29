## Script for MethylKit by Kudakwashe Nyamupangedengu - 11/06/2019
#===========================================================================
## Make sure you are working in the correct directory

setwd("/home/studentsgh129/Kudakwashe/DNA_Methylation_Analysis/Scripts/MethylKit")
setwd("C:/Users/kuda1/Documents/Postgraduate/Masters/Experimental Work/Bioinformatics/Scripts/methylKit")
## Install Methylkit package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("methylKit")

## Load package

library("methylKit")
#=====================Reading in and cleaning up input files====================================================
## Read input files, IDs them as test/control and declares them as treatment or control

InputFileList <- list("R1_CTRL_24hr.txt", "R2_CTRL_24hr.txt", "R3_CTRL_24hr.txt", ##Make a R list of methylation call files to be analysed.
                      "R1_D_24hr.txt", "R2_D_24hr.txt", "R3_D_24hr.txt",          ##The files have to be in the directory. So I need to input the paths to each methylation call file.
                      "R1_CTRL_96hr.txt", "R2_CTRL_96hr.txt", "R3_CTRL.txt",
                      "R1_MONO_96hr.txt", "R1_MONO_96hr.txt", "R1_MONO_96hr.txt",
                      "R1_C_HEK.txt", "R2_C_HEK.txt", "R1_T_HEK.txt", "R2_T_HEK.txt")


MyMethObj = methRead("R1_CTRL_24hr.txt" ,  ## Making a object containing the input methylation call files.
               sample.id = List("R1C24", "R2C24", "R3C24" ##Assign an identifier for each file
                              , "R1T24", "R2T24", "R3T24"
                              , "R1C96", "R2C96", "R3C96"
                              , "R1T96", "R2T96", "R3T96"
                              , "R1CHEK", "R2CHEK"
                              , "R1THEK", "R2THEK") ,
                 assembly = "Hg19" , ## Defines the genome assembly. Should be uniform for each call file.
                     ##Pipeline = list(fraction=FALSE, chr.col = , start.col = , end.col, coverage.col = , strand.col = , freqC.col = ) ## If the input file is a generic text file, this is NB!
                        treatment = c(0,0,0,1,1,1,0,0,0,1,1,1,0,0,1,1) ,      ## Denotes which samples are control and which are test.
                           context = "CpG", ##Denotes methylation context
                             dbtype = "tabix",
                              dbdir =     )  ##Directory where tabix file will be saved


#===================First look at raw data statistics===================================================
## First look at the methylation data
## % methylation statistics per base
GetMethylationStats(MyMethObj[], plot = TRUE, both.strands = FALSE)   ##% Methylation plot should have two peaks on both ends: C's should be methylated or not methylated
                                                                      ##Both.strands (If true) looks at each strand separately
                                                                      ##Square brackets to look at a sample in the object individually
## % Coverage statistics per base
GetCoverageStats(MyMethObj[], plot = TRUE, both.strands = FALSE)      ##PCR Bias: A second peak on the right side of the peak
                                                                      ##Both.strands (If true) looks at each strand separately
## Filter the methylation data: Set coverage minimum. Remove PCR biased bases.
FilteredMyMethObj = filterbycoverage(MyMethObj, lo.count = , lo.percent = NULL,
                                     hi.count = NULL, hi.percent = 99.9)      ##lo.count/hi.count - Coverage minimum/maximums
                                                                              ##lo.percent/hi.percent - %methylation of that base (Used to filter out biased bases)
#====================Merging samples for grouped comparisons===================================================
## We merge all samples to one object, Destrand merges reads on both strands for better coverage

MergedMeth = unite(FilteredMyMethObj, destrand = TRUE)
#===================Correlation and PCA of merged data======================================================
## Correlation - Shows the relationship between the samples. Returns correlation coefficients or scatterplots.

GetCorrelation(MergedMeth, plot = TRUE)

##We can cluster samples based on the similarity of the methylation profiles.

ClusterSamples(MergedMeth, dist = "correlation", method = "ward.D", plot = TRUE)
PCASamples(MergedMeth)

#====================Differentially Methylated Regions========================
DiffMergedMeth = calculateDiffMeth(MergedMeth)                    ##covariates - a data frame containing covariates
                                                                  ##adjust - corrects p-values, default is SLIM
                                                                  ##test - statistical test for determining methylation differences
                                                                  ##Overdispersion = "MN". Statistical state automatically changes to F-test

Hyper.DiffMergedMeth = get.methylDiff(DiffMergedMeth,
                                        difference = ,            ##Cutoff for absolute value of methylation percentage change between test & control
                                          qvalue = ,              ##Cutoff for qvalue of differential methylation satistic
                                              type = "")          ##Species type of differentially methylated regions to be returned: hyper/hypo/all

diffMethPerChr(DiffMergedMeth,                                    ##Gets number of hyper/hypo methylated regions
                plot = TRUE ,                                     ##Horizontal barplots for proportion of diff methylated bases
                  qvalue.cutoff = ,                               ##Cutoff for q-value
               meth.cutoff = )                                    ##Cutoff for percent methylation difference
#======================Annotating Differentially methylated regions==============
BiocManager::install('genomation')
library(genomation)                                       ##Package used for annotation

Genome_Obj = readGeneric("",                              ##Bed file containing gene information
                          header = TRUE,                  ##If the file has column headings, requires most basic: chr, start,end
                              keep.all.metadata = TRUE,   ##Determines if the extra columns (outside of the core) should be kept or not
                                 zero.based = TRUE)       ##Tells if the ranges started from 0 or 1

annotateWithGeneParts(as(DiffMergedMeth, "GRanges"), Genome_Obj) ##Annotate the DMRs with promoter/intron/exon data
