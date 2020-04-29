##Script for eDMR by Kudakwashe Nyamupangedengu - 04 July 2019
#########################Setting the correct working directory#######################################################################################

setwd("C:/Users/kuda1/Documents/Postgraduate/Masters/Experimental WOrk/Bioinformatics/Scripts/eDMR")
##########################Installation and Loading###################################################################################################

install.packages( c("data.table", "mixtools", "devtools"))
BiocManager::install("GenomicRanges")
BiocManager::install("IRanges")
install.packages("remotes")
remotes::install_github("ShengLi/edmr")
## Load all packages

library(edmr)
library(GenomicRanges)
library(IRanges)
library(mixtools)
library(data.table)
#################################Load necessary files for analysis##########################################################################
                                                                                       ##BEDfiles for annotation
cpgi.anno("Example_CpGIslands.bed", shore.width = 2000, shelf.width = 2000)
genebody.anno("Example_Gene.bed")

myMethylatedRegions = source("calculatedMethylDiffobject.txt")                         ##Load MethylKit object as input file
##########################################Get Mixtools models of CpGs##############################################################################################
Mixmodel_myMethylatedRegions = myDiff.to.mixmdl(myMethylatedRegions, plot = FALSE, main = "Sample Bimordal Normal Distribution")   ##Mixtools model of the object

##########################################Weighted Cost Function and Bimordal normal distribution plots########################################################################
plotMdl1(Mixmodel_myMethylatedRegions, subtitle = "Bimordal normal distribution")
plotCost(Mixmodel_myMethylatedRegions, main = "Weighted Cost Function")
##########################################Filter the DMRs#######################################################################
Filtered_myMethylatedRegions = filter.dmr(Mixmodel_myMethylatedRegions,                                                    ##Filters significant DMRs
                                            DMR.qvalue = 0.01,                                                    ##Qvalue cutoff for DMC definition
                                                mean.meth.diff = 20,                                              ##Cutoff of the DMR mean methylation difference
                                                    num.CpGs = 3,                                                 ##Cutoff no. of CpGs in each region to call DMR
                                                      num.DMCs = 3)                                               ##Cutoff no. of DMCs
#################################Create eDMR object of class GRanges#################################################################################

myDiffMethylatedRegions = edmr(myDiff, step = 100,                                       ##Numerical value for calculating auto-correlation
          dist = "none",                                                                 ##Distance cut-off to call a gap for DMR. Default none=automatically determined by bimodal normal distribution
              DMC.qvalue = 0.01,                                                         ##Qvalue cutoff for DMC definition
                  DMC.methdiff = 25,                                                     ##Methylation difference cutoff for DMC definition
                      num.DMCs = 1,                                                      ##Cutoff of the number of DMCs in each region to call DMR
                          num.CpGs = 3,                                                  ##Cutoff number of CpGs
                              DMR.methdiff = 20,                                         ##Cutoff of the DMR mean methylation difference
                                  plot = TRUE,
                                      main = "Binomial Distribution of Example data",
                                          mode = 1,                                      ##Mode of call DMRS: 1. Using all CpGs together. 2. Use unidirectional CpGs to call DMRs
                                              ACF = TRUE,
                                                  fuzzypval = 0.04)

##########################################Annotate gene bodies using a bed file###############################################################################################

genebody.anno()
############################################Get genes in DMRs##########################################################################################################################

get.dmr.genes(myDMR = myDMR,                                                             ##Get gene list based on gene body granges
                    subject,                                                             ##GRanges used to annotate DMRs.
                        id.type = "Gene.symbol")                                         ##The col names that will be used to annotate the DMR
############################################Get hyper/hypo methylated regions##########################################################################################################3

myHyperMeth = get.hyper.dmr(myFilteredDMR)
myHypoMeth = get.hypo.dmr(myFilteredDMR)
