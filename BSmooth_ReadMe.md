Step 1: Make a BSseq object
	makeBSseqData(dat, sampleNames) function. Function has to be sourced.
		dat is a list of multiple data frames from biological replicates. One element = data from 1 rep.
		data frame must be = 1 chromosome number ("chr")
				     2 genomic co-ordinates ("pos")
				     3 read coverage ("N")
				     4 No. of reads showing methylation ("X")
		sampleNames is a vector of characters for the sample names. Length of vector should match the length of the input list.

Step 2: Perform smoothing
	BSmooth() function
		ns - min. number of loci in smoothing window. 1 means at least one loci (CpG per window)
		h - min. number of bases in window. Shortest gene (Tests Determining Factor) is 14 bases.

Step 3: Plot regions of DMRs
	PlotRegions()
	regions - data frame GRanges (start, end & chr)
	main - Plot title
	col - Colour of methylation estimates
	addticks - Ticks showing location of methylation loci

Step 4: Distribution plots
	Statistics for BSseq objects including plotting capabilities
	BinomialGoodnessOfFit()
	PoissonGoodnessOfFit()         ** $ refers to a specific column is a specific data frame
