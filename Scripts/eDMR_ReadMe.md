Input file: MethylKit object (created by CalculateDiffMeth function). chr start end strand pvalue qvalue methdiff
	    The object has to include all the CpG sites for the signficant test, not just the DMCs.
	    I need to export a MethylKit object for each sample

Step 1: Load input objects from directory
	source() function

Step 2: Fit binomial normal distribution to CpGs. Mixtools is used to fit models to the data.
	myDiff.to.mixmdl() function

Step 3: Plot the distribution and weighted cost function of the models
	plotCost() function
	plotMdl1() function

Step 4: Filter the methylated regions for statistically different regions. Stats are the models from the previous step.
	filter.dmr() function

Step 5: Calculate DMR candidates. Result is an object of class GRanges
	edmr() function
