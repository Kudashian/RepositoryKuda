Input file is a methylation call file from aligner pipeline or generic text file with aligner output format:
chrBase chr base strand coverage freqC freqT

Step 1: Create a list containing all the input files which will be fed into the script.
	list() function.

Tabix is the machine-readable version of an input file. Tabix file is usually indexed zipped files.
If memory is not enough then methylKit objects can be saved as tabix files.

Step 2: All samples need to be merged in one object for comparison analysis.
	Unite() function.
		min.per.group - integer which denotes minimum no. of samples per replicate  which need to cover a base/region. If set to 2, then the bases covered by at least 2 replicates will be united.
	(Missing data appears as NAs)
		save.db - determines whether resultant object must be saved

Step 3: Correlate and perform PCA of grouped samples
	Correlation() function.
		method - string which indicated method
 	PCASamples() function
	AssocComp() function

Step 4: Tiling window analysis
	Summarizes methylation information over tiling windows vs base-pair resolution analysis.
	TileMethylCounts() functions
		win.size
		step.size        **We won't know where exactly each tile is.
	These tiles can be loaded into Unite & CalculateDiffMeth functions

Step 5: Finding differentially methylated bases/regions
	CalculateDiffMeth() function
	Fisher's exact (single samples) & logistic regression (multiple replicates) calculate P-values
	You can pool() all samples to one test/control to force Fisher's test
	P-values adjusted to Q-values using SLIM method

	Covariates ?
	Chi-squared tests used as default
	Correct for overdispersion
		Variability then assumed by the distribution

Step 6: Annotation
	"Genomation" package
	Read in gene BED file
	ReadGeneric() function
		inputs genomic data and converts it to a GRanges object

	
