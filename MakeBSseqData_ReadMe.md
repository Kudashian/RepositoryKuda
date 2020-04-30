Utility function for creating BSseq object from multiple data frames.

After sequence alignment and proper processing, the BS-seq data can be
summarized by following information at each C position (mostly CpG sites, with
some CH): chromosome number, genomic coordinate, total number of reads covering
the position, and number of reads showing methylation at this position.

For replicated samples, the data need to be merged based on the chromosome
number and genomic coordinates. This function provide such functionality.
It takes replicated data as a list of data frames, merged them
and create a BSseq object.

makeBSseqData(dat, sampleNames)

dat: A list of multiple data frames from biological replicates. Each element
represents data from one replicate. The data frame MUST contain following
columns in correct order: (1) Chromosome number; (2) Genomic coordinates;
(3) Read coverage of the position from BS-seq data; (4) Number of reads
showing methylation of the position.
The colnames MUST BE "chr", "pos", "N", "X".
