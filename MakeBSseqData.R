###### make an object of BSseq given count data from several replicates
## The input is a list of  data frames with columns: chr, pos, N, X.
makeBSseqData <- function(dat, sampleNames) {
  n0 <- length(dat)                                                              ##n0 is the number of replicates

  if(missing(sampleNames))
    sampleNames <- paste("sample", 1:n0, sep="")                                 ##Gives idenitifiers to samples which don't have names

  alldat <- dat[[1]]
  if(any(alldat[,"N"] < alldat[,"X"], na.rm=TRUE))                               ##Any is the function which checks if any condition is true
    stop("Some methylation counts are greater than coverage.\n")                 ##Checks whether the counts match the coverage

  ix.X <- which(colnames(alldat) == "X")
  ix.N <- which(colnames(alldat) == "N")
  colnames(alldat)[ix.X] <- "X1"
  colnames(alldat)[ix.N] <- "N1"

  if(n0 > 1) {                                                                   ##Function within a function which merges multiple replicates
    for(i in 2:n0) {
      thisdat <- dat[[i]]
      if(any(thisdat[,"N"] < thisdat[,"X"], na.rm=TRUE))
        stop("Some methylation counts are greater than coverage.\n")

      ix.X <- which(colnames(thisdat) == "X")
      ix.N <- which(colnames(thisdat) == "N")
      colnames(thisdat)[c(ix.X,ix.N)] <- paste(c("X", "N"),i, sep="")
      alldat <- merge(alldat, thisdat, all=TRUE)
    }
  }

  ## make BSseq object
  ix.X <- grep("X", colnames(alldat))
  ix.N <- grep("N", colnames(alldat))
  alldat[is.na(alldat)] <- 0
  M <- as.matrix(alldat[,ix.X, drop=FALSE])
  Cov <- as.matrix(alldat[,ix.N, drop=FALSE])
  colnames(M) <- colnames(Cov) <- sampleNames

  ## order CG sites according to positions
  idx <- split(1:length(alldat$chr), alldat$chr)
  M.ordered <- M
  Cov.ordered <- Cov
  pos.ordered <- alldat$pos

  for( i in seq(along=idx) ) {
    thisidx = idx[[i]]
    thispos = alldat$pos[ thisidx ]
    dd = diff(thispos)
    if( min(dd)<0 ) { # not ordered
      warning( paste0("CG positions in chromosome ",  names(idx)[i], " is not ordered. Reorder CG sites.\n") )
      iii = order(thispos)
      M.ordered[thisidx, ] <- M[thisidx, ][iii,]
      Cov.ordered[thisidx, ] <- Cov[thisidx, ][iii,]
      pos.ordered[thisidx] <- alldat$pos[thisidx][iii]
    }
  }

  result <- BSseq(chr=alldat$chr, pos=pos.ordered, M=M.ordered, Cov=Cov.ordered)

  ##    result <- BSseq(chr=alldat$chr, pos=alldat$pos, M=M, Cov=Cov)

  result
}
