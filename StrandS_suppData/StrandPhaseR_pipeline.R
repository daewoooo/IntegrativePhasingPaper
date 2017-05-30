#!/usr/bin/Rscript

args=commandArgs(TRUE)

#add user defined path to load needed libraries
.libPaths( c( .libPaths(), args[6]) )

suppressPackageStartupMessages(library(StrandPhaseR))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))

strandPhaseR(inputfolder=args[1], outputfolder=args[2], WCregions = args[3] , positions=args[4], chromosomes=args[5], numCPU = 1, pairedEndReads=TRUE, min.mapq=10, min.baseq=20, num.iterations=2, translateBases=TRUE, fillMissAllele=args[7], splitPhasedReads=FALSE, exportVCF='NA12878',callBreaks=FALSE, bsGenome=BSgenome.Hsapiens.UCSC.hg19)
