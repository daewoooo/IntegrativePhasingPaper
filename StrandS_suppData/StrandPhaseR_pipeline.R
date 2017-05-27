#!/usr/bin/Rscript

require('StrandPhaseR')
require('BSgenome.Hsapiens.UCSC.hg19')

args=commandArgs(TRUE)

strandPhaseR(inputfolder=args[1], outputfolder=args[2], WCregions = args[3] , positions=args[4], numCPU = 1, pairedEndReads=TRUE, chromosomes=c(1:22), min.mapq=10, min.baseq=20, num.iterations=2, translateBases=TRUE, fillMissAllele=NULL, splitPhasedReads=FALSE, exportVCF='NA12878',callBreaks=FALSE, bsGenome=BSgenome.Hsapiens.UCSC.hg19)
