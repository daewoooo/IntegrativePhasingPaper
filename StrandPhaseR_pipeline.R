
require('StrandPhaseR')
require('BSgenome.Hsapiens.UCSC.hg19')
strandPhaseR(inputfolder=<StrandSeq_bams> , outputfolder=<user_defined> , numCPU = 1, WCregions = <WCregions_list> , positions='Platinum_NA12878_SNP_allChroms', pairedEndReads=TRUE, chromosomes=c(1:22), min.mapq=10, min.baseq=20, num.iterations=2, translateBases=TRUE, fillMissAllele='NA12878_merged.bam', splitPhasedReads=TRUE, exportVCF='NA12878',callBreaks=FALSE, bsGenome=BSgenome.Hsapiens.UCSC.hg19)
