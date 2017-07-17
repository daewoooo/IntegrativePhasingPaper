import pysam
from os.path import splitext

"""
R shell command to install StrandPhaseR: 
> library(devtools)
> withr::with_libpaths(new = "~/", install_git("git://github.com/daewoooo/StrandPhaseR.git", branch = "master"))
Download Strand-Seq BAMs in the folder StrandS_BAMs using link: https://zenodo.org/record/583682#.WVusQPF95hG

# Install Miniconda
RUN echo 'export PATH=/opt/miniconda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH /opt/miniconda/bin:$PATH
RUN conda config --add channels bioconda --add channels r --add channels conda-forge

#TODO: update whatshap version
# Install all dependencies available as conda packages
RUN conda install -y python=3.5.2 snakemake=3.7.1 samtools=1.2 picard=1.126 \
        vcftools=0.1.14 bwa=0.7.12 whatshap=0.13 \
        biopython=1.68 htslib=1.4 bcftools=1.5
"""
strandseq_bams = ['NW130711_291.bam','NW130711_293.bam','NW130711_295.bam','NW130711_296.bam','NW130711_300.bam','NW130711_302.bam','NW130711_304.bam','NW130711_307.bam','NW130711_317.bam','NW130711_319.bam','NW130711_320.bam','NW130711_321.bam','NW130711_322.bam','NW130711_323.bam','NW130711_324.bam','NW130711_325.bam','NW130711_327.bam','NW130711_329.bam','NW130711_331.bam','NW130711_332.bam','NW130711_334.bam','NW130711_345.bam','NW130711_350.bam','NW130711_354.bam','NW130711_355.bam','NW130711_362.bam','NW130711_365.bam','NW130711_369.bam','NW130711_370.bam','NW130711_371.bam','NW130711_373.bam','NW130711_374.bam','NW130711_376.bam','NW130711_377.bam','NW130711_382.bam','NW150212-III_03.bam','NW150212-III_04.bam','NW150212-III_05.bam','NW150212-III_06.bam','NW150212-III_07.bam','NW150212-III_08.bam','NW150212-III_09.bam','NW150212-III_11.bam','NW150212-III_12.bam','NW150212-III_13.bam','NW150212-III_16.bam','NW150212-III_17.bam','NW150212-III_18.bam','NW150212-III_21.bam','NW150212-III_23.bam','NW150212-III_25.bam','NW150212-III_26.bam','NW150212-III_27.bam','NW150212-III_28.bam','NW150212-III_29.bam','NW150212-III_30.bam','NW150212-III_32.bam','NW150212-III_33.bam','NW150212-III_35.bam','NW150212-III_36.bam','NW150212-III_39.bam','NW150212-III_42.bam','NW150212-III_45.bam','NW150212-III_48.bam','NW150212-III_51.bam','NW150212-III_52.bam','NW150212-III_55.bam','NW150212-III_57.bam','NW150212-III_59.bam','NW150212-III_63.bam','NW150212-III_66.bam','NW150212-III_68.bam','NW150212-III_69.bam','NW150212-III_75.bam','NW150212-III_76.bam','NW150212-III_77.bam','NW150212-III_80.bam','NW150212-III_81.bam','NW150212-III_84.bam','NW150212-III_85.bam','NW150212-III_87.bam','NW150212-III_90.bam','NW150212-III_92.bam','NW150212-III_93.bam','NW150212-III_94.bam','NW150212-IV_04.bam','NW150212-IV_08.bam','NW150212-IV_09.bam','NW150212-IV_10.bam','NW150212-IV_11.bam','NW150212-IV_12.bam','NW150212-IV_15.bam','NW150212-IV_17.bam','NW150212-IV_21.bam','NW150212-IV_24.bam','NW150212-IV_25.bam','NW150212-IV_26.bam','NW150212-IV_28.bam','NW150212-IV_31.bam','NW150212-IV_34.bam','NW150212-IV_39.bam','NW150212-IV_41.bam','NW150212-IV_43.bam','NW150212-IV_44.bam','NW150212-IV_45.bam','NW150212-IV_46.bam','NW150212-IV_47.bam','NW150212-IV_50.bam','NW150212-IV_53.bam','NW150212-IV_56.bam','NW150212-IV_58.bam','NW150212-IV_60.bam','NW150212-IV_62.bam','NW150212-IV_63.bam','NW150212-IV_64.bam','NW150212-IV_65.bam','NW150212-IV_66.bam','NW150212-IV_68.bam','NW150212-IV_69.bam','NW150212-IV_72.bam','NW150212-IV_73.bam','NW150212-IV_74.bam','NW150212-IV_75.bam','NW150212-IV_77.bam','NW150212-IV_79.bam','NW150212-IV_81.bam','NW150212-IV_83.bam','NW150212-IV_84.bam','NW150212-IV_85.bam','NW150212-IV_87.bam','NW150212-IV_89.bam','NW150212-IV_91.bam','NW150212-IV_92.bam','NW150212-IV_93.bam']
#tools
picard = 'picard'
samtools = 'samtools'
whatshap = 'whatshap'
bcftools = 'bcftools'
# PATH to the directory where StrandPhaseR is installed
strandphaser = '~/'

# parameters
coverage = [2]
chromosome = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
strandseqcoverage= [134]
trials=[3]


rule master:
	input: 'summary.eval'
	message: 'MASTER rule'

rule download_SS:
	output: 'download_ss/'
	run:
		shell("git clone https://github.com/daewoooo/IntegrativePhasingPaper.git {output}")

rule create_dummy_file:
	output:'download_ss/StrandS_suppData/TRIAL{trials,[0-9]+}_downsampled/WCregions/NA12878_WC_regions_hg19_0cellsSample',
	shell: 'touch {output}'

rule download_strandseq_bams:
	threads: 10
	output: 'StrandS_BAMs/{file,NW.*}.bam'
	log: 'StrandS_BAMs/{file}.log'
	shell: 'wget -O {output} -o {log} https://zenodo.org/record/830278/files/{wildcards.file}.bam'

rule merge_strandseq_bams:
	input:
		bams=expand('StrandS_BAMs/{bam}', bam=strandseq_bams),
	output:
		bam='StrandS_BAMs/NA12878_merged.bam',
	shell:
		'samtools merge {output.bam} {input.bams}'

rule run_SS_pipeline:
	input:
		strandseqbams=expand('StrandS_BAMs/{bam}', bam=strandseq_bams),
		dontknow= 'download_ss/',
		mergedbam='StrandS_BAMs/NA12878_merged.bam',
	output:
		vcf='StrandPhaseR_TRIAL_{trials,[0-9]+}_{strandseqcoverage,[0-9]+}cells/VCFfiles/chr{chromosome,[0-9]+}_phased.vcf', 
	run:
		if wildcards.strandseqcoverage=='0':
			shell('awk \'($0 ~ /^#/)\' {input.unphasedvcf} > {output.vcf}')
		else:
			shell('Rscript download_ss/StrandS_suppData/StrandPhaseR_pipeline.R StrandS_BAMs StrandPhaseR_TRIAL_{wildcards.trials}_{wildcards.strandseqcoverage}cells download_ss/StrandS_suppData/TRIAL{wildcards.trials}_downsampled/WCregions/NA12878_WC_regions_hg19_{wildcards.strandseqcoverage}cellsSample download_ss/StrandS_suppData/Platinum_NA12878_SNP_allChroms.txt {wildcards.chromosome} {strandphaser} {input.mergedbam}')

rule download_pacbio:
	output:
		protected("bam/NA12878.pacbio.chrall.covall.{ext,bam}")
	shell:
		"wget -O {output} ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.{wildcards.ext}"

rule download_illumina:
	output:
		protected("download/NA12878.illumina.chrall.covall.{ext,bam}")
	shell:
		"wget -O {output} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/high_coverage_alignment/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.{wildcards.ext}"

rule index_illumina:
	input: 'download/NA12878.illumina.chrall.covall.bam'
	output: 'download/NA12878.illumina.chrall.covall.bam.bai'
	shell: 'samtools index {input}'

rule download_reference:
	output:
		'reference/human_g1k_v37.fasta'
	shell:
		"""
		wget -O {output}.gz.incomplete ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
		mv {output}.gz.incomplete {output}.gz
		# gunzip fails with "decompression OK, trailing garbage ignored" because
		# the file is razf-compressed (gzip-compatible, but with an index at end)
		gunzip -f {output}.gz || true
		# cd reference && md5sum -c MD5SUM
		"""

rule update_ref:
	input: 'reference/human_g1k_v37.fasta'
	output: 'reference/human_g1k_v37.notation.fasta'
	shell: 'sed -e \'s/>/>chr/g\' {input} > {output}'


rule download_10xG_vcf:
	output: 
		protected("download/NA12878.10xG.phased.chrall.vcf.gz")
	shell:
		"wget -O {output} https://s3-us-west-2.amazonaws.com/10x.files/samples/genome/NA12878_WGS_210/NA12878_WGS_210_phased_variants.vcf.gz"

rule reheader_illumina:
	input:'download/NA12878.illumina.chrall.covall.bam', 'download/NA12878.illumina.chrall.covall.bam.bai'
	output: 'bam/NA12878.illumina.chrall.reheader.covall.sam'
	shell: "samtools view -h {input} | sed -e \'/^@SQ/s/SN:/SN:chr/\' -e \'/^\[^@\]/s/\t/\tchr/2\'  > {output}"

rule filter_illumina:
	input:'download/NA12878.illumina.chrall.covall.bam', 'download/NA12878.illumina.chrall.covall.bam.bai', 'bam/NA12878.illumina.chrall.reheader.covall.sam'
	output:'bam/NA12878.illumina.chrall.covall.bam'
	shell: 'samtools reheader {input[2]} {input[0]} > {output}'

rule index_bams:
	input:'bam/NA12878.{data, (illumina|pacbio|10xG)}.chrall.covall.bam'
	output: 'bam/NA12878.{data, (illumina|pacbio|10xG)}.chrall.covall.bam.bai'
	shell: 'samtools index {input}'

rule consensus_VCFs:
	threads: 100
	output:
		protected("download/NA12878.benchmark.phased.chrall.vcf.gz")
	shell:
		"wget -O {output} ftp://platgene_ro:""@ussd-ftp.illumina.com/2016-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz"

rule filter_10xG_vcf:
	input: 'download/NA12878.10xG.phased.chrall.vcf.gz',
	output: 'vcf/NA12878.10xG.phased.chrall.vcf'
	shell: 'zcat {input} | awk -vOFS=\'\\t\' \'($0 ~ /^#/) || (($7=".") && ($8="."))\' -  > {output}'

rule unzip_vcf:
	input: 'download/NA12878.benchmark.phased.chrall.vcf.gz'
	output: 'vcf/NA12878.benchmark.phased.chrall.vcf'
	shell: 'zcat {input} > {output}'

rule unphase_vcfs:
	input: 'vcf/NA12878.{data,(benchmark|10xG)}.phased.chrall.vcf',
	output: 'vcf/NA12878.{data,(benchmark|10xG)}.unphased.chrall.vcf',
	shell: '{whatshap} unphase {input} > {output}'

rule split_vcf:
	input: 'vcf/NA12878.{data,(benchmark|10xG)}.{type,(unphased|phased)}.chrall.vcf',
	output: 'vcf/NA12878.{data,(benchmark|10xG)}.{type,(unphased|phased)}.chr{chromosome,[0-9]+}.vcf',
	message: 'Extracting chromosome {wildcards.chromosome} from {input}'
	shell: """awk '/^#/ || ($1 == "chr{wildcards.chromosome}")' {input} > {output}"""

rule extract_chromosome:
	input:
		bam='bam/NA12878.{data, (illumina|pacbio|10xG)}.chrall.covall.bam',
		bai='bam/NA12878.{data, (illumina|pacbio|10xG)}.chrall.covall.bam.bai'
	output:
		bam='bam/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.covall.bam'
	log: 'bam/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.covall.bam.log'
	shell: '(samtools view -h {input.bam} chr{wildcards.chromosome} | samtools view -Sb - > {output.bam}) 2>{log}'


rule calculate_coverage:
	input: 'bam/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.covall.bam'
	output: 'bam/stats/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.covall.bam.coverage'
	message: 'Computing coverage for {input}'
	run: 
		bam = pysam.Samfile(input[0])
		length = None
		for e in bam.header.get('SQ'):
			if not any(c.isalpha() for c in e['SN'][3:]):
				if e['SN'][3:] == wildcards.chromosome:
					length = e['LN']
		assert length != None
		shell("samtools depth {input} | awk '{{sum+=$3}} END {{ print sum/{length} }}' > {output}")

rule downsampling:
	input:
		bam='bam/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.covall.bam',
		coverage='bam/stats/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.covall.bam.coverage'
	output: 
		bam='bam/TRIAL-{trials,[0-9]+}/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}.bam',
		bai='bam/TRIAL-{trials,[0-9]+}/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}.bai'
	log: 'bam/TRIAL-{trials,[0-9]+}/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}.bam.log'
	message: 'Downsampling {input} to {wildcards.coverage}x'
	run:
		input_coverage = float(open(input.coverage).readline())
		p=1
		if wildcards.coverage!='all':	
			p = float(wildcards.coverage) / input_coverage
		seed = hash(output)
		shell('{picard} DownsampleSam INPUT={input.bam} RANDOM_SEED=null CREATE_INDEX=true OUTPUT={output.bam} PROBABILITY= {p} VALIDATION_STRINGENCY=SILENT > {log} 2>&1')

rule add_read_groups:
	input: 'bam/TRIAL-{trials,[0-9]+}/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}.bam'
	output: 
		bam='bam/TRIAL-{trials,[0-9]+}/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}.readgroup.bam',
		bai='bam/TRIAL-{trials,[0-9]+}/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}.readgroup.bai'
	log: 'bam/TRIAL-{trials,[0-9]+}/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}.readgroup.bam.log'
	message: 'Adding read groups'
	run: 
		shell('{picard} AddOrReplaceReadGroups CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT I={input} O={output.bam} ID=NA12878 LB=UNKNOWN PL={wildcards.data} PU=UNKNOWN SM=NA12878 > {log} 2>&1')

rule sort_bams:
	input: 'bam/TRIAL-{trials,[0-9]+}/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}.readgroup.bam'
	output: 'bam/TRIAL-{trials,[0-9]+}/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}.readgroup.sorted.bam', 'bam/TRIAL-{trials,[0-9]+}/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}.readgroup.sorted.bai', 'bam/TRIAL-{trials,[0-9]+}/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}.readgroup.sorted.bam.md5'
	log: 'bam/TRIAL-{trials,[0-9]+}/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}.readgroup.sorted.bam.log'
	message: 'Sorting {input}'
	shell: 'time ({picard} SortSam VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=50000 SORT_ORDER=coordinate CREATE_INDEX=true CREATE_MD5_FILE=true I={input[0]} O={output[0]}) > {log} 2>&1'


# run different combinations
rule integrative_whatshap_pacbio_SS:
	input: 
		bam1='bam/TRIAL-{trials,[0-9]+}/NA12878.pacbio.chr{chromosome,[0-9]+}.cov{pcoverage,(all|[0-9]+)}.readgroup.sorted.bam',
		vcf1='vcf/NA12878.benchmark.unphased.chr{chromosome,[0-9]+}.vcf',
		vcf2='StrandPhaseR_TRIAL_{trials,[0-9]+}_{strandseqcoverage,[0-9]+}cells/VCFfiles/chr{chromosome,[0-9]+}_phased.vcf',
		ref='reference/human_g1k_v37.notation.fasta',
	output: 
		vcf= 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chr{chromosome,[0-9]+}.noindels.vcf',
	log: 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chr{chromosome,[0-9]+}.noindels.vcf.log'
	message: 'Running WHATSHAP on {input.bam1}'
	shell: '{whatshap} phase --max-coverage 15 --distrust-genotypes --sample NA12878 --chromosome chr{wildcards.chromosome} --reference {input.ref} {input.vcf1} {input.vcf2} {input.bam1} --output {output.vcf} 2> {log}'

rule integrative_whatshap_pacbio_SS_indels:
	input: 
		bam1='bam/TRIAL-{trials,[0-9]+}/NA12878.pacbio.chr{chromosome,[0-9]+}.cov{pcoverage,(all|[0-9]+)}.readgroup.sorted.bam',
		vcf1='vcf/NA12878.benchmark.unphased.chr{chromosome,[0-9]+}.vcf',
		vcf2='StrandPhaseR_TRIAL_{trials,[0-9]+}_{strandseqcoverage,[0-9]+}cells/VCFfiles/chr{chromosome,[0-9]+}_phased.vcf',
		ref='reference/human_g1k_v37.notation.fasta',
	output: 
		vcf= 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chr{chromosome,[0-9]+}.indels.vcf',
	log: 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chr{chromosome,[0-9]+}.indels.vcf.log'
	message: 'Running WHATSHAP on {input.bam1}'
	shell: '{whatshap} phase --max-coverage 15 --indels --sample NA12878 --chromosome chr{wildcards.chromosome} --reference {input.ref} {input.vcf1} {input.vcf2} {input.bam1} --output {output.vcf} 2> {log}'

rule integrative_whatshap_illumina_SS:
	input: 
		bam1='bam/TRIAL-{trials,[0-9]+}/NA12878.illumina.chr{chromosome,[0-9]+}.cov{icoverage,(all|[0-9]+)}.readgroup.sorted.bam',
		vcf1='vcf/NA12878.benchmark.unphased.chr{chromosome,[0-9]+}.vcf',
		vcf2='StrandPhaseR_TRIAL_{trials,[0-9]+}_{strandseqcoverage,[0-9]+}cells/VCFfiles/chr{chromosome,[0-9]+}_phased.vcf',
		ref='reference/human_g1k_v37.notation.fasta',
	output: 
		vcf= 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio0.illumina{icoverage,(all|[0-9]+)}.10x0.chr{chromosome,[0-9]+}.noindels.vcf',
	log: 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio0.illumina{icoverage,(all|[0-9]+)}.10x0.chr{chromosome,[0-9]+}.noindels.vcf.log'
	message: 'Running WHATSHAP on {input.bam1}'
	shell: '{whatshap} phase --max-coverage 15 --indels --sample NA12878 --chromosome chr{wildcards.chromosome} --reference {input.ref} {input.vcf1} {input.vcf2} {input.bam1} --output {output.vcf} 2> {log}'

rule integrative_whatshap_pacbio_xg:
	input: 
		bam1='bam/TRIAL-{trials,[0-9]+}/NA12878.pacbio.chr{chromosome,[0-9]+}.cov{pcoverage,(all|[0-9]+)}.readgroup.sorted.bam',
		vcf1='vcf/NA12878.benchmark.unphased.chr{chromosome,[0-9]+}.vcf',
		vcf2='vcf/NA12878.10xG.phased.chr{chromosome,[0-9]+}.vcf',
		ref='reference/human_g1k_v37.notation.fasta',
	output: 
		vcf= 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x{xcoverage,(all|[0-9]+)}.chr{chromosome,[0-9]+}.noindels.vcf',
	log: 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x{xcoverage,(all|[0-9]+)}.chr{chromosome,[0-9]+}.noindels.vcf.log'
	message: 'Running WHATSHAP on {input.bam1}'
	shell: '{whatshap} phase --max-coverage 15 --distrust-genotypes --chromosome chr{wildcards.chromosome} --sample NA12878 --reference {input.ref} {input.vcf1} {input.vcf2} {input.bam1} --output {output.vcf} 2> {log}'

rule integrative_whatshap_pacbio_only:
	input: 
		bam1='bam/TRIAL-{trials,[0-9]+}/NA12878.pacbio.chr{chromosome,[0-9]+}.cov{pcoverage,(all|[0-9]+)}.readgroup.sorted.bam',
		vcf1='vcf/NA12878.benchmark.unphased.chr{chromosome,[0-9]+}.vcf',
		ref='reference/human_g1k_v37.notation.fasta',
	output: 
		vcf= 'whatshap_pacbio_only/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chr{chromosome,[0-9]+}.noindels.vcf',
	log: 'whatshap_pacbio_only/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chr{chromosome,[0-9]+}.noindels.vcf.log'
	message: 'Running WHATSHAP on {input.bam1}'
	shell: '{whatshap} phase --max-coverage 15 --distrust-genotypes --chromosome chr{wildcards.chromosome} --sample NA12878 --reference {input.ref} {input.vcf1} {input.bam1} --output {output.vcf} 2> {log}'

rule integrative_whatshap_SS_xg:
	input: 
		vcf3='vcf/NA12878.10xG.phased.chr{chromosome,[0-9]+}.vcf',
		vcf1='vcf/NA12878.benchmark.unphased.chr{chromosome,[0-9]+}.vcf',
		vcf2='StrandPhaseR_TRIAL_{trials,[0-9]+}_{strandseqcoverage,[0-9]+}cells/VCFfiles/chr{chromosome,[0-9]+}_phased.vcf',
	output: 
		vcf= 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio0.illumina0.10x{xcoverage,(all|[0-9]+)}.chr{chromosome,[0-9]+}.noindels.vcf',
	log: 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio0.illumina0.10x{xcoverage,(all|[0-9]+)}.chr{chromosome,[0-9]+}.noindels.vcf.log'
	shell: '{whatshap} phase --max-coverage 15 --distrust-genotypes --sample NA12878 {input.vcf1} {input.vcf2} {input.vcf3} --output {output.vcf} 2> {log}'

rule merge_vcfs_pacbio_SS:
	input:expand('whatshap_integrative_phasing/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio{pcoverage}.illumina0.10x0.chr{chromosome}.noindels.vcf', chromosome=chromosome, trials=trials, pcoverage=coverage, strandseqcoverage=strandseqcoverage), 
	output:'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chrall.noindels.vcf', 
	shell: '({bcftools} concat {input} | vcf-sort > {output[0]})'

rule merge_vcfs_illumina_SS:
	input:expand('whatshap_integrative_phasing/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio0.illumina{icoverage}.10x0.chr{chromosome}.noindels.vcf', chromosome=chromosome, trials=trials, icoverage =coverage, strandseqcoverage=strandseqcoverage),
	output:'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio0.illumina{icoverage,(all|[0-9]+)}.10x0.chrall.noindels.vcf',
	shell: '({bcftools} concat {input} | vcf-sort > {output[0]})'

rule merge_vcfs_SS_xg:
	input:expand('whatshap_integrative_phasing/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio0.illumina0.10xall.chr{chromosome}.noindels.vcf', chromosome=chromosome, trials=trials, strandseqcoverage=strandseqcoverage),
	output:'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio0.illumina0.10xall.chrall.noindels.vcf',
	shell: '({bcftools} concat {input} | vcf-sort > {output[0]})'

rule merge_vcfs_pacbio_xg:
	input:expand('whatshap_integrative_phasing/TRIAL-{trials}/strandseqcells0.pacbio{pcoverage}.illumina0.10xall.chr{chromosome}.noindels.vcf', chromosome=chromosome, trials=trials, pcoverage=coverage),
	output:'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10xall.chrall.noindels.vcf',
	shell: '({bcftools} concat {input} | vcf-sort > {output[0]})'

rule merge_vcfs_pacbio_only:
	input: expand('whatshap_pacbio_only/TRIAL-{trials}/strandseqcells0.pacbio{pcoverage}.illumina0.10x0.chr{chromosome}.noindels.vcf',chromosome=chromosome, trials=trials, pcoverage=coverage), 
	output: 'whatshap_pacbio_only/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chrall.noindels.vcf'
	shell: '({bcftools} concat {input} | vcf-sort > {output[0]})'

# bcftools error: The sequence "chr2" not defined in the header of StrandSeq VCF files. Not running evaluation for it.
rule merge_vcfs_SS_only:
	input: expand('StrandPhaseR_TRIAL_{trials}_{strandseqcoverage}cells/VCFfiles/chr{chromosome}_phased.vcf', trials=trials, strandseqcoverage=strandseqcoverage, chromosome=chromosome)
	output: 'StrandPhaseR_TRIAL_{trials,[0-9]+}_{strandseqcoverage,[0-9]+}cells/VCFfiles/chrall_phased.vcf'
	shell: '({bcftools} concat {input} | vcf-sort > {output[0]})'

rule evaluate_whatshap:
	input:
		truth='vcf/NA12878.benchmark.phased.chrall.vcf',
		phased='whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina{icoverage,(all|[0-9]+)}.10x{xcoverage,(all|[0-9]+)}.chrall.noindels.vcf',
	output:	'eval/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina{icoverage,(all|[0-9]+)}.10x{xcoverage,(all|[0-9]+)}.chrall.noindels.eval',
	log: 'eval/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina{icoverage,(all|[0-9]+)}.10x{xcoverage,(all|[0-9]+)}.chrall.noindels.log',
	shell: '{whatshap} compare --names benchmark,whatshap --tsv-pairwise {output} --only-snvs {input.truth} {input.phased} > {log} 2>&1'

rule evaluate_whatshap_indels:
	input:
		truth='vcf/NA12878.benchmark.phased.chrall.vcf',
		phased='whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chrall.indels.vcf',
	output:	'eval/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina{icoverage,(all|[0-9]+)}.10x{xcoverage,(all|[0-9]+)}.chrall.indels.eval',
	log: 'eval/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina{icoverage,(all|[0-9]+)}.10x{xcoverage,(all|[0-9]+)}.chrall.indels.log',
	shell: '{whatshap} compare --names benchmark,whatshap --tsv-pairwise {output} {input.truth} {input.phased} > {log} 2>&1'

rule evaluate_whatshap_pacbio_only:
	input:
		truth='vcf/NA12878.benchmark.phased.chrall.vcf',
		phased='whatshap_pacbio_only/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chrall.noindels.vcf',
	output:	'eval/whatshap_pacbio_only/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chrall.noindels.eval',
	log: 'eval/whatshap_pacbio_only/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chrall.noindels.log',
	shell: '{whatshap} compare --names benchmark,whatshap --only-snvs --tsv-pairwise {output} {input.truth} {input.phased} > {log} 2>&1'

rule evaluate_whatshap_10xg_only:
	input:
		truth='vcf/NA12878.benchmark.phased.chrall.vcf',
		phased='vcf/NA12878.10xG.phased.chrall.vcf',
	output:	'eval/whatshap_10xG_only/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio0.illumina0.10xall.chrall.noindels.eval',
	log: 'eval/whatshap_10xG_only/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio0.illumina0.10xall.chrall.noindels.log',
	shell: '{whatshap} compare --names benchmark,whatshap --only-snvs --tsv-pairwise {output} {input.truth} {input.phased} > {log} 2>&1'


rule evaluate_whatshap_SS_only:
	input:
		truth='vcf/NA12878.benchmark.phased.chrall.vcf',
		phased='StrandPhaseR_TRIAL_{trials,[0-9]+}_{strandseqcoverage,[0-9]+}cells/VCFfiles/chrall_phased.vcf',
	output:	'eval/whatshap_SS_only/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio0.illumina0.10x0.chrall.noindels.eval',
	log: 'eval/whatshap_SS_only/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio0.illumina0.10x0.chrall.noindels.log',
	shell: '{whatshap} compare --names benchmark,whatshap --only-snvs --tsv-pairwise {output} {input.truth} {input.phased} > {log} 2>&1'

rule summary:
	output:'summary.eval',
	input:
		expand('eval/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio{pcoverage}.illumina0.10x0.chrall.noindels.eval', strandseqcoverage=strandseqcoverage, pcoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio0.illumina{icoverage}.10x0.chrall.noindels.eval', strandseqcoverage=strandseqcoverage, icoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio0.illumina0.10x{xcoverage}.chrall.noindels.eval', strandseqcoverage=strandseqcoverage, xcoverage=['all'], trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells0.pacbio{pcoverage}.illumina0.10x{xcoverage}.chrall.noindels.eval', xcoverage=['all'], pcoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells0.pacbio0.illumina{icoverage}.10x0.chrall.noindels.eval',  icoverage=coverage, trials=trials),
		expand('eval/whatshap_pacbio_only/TRIAL-{trials}/strandseqcells0.pacbio{pcoverage}.illumina0.10x0.chrall.noindels.eval',  pcoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells0.pacbio0.illumina{icoverage}.10x0.chrall.noindels.eval', icoverage=coverage, trials=trials),
		expand('eval/whatshap_10xG_only/TRIAL-{trials}/strandseqcells0.pacbio0.illumina0.10xall.chrall.noindels.eval',  trials=trials),
		#expand('eval/whatshap_SS_only/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio0.illumina0.10x0.chrall.noindels.eval', trials=trials, strandseqcoverage=strandseqcoverage),
	message: 'Aggregating statistics to {output}'
	run:
		first = input[0]
		rest = ' '.join(input[1:]) 
		shell('(cat {first} && ( cat {rest} | grep -v \'^#\' ) ) > {output}')
