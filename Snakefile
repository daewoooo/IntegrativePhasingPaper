import pysam
from os.path import splitext

"""
Dependencies: samtools, picard, whatshap, and all these binaries in tools directory.
whatshap branch: read_selection_fix
StrandPhaseR: withr::with_libpaths(new = "tools/StrandPhaseR", install_git("git://github.com/daewoooo/StrandPhaseR.git", branch = "master"))
"""

#tools
picard = 'tools/picard'
samtools = 'tools/samtools'
whatshap = 'tools/whatshap'
reference = 'reference/human_g1k_v37.notation.fasta'
# PATH to the directory where StrandPhaseR is installed
strandphaser = 'tools/StrandPhaseR'

# parameters
coverage = [2,3,4,5,10,15,25,30]
chromosome = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
strandseqcoverage= [5,10,20,40,60,80,100,120,134]
trials=[2,3,4,5,6]

rule master:
	input: 'summary.eval'
	message: 'MASTER rule'

rule download_SS:
	threads: 100
	output: 'download/'
	run:
		shell("git clone https://github.com/daewoooo/IntegrativePhasingPaper.git {output}")

rule create_dummy_file:
	output:'download/StrandS_suppData/TRIAL{trials,[0-9]+}_downsampled/WCregions/NA12878_WC_regions_hg19_0cellsSample',
	shell: 'touch {output}'

#TODO: add rule to download SS BAMs
rule run_SS_pipeline:
	input: '../StrandS_BAMs', 'download/', 'vcf/NA12878.benchmark.unphased.chr{chromosome,[0-9]+}.vcf', '../StrandS_BAMs/NA12878_merged.bam'
	output:'StrandPhaseR_TRIAL_{trials,[0-9]+}_{strandseqcoverage,[0-9]+}cells/VCFfiles/chr{chromosome,[0-9]+}_phased.vcf', 
	run: 
		if wildcards.strandseqcoverage=='0':
			shell('awk \'($0 ~ /^#/)\' {input[2]} > {output}')
		else:
			shell('Rscript download/StrandS_suppData/StrandPhaseR_pipeline.R {input[0]} StrandPhaseR_TRIAL_{wildcards.trials}_{wildcards.strandseqcoverage}cells download/StrandS_suppData/TRIAL{wildcards.trials}_downsampled/WCregions/NA12878_WC_regions_hg19_{wildcards.strandseqcoverage}cellsSample download/StrandS_suppData/Platinum_NA12878_SNP_allChroms.txt {wildcards.chromosome} {strandphaser} {input[3]}')


rule download_pacbio:
	threads: 100
	output:
		protected("bam/NA12878.pacbio.chrall.{ext,bam}")
	shell:
		"wget -O {output} ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.{wildcards.ext}"

rule download_illumina:
	threads: 100
	output:
		protected("download/NA12878.illumina.chrall.{ext,bam}")
	shell:
		"wget -O {output} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/high_coverage_alignment/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.{wildcards.ext}"

rule index_illumina:
	input: 'download/NA12878.illumina.chrall.bam'
	output: 'download/NA12878.illumina.chrall.bam.bai'
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

rule download_10xG:
	output:
		protected("download/NA12878.10xG.chrall.{ext,(bam|bam.bai)}")
	shell:
		"wget -O {output} https://s3-us-west-2.amazonaws.com/10x.files/samples/genome/NA12878_WGS_210/NA12878_WGS_210_phased_possorted_bam.{wildcards.ext}"

rule download_10xG_vcf:
	output: 
		protected("download/NA12878.10xG.phased.chrall.vcf.gz")
	shell:
		"wget -O {output} https://s3-us-west-2.amazonaws.com/10x.files/samples/genome/NA12878_WGS_210/NA12878_WGS_210_phased_variants.vcf.gz"

rule reheader_illumina:
	input:'download/NA12878.illumina.chrall.bam', 'download/NA12878.illumina.chrall.bam.bai'
	output: 'bam/NA12878.illumina.chrall.reheader.sam'
	shell: "samtools view -h {input} | sed -e \'/^@SQ/s/SN:/SN:chr/\' -e \'/^\[^@\]/s/\t/\tchr/2\'  > {output}"

rule filter_illumina:
	input:'download/NA12878.illumina.chrall.bam', 'download/NA12878.illumina.chrall.bam.bai', 'bam/NA12878.illumina.chrall.reheader.sam'
	output:'bam/NA12878.illumina.chrall.bam'
	shell: 'samtools reheader {input[2]} {input[0]} > {output}'

rule filter_10xG:
	input: 'download/NA12878.10xG.chrall.bam',
	output: 'bam/NA12878.10xG.chrall.bam'
	shell: 'mv {input} {output}'

rule index_bams:
	input:'bam/NA12878.{data, (illumina|pacbio|10xG)}.chrall.bam'
	output: 'bam/NA12878.{data, (illumina|pacbio|10xG)}.chrall.bam.bai'
	shell: 'samtools index {input}'

#TODO: update the link
rule consensus_VCFs:
	threads: 100
	output:
		protected("download/NA12878.benchmark.phased.chrall.vcf.gz")
	shell:
		"wget -O {output} ftp://platgene_ro:""@ussd-ftp.illumina.com/2016-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz"

rule unzip_vcf:
	input: 'download/NA12878.benchmark.phased.chrall.vcf.gz'
	output: 'vcf/NA12878.benchmark.phased.chrall.vcf'
	shell: 'zcat {input} > {output}'

rule filter_10xG_vcf:
	input: 'download/NA12878.10xG.phased.chrall.vcf.gz',
	output: 'vcf/NA12878.10xG.phased.chrall.vcf'
	shell: 'zcat {input} | awk -vOFS=\'\\t\' \'($0 ~ /^#/) || (($7=".") && ($8="."))\' -  > {output}'

rule unphase_vcfs:
	input: 'vcf/NA12878.{data,(benchmark|10xG)}.phased.chrall.vcf',
	output: 'vcf/NA12878.{data,(benchmark|10xG)}.unphased.chrall.vcf'
	shell: '{whatshap} unphase {input} > {output}'

rule split_vcf:
	input: 'vcf/NA12878.{data,(benchmark|10xG)}.{type,(unphased|phased)}.chrall.vcf'
	output: 'vcf/NA12878.{data,(benchmark|10xG)}.{type,(unphased|phased)}.chr{chromosome,[0-9]+}.vcf'
	message: 'Extracting chromosome {wildcards.chromosome} from {input}'
	shell: """awk '/^#/ || ($1 == "chr{wildcards.chromosome}")' {input} > {output}"""

rule extract_chromosome:
	input:
		bam='bam/NA12878.{data, (illumina|pacbio|10xG)}.chrall.bam',
		bai='bam/NA12878.{data, (illumina|pacbio|10xG)}.chrall.bam.bai'
	output:
		bam='bam/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.bam'
	log: 'bam/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.bam.log'
	shell: '(samtools view -h {input.bam} chr{wildcards.chromosome} | samtools view -Sb - > {output.bam}) 2>{log}'


rule calculate_coverage:
	input: 'bam/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.bam'
	output: 'bam/stats/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.bam.coverage'
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
		bam='bam/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.bam',
		coverage='bam/stats/NA12878.{data, (illumina|pacbio|10xG)}.chr{chromosome,[0-9]+}.bam.coverage'
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
		ref=reference,
	output: 
		vcf= 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chr{chromosome,[0-9]+}.noindels.vcf',
	log: 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chr{chromosome,[0-9]+}.noindels.vcf.log'
	message: 'Running WHATSHAP on {input.bam1}'
	shell: '{whatshap} phase --max-coverage 15 --include-VCF-phasing --distrust-genotypes --sample NA12878 --reference {input.ref} {input.vcf1} {input.vcf2} {input.bam1} --output {output.vcf} 2> {log}'

rule integrative_whatshap_pacbio_SS_indels:
	input: 
		bam1='bam/TRIAL-{trials,[0-9]+}/NA12878.pacbio.chr{chromosome,[0-9]+}.cov{pcoverage,(all|[0-9]+)}.readgroup.sorted.bam',
		vcf1='vcf/NA12878.benchmark.unphased.chr{chromosome,[0-9]+}.vcf',
		vcf2='StrandPhaseR_TRIAL_{trials,[0-9]+}_{strandseqcoverage,[0-9]+}cells/VCFfiles/chr{chromosome,[0-9]+}_phased.vcf',
		ref=reference,
	output: 
		vcf= 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chr{chromosome,[0-9]+}.indels.vcf',
	log: 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chr{chromosome,[0-9]+}.indels.vcf.log'
	message: 'Running WHATSHAP on {input.bam1}'
	shell: '{whatshap} phase --max-coverage 15 --include-VCF-phasing --distrust-genotypes --indels --sample NA12878 --reference {input.ref} {input.vcf1} {input.vcf2} {input.bam1} --output {output.vcf} 2> {log}'

rule integrative_whatshap_illumina_SS:
	input: 
		bam1='bam/TRIAL-{trials,[0-9]+}/NA12878.illumina.chr{chromosome,[0-9]+}.cov{icoverage,(all|[0-9]+)}.readgroup.sorted.bam',
		vcf1='vcf/NA12878.benchmark.unphased.chr{chromosome,[0-9]+}.vcf',
		vcf2='StrandPhaseR_TRIAL_{trials,[0-9]+}_{strandseqcoverage,[0-9]+}cells/VCFfiles/chr{chromosome,[0-9]+}_phased.vcf',
		ref=reference,
	output: 
		vcf= 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio0.illumina{icoverage,(all|[0-9]+)}.10x0.chr{chromosome,[0-9]+}.noindels.vcf',
	log: 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio0.illumina{icoverage,(all|[0-9]+)}.10x0.chr{chromosome,[0-9]+}.noindels.vcf.log'
	message: 'Running WHATSHAP on {input.bam1}'
	shell: '{whatshap} phase --max-coverage 15 --include-VCF-phasing --distrust-genotypes --sample NA12878 {input.vcf1} {input.vcf2} {input.bam1} --output {output.vcf} 2> {log}'

rule integrative_whatshap_pacbio_xg:
	input: 
		bam1='bam/TRIAL-{trials,[0-9]+}/NA12878.pacbio.chr{chromosome,[0-9]+}.cov{pcoverage,(all|[0-9]+)}.readgroup.sorted.bam',
		vcf1='vcf/NA12878.benchmark.unphased.chr{chromosome,[0-9]+}.vcf',
		vcf2='vcf/NA12878.10xG.phased.chr{chromosome,[0-9]+}.vcf',
		ref=reference,
	output: 
		vcf= 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x{xcoverage,(all|[0-9]+)}.chr{chromosome,[0-9]+}.noindels.vcf',
	log: 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x{xcoverage,(all|[0-9]+)}.chr{chromosome,[0-9]+}.noindels.vcf.log'
	message: 'Running WHATSHAP on {input.bam1}'
	shell: '{whatshap} phase --max-coverage 15 --include-VCF-phasing --distrust-genotypes --sample NA12878 --reference {input.ref} {input.vcf1} {input.vcf2} {input.bam1} --output {output.vcf} 2> {log}'

rule integrative_whatshap_pacbio_only:
	input: 
		bam1='bam/TRIAL-{trials,[0-9]+}/NA12878.pacbio.chr{chromosome,[0-9]+}.cov{pcoverage,(all|[0-9]+)}.readgroup.sorted.bam',
		vcf1='vcf/NA12878.benchmark.unphased.chr{chromosome,[0-9]+}.vcf',
		ref=reference,
	output: 
		vcf= 'whatshap_pacbio_only/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chr{chromosome,[0-9]+}.noindels.vcf',
	log: 'whatshap_pacbio_only/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chr{chromosome,[0-9]+}.noindels.vcf.log'
	message: 'Running WHATSHAP on {input.bam1}'
	shell: '{whatshap} phase --max-coverage 15 --distrust-genotypes --sample NA12878 --reference {input.ref} {input.vcf1} {input.bam1} --output {output.vcf} 2> {log}'

rule evaluate_whatshap_pacbio_only:
	input:
		truth='vcf/NA12878.benchmark.phased.chr{chromosome,[0-9]+}.vcf',
		phased='whatshap_pacbio_only/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chr{chromosome,[0-9]+}.noindels.vcf',
	output:	'eval/whatshap_pacbio_only/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chr{chromosome,[0-9]+}.noindels.eval',
	log: 'eval/whatshap_pacbio_only/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chr{chromosome,[0-9]+}.noindels.log',
	shell: '{whatshap} compare --names benchmark,whatshap --tsv-pairwise {output} {input.truth} {input.phased} |& tee {log}'

rule integrative_whatshap_SS_xg:
	input: 
		vcf3='vcf/NA12878.10xG.phased.chr{chromosome,[0-9]+}.vcf',
		vcf1='vcf/NA12878.benchmark.unphased.chr{chromosome,[0-9]+}.vcf',
		vcf2='StrandPhaseR_TRIAL_{trials,[0-9]+}_{strandseqcoverage,[0-9]+}cells/VCFfiles/chr{chromosome,[0-9]+}_phased.vcf',
	output: 
		vcf= 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio0.illumina0.10x{xcoverage,(all|[0-9]+)}.chr{chromosome,[0-9]+}.noindels.vcf',
	log: 'whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio0.illumina0.10x{xcoverage,(all|[0-9]+)}.chr{chromosome,[0-9]+}.noindels.vcf.log'
	message: 'Running WHATSHAP on {input.bam1}'
	shell: '{whatshap} phase --max-coverage 15 --include-VCF-phasing --distrust-genotypes --sample NA12878 {input.vcf1} {input.vcf2} {input.vcf3} --output {output.vcf} 2> {log}'

rule evaluate_whatshap:
	input:
		truth='vcf/NA12878.benchmark.phased.chr{chromosome,[0-9]+}.vcf',
		phased='whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina{icoverage,(all|[0-9]+)}.10x{xcoverage,(all|[0-9]+)}.chr{chromosome,[0-9]+}.noindels.vcf',
	output:	'eval/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina{icoverage,(all|[0-9]+)}.10x{xcoverage,(all|[0-9]+)}.chr{chromosome,[0-9]+}.noindels.eval',
	log: 'eval/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina{icoverage,(all|[0-9]+)}.10x{xcoverage,(all|[0-9]+)}.chr{chromosome,[0-9]+}.noindels.log',
	shell: '{whatshap} compare --names benchmark,whatshap --tsv-pairwise {output} --only-snvs {input.truth} {input.phased} |& tee {log}'

rule evaluate_whatshap_indels:
	input:
		truth='vcf/NA12878.benchmark.phased.chr{chromosome,[0-9]+}.vcf',
		phased='whatshap_integrative_phasing/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina0.10x0.chr{chromosome,[0-9]+}.indels.vcf',
	output:	'eval/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina{icoverage,(all|[0-9]+)}.10x{xcoverage,(all|[0-9]+)}.chr{chromosome,[0-9]+}.indels.eval',
	log: 'eval/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio{pcoverage,(all|[0-9]+)}.illumina{icoverage,(all|[0-9]+)}.10x{xcoverage,(all|[0-9]+)}.chr{chromosome,[0-9]+}.indels.log',
	shell: '{whatshap} compare --names benchmark,whatshap --tsv-pairwise {output} {input.truth} {input.phased} |& tee {log}'

rule summary:
	output:'summary.eval',
	input:
		expand('eval/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio{pcoverage}.illumina0.10x0.chr{chromosome}.noindels.eval', chromosome=chromosome, strandseqcoverage=strandseqcoverage, pcoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio{pcoverage}.illumina0.10x0.chr{chromosome}.indels.eval', chromosome=chromosome, strandseqcoverage=strandseqcoverage, pcoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio0.illumina{icoverage}.10x0.chr{chromosome}.noindels.eval',chromosome=chromosome, strandseqcoverage=strandseqcoverage, icoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio0.illumina0.10x{xcoverage}.chr{chromosome}.noindels.eval',chromosome=chromosome, strandseqcoverage=strandseqcoverage, xcoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells0.pacbio{pcoverage}.illumina0.10x{xcoverage}.chr{chromosome}.noindels.eval',chromosome=chromosome, xcoverage=coverage, pcoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells0.pacbio0.illumina{icoverage}.10x0.chr{chromosome}.noindels.eval', chromosome=chromosome, icoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio{pcoverage}.illumina0.10x0.chr{chromosome}.indels.eval', chromosome=chromosome, strandseqcoverage=strandseqcoverage, pcoverage=coverage, trials=trials),
		expand('eval/whatshap_pacbio_only/TRIAL-{trials}/strandseqcells0.pacbio{pcoverage}.illumina0.10x0.chr{chromosome}.noindels.eval', chromosome=chromosome, pcoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells0.pacbio{pcoverage}.illumina0.10x0.chr{chromosome}.indels.eval', chromosome=chromosome, pcoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells0.pacbio0.illumina{icoverage}.10x0.chr{chromosome}.noindels.eval',chromosome=chromosome, icoverage=coverage, trials=trials),
	message: 'Aggregating statistics to {output}'
	run:
		first = input[0]
		rest = ' '.join(input[1:]) 
		shell('(cat {first} && ( cat {rest} | grep -v \'^#\' ) ) > {output}')
