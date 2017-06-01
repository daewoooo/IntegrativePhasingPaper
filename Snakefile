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
coverage = [2,3,4,5,10,15,25,30,'all']
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

rule evaluate_whatshap_10xg_only:
	input:
		truth='vcf/NA12878.benchmark.phased.chr{chromosome,[0-9]+}.vcf',
		phased='vcf/NA12878.10xG.phased.chr{chromosome,[0-9]+}.vcf',
	output:	'eval/whatshap_10xG_only/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio0.illumina0.10xall.chr{chromosome,[0-9]+}.noindels.eval',
	log: 'eval/whatshap_10xG_only/TRIAL-{trials,[0-9]+}/strandseqcells0.pacbio0.illumina0.10xall.chr{chromosome,[0-9]+}.noindels.log',
	shell: '{whatshap_cmp} compare --names benchmark,whatshap --only-snvs --tsv-pairwise {output} {input.truth} {input.phased} |& tee {log}'


rule evaluate_whatshap_SS_only:
	input:
		truth='vcf/NA12878.benchmark.phased.chr{chromosome,[0-9]+}.vcf',
		phased='StrandPhaseR_TRIAL_{trials,[0-9]+}_{strandseqcoverage,[0-9]+}cells/VCFfiles/chr{chromosome,[0-9]+}_phased.vcf',
	output:	'eval/whatshap_SS_only/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio0.illumina0.10x0.chr{chromosome,[0-9]+}.noindels.eval',
	log: 'eval/whatshap_SS_only/TRIAL-{trials,[0-9]+}/strandseqcells{strandseqcoverage,[0-9]+}.pacbio0.illumina0.10x0.chr{chromosome,[0-9]+}.noindels.log',
	shell: '{whatshap_cmp} compare --names benchmark,whatshap --only-snvs --tsv-pairwise {output} {input.truth} {input.phased} |& tee {log}'

# Added rules from Tobias Snakefile
rule download_bam:
	output: 'bam/{sample}.bam'
	log: 'bam/{sample}.wgetlog'
	resources: download=1
	shell: 'wget --output-file={log} -O {output} ftp://ftp.sra.ebi.ac.uk/vol1/ERA172/ERA172924/bam/{wildcards.sample}_S1.bam'
		

rule download_giab_vcf:
	input: 'names/NA12878.txt'
	output: 'vcf/giab/NA12878.vcf.gz'
	log: 'vcf/giab/NA12878.vcf.gz.wgetlog'
	resources: download=1
	shell: 'wget --output-file={log} -O - ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz | bcftools reheader -s {input} > {output}'

rule download_platinum_vcf:
	output: 'vcf/platinum/NA12878.vcf.gz'
	log: 'vcf/platinum/NA12878.vcf.gz.wgetlog'
	resources: download=1
	shell: 'wget --output-file={log} --password "" -O {output} ftp://platgene_ro@ussd-ftp.illumina.com/2016-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz'

rule vcf_sitesonly:
	input: 
		vcf='vcf/{what}/{sample}.vcf.gz',
		tbi='vcf/{what}/{sample}.vcf.gz.tbi',
	output:
		vcf='vcf-sitesonly/{what}/{sample}.{chromosome}.vcf.gz'
	log: 'vcf-sitesonly/{what}/{sample}.{chromosome}.log'
	shell: 'bcftools view -r {wildcards.chromosome} -G {input.vcf} | bgzip > {output.vcf}'

rule sample_name:
	output: 'names/{sample}.txt'
	shell: 'echo {wildcards.sample} > {output}'

rule compare_giab_platinum:
	input:
		platinum='vcf/platinum/{sample}.vcf.gz',
		giab='vcf/giab/{sample}.vcf.gz',
	output: 'comparison/giab-platinum.{sample}.txt'
	log: 'comparison/giab-platinum.{sample}.log'
	shell: './compare-vcfs.py --names GIAB,PLATINUM {wildcards.sample} {input.giab} {input.platinum} > {output} 2> {log}'

rule compare_platinum_freebayes:
	input:
		platinum='vcf/platinum/{sample}.vcf.gz',
		freebayes='freebayes-retype-merged/all-filtered.vcf.gz',
	output: 'comparison/platinum-freebayes.{sample}.txt'
	log: 'comparison/platinum-freebayes.{sample}.log'
	shell: './compare-vcfs.py --names PLATINUM,FREEBAYES {wildcards.sample} {input.platinum} {input.freebayes} > {output} 2> {log}'

rule index_bam:
	input: 
		bam='{file}.bam'
	output:
		bai='{file}.bam.bai'
	shell:
		'samtools index {input.bam}'

rule bgzip:
	input: '{file}.vcf'
	output: '{file}.vcf.gz'
	shell: 'bgzip -c {input} > {output}'

rule tabix:
	input: '{file}.vcf.gz'
	output: '{file}.vcf.gz.tbi'
	shell: 'tabix {input}'
			
rule freebayes_genotype:
	input:
		vcf='vcf-sitesonly/platinum/NA12878.{chromosome}.vcf.gz',
		vcfidx='vcf-sitesonly/platinum/NA12878.{chromosome}.vcf.gz.tbi',
		bams=expand('bam/{sample}.bam', sample=samples),
		bais=expand('bam/{sample}.bam.bai', sample=samples),
	output:
		vcf='freebayes-retype/{chromosome}.vcf.gz'
	log: 'freebayes-retype/{chromosome}.log'
	shell:
		'(time samtools merge -R {wildcards.chromosome} - {input.bams} |  freebayes -f {ref} --haplotype-basis-alleles {input.vcf} -@ {input.vcf} --stdin | bgzip > {output.vcf} ) 2> {log}'

rule merge_freebayes:
	input:
		vcfs=expand('freebayes-retype/{chromosome}.vcf.gz', chromosome=chromosomes),
		vcfidx=expand('freebayes-retype/{chromosome}.vcf.gz.tbi', chromosome=chromosomes),
	output:
		'freebayes-retype-merged/all.vcf.gz'
	log:
		'freebayes-retype-merged/all.vcf.gz.log'
	shell:
		'(bcftools concat {input.vcfs} | vcf-sort | bgzip > {output}) 2> {log}'

rule filter_freebayes_vcf:
	input:
		vcf='freebayes-retype-merged/all.vcf.gz'
	output:
		vcf='freebayes-retype-merged/all-filtered.vcf.gz'
	shell:
		'zcat {input.vcf} | awk \'($0 ~ /^#/) || (($6 != ".") && ($6 >= 30))\' | bgzip > {output.vcf}'

rule create_ped:
	output: 'ped/family.ped'
	shell: 'echo "FAM1 NA12878 NA12891 NA12892 0 1" > {output}'
rule whatshap_trio_pure_genetic:
	input:
		vcf='freebayes-retype-merged/all-filtered.vcf.gz',
		ped='ped/family.ped',
	output:
		vcf='whatshap/trio-pure-genetic-phasing.vcf.gz',
	log: 'whatshap/trio-pure-genetic-phasing.log'
	shell: '({whatshap} phase --chromosome chr1 --ped {input.ped} {input.vcf} | bgzip > {output.vcf}) 2> {log}'

# Some hard-coded wildcards to restrict the analysis.
rule whatshap_trio_pure_genetic:
	input:
		vcf='freebayes-retype-merged/all-filtered.vcf.gz',
		ped='ped/family.ped',
	output:
		vcf='whatshap/trio-pure-genetic-phasing.vcf.gz',
	log: 'whatshap/trio-pure-genetic-phasing.log'
	shell: '({whatshap} phase --chromosome chr1 --ped {input.ped} {input.vcf} | bgzip > {output.vcf}) 2> {log}'

rule whatshap_trio:
	input:
		vcf='freebayes-retype-merged/all-filtered.vcf.gz',
		ped='ped/family.ped',
		bam = 'bam/TRIAL-3/NA12878.pacbio.chr1.covall.readgroup.sorted.bam'
	output:
		vcf='whatshap/trio-phasing.vcf.gz',
	log: 'whatshap/trio-phasing.log'
	shell: '({whatshap} phase --chromosome chr1 --max-coverage 45 --ped {input.ped} {input.vcf} {input.bam}| bgzip > {output.vcf}) 2> {log}'


rule extract_child:
	input: 'whatshap/trio-phasing.vcf.gz', 'whatshap/trio-pure-genetic-phasing.vcf.gz'
	output: 'whatshap/trio-phasing.child.vcf', 'whatshap/trio-pure-genetic-phasing.child.vcf'
	shell: 'vcf-subset -c NA12878 {input[0]} > {output[0]} && vcf-subset -c NA12878 {input[1]} > {output[1]}'

rule eval_whatshap_trio:
	input:
		truth='vcf/NA12878.benchmark.phased.chr1.vcf',
		phased='whatshap/trio-phasing.child.vcf',
	output:	'eval/whatshap/trio-phasing.child.eval',
	log: 'eval/whatshap/trio-phasing.child.log',
	shell: '{whatshap} compare --names benchmark,whatshap --only-snvs --tsv-pairwise {output} {input.truth} {input.phased} |& tee {log}'

rule eval_whatshap_trio_pure_genetic:
	input:
		truth='vcf/NA12878.benchmark.phased.chr1.vcf',
		phased='whatshap/trio-pure-genetic-phasing.child.vcf',
	output:	'eval/whatshap/whatshap/trio-pure-genetic-phasing.child.eval',
	log: 'eval/whatshap/whatshap/trio-pure-genetic-phasing.child.log',
	shell: '{whatshap} compare --names benchmark,whatshap --only-snvs --tsv-pairwise {output} {input.truth} {input.phased} |& tee {log}'
	

rule summary:
	output:'summary.eval',
	input:
		expand('eval/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio{pcoverage}.illumina0.10x0.chr{chromosome}.noindels.eval', chromosome=chromosome, strandseqcoverage=strandseqcoverage, pcoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio{pcoverage}.illumina0.10x0.chr{chromosome}.indels.eval', chromosome=chromosome, strandseqcoverage=strandseqcoverage, pcoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio0.illumina{icoverage}.10x0.chr{chromosome}.noindels.eval',chromosome=chromosome, strandseqcoverage=strandseqcoverage, icoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio0.illumina0.10x{xcoverage}.chr{chromosome}.noindels.eval',chromosome=chromosome, strandseqcoverage=strandseqcoverage, xcoverage=['all'], trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells0.pacbio{pcoverage}.illumina0.10x{xcoverage}.chr{chromosome}.noindels.eval',chromosome=chromosome, xcoverage=['all'], pcoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells0.pacbio0.illumina{icoverage}.10x0.chr{chromosome}.noindels.eval', chromosome=chromosome, icoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio{pcoverage}.illumina0.10x0.chr{chromosome}.indels.eval', chromosome=chromosome, strandseqcoverage=strandseqcoverage, pcoverage=coverage, trials=trials),
		expand('eval/whatshap_pacbio_only/TRIAL-{trials}/strandseqcells0.pacbio{pcoverage}.illumina0.10x0.chr{chromosome}.noindels.eval', chromosome=chromosome, pcoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells0.pacbio{pcoverage}.illumina0.10x0.chr{chromosome}.indels.eval', chromosome=chromosome, pcoverage=coverage, trials=trials),
		expand('eval/TRIAL-{trials}/strandseqcells0.pacbio0.illumina{icoverage}.10x0.chr{chromosome}.noindels.eval',chromosome=chromosome, icoverage=coverage, trials=trials),
		expand('eval/whatshap_10xG_only/TRIAL-{trials}/strandseqcells0.pacbio0.illumina0.10xall.chr{chromosome}.noindels.eval', chromosome=chromosome, trials=trials),
		expand('eval/whatshap_SS_only/TRIAL-{trials}/strandseqcells{strandseqcoverage}.pacbio0.illumina0.10x0.chr{chromosome}.noindels.eval', chromosome=chromosome, trials=trials, strandseqcoverage=strandseqcoverage),
		expand('eval/whatshap/whatshap/trio-pure-genetic-phasing.child.eval'),
		expand('eval/whatshap/trio-phasing.child.eval')
	message: 'Aggregating statistics to {output}'
	run:
		first = input[0]
		rest = ' '.join(input[1:]) 
		shell('(cat {first} && ( cat {rest} | grep -v \'^#\' ) ) > {output}')
