#!/usr/bin/env python3
import sys
import subprocess
from collections import namedtuple, defaultdict
from bisect import bisect_left
from whatshap.args import HelpfulArgumentParser as ArgumentParser

Variant = namedtuple('Variant', ['chromosome', 'position', 'ref', 'alt'])

valid_genotypes = set(['.', '0/0', '0/1', '1/1', '0|1', '1|0', '1|1', '0|0'])
unphase_genotype = {
	'.':'.', 
	'0/0':'0/0', 
	'0/1':'0/1',
	'1/1':'1/1',
	'1|0':'0/1',
	'0|1':'0/1',
	'1|1':'1/1',
	'0|0':'0/0'
}


def read_variants(filename, sample, only_genotype=None, snvs=True, indels=False):
	'''Returns a map (chr, pos, ref, alt) --> genotype'''
	command = "bcftools query -s {}  -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n' {}".format(sample, filename)
	result = {}
	skipped_invalid = 0
	skipped_multiallelic = 0
	skipped_n = 0
	skipped_snps = 0
	skipped_nonsnps = 0
	bcftools = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
	for line in (s.decode('ascii') for s in bcftools.stdout.splitlines()):
	#for line in subprocess.getoutput(command).split('\n'):
		if line.startswith('Warning:') or line.startswith('[W::'):
			continue
		fields = line.split('\t')
		#if len(fields) != 5:
			#print(fields)
			#continue
		assert len(fields) == 5, 'line="{}", fields={}'.format(line, fields)
		alts = fields[3].split(',')
		if len(alts) > 1:
			skipped_multiallelic += 1
			continue
		variant = Variant(fields[0], int(fields[1]), fields[2], alts[0])
		genotype = fields[4]
		if ('N' in variant.ref) or ('N' in variant.alt):
			skipped_n += 1
			continue
		is_snv = (len(variant.ref) == 1) and (len(variant.alt) == 1)
		if is_snv and not snvs:
			skipped_snps += 1
			continue
		if (not is_snv) and (not indels):
			skipped_nonsnps += 1
			continue
		if genotype not in valid_genotypes:
			skipped_invalid += 1
			continue
		genotype = unphase_genotype[genotype]
		if only_genotype is not None:
			if genotype != only_genotype:
				continue
		result[variant] = genotype
	if skipped_invalid > 0:
		print('Skipped {} invalid genotypes when reading {}'.format(skipped_invalid, filename))
	if skipped_multiallelic > 0:
		print('Skipped {} variants that are multi-allelic {}'.format(skipped_multiallelic, filename))
	if skipped_n > 0:
		print('Skipped {} variants because REF or ALT contains an N character when reading {}'.format(skipped_n, filename))
	if skipped_snps > 0:
		print('Skipped {} variants that are SNVs when reading {}'.format(skipped_snps, filename))
	if skipped_nonsnps > 0:
		print('Skipped {} variants that are not SNVs when reading {}'.format(skipped_nonsnps, filename))
	return result

def read_annotation_track(filename):
	result = defaultdict(list)
	length = 0
	n = 0
	for line in open(filename):
		fields = line.split()
		assert len(fields) == 3
		chromosome, start, end = fields[0], int(fields[1]), int(fields[2]) + 1
		assert start < end
		if len(result[chromosome]) > 0:
			assert result[chromosome][-1] <= start 
		result[chromosome] += [start, end]
		n += 1
		length += end-start
	print('Read {} intervals with a total length of {} from {}'.format(n, length, filename))
	return result

def in_annotated_region(annotation, chromosome, position):
	try:
		l = annotation[chromosome]
	except KeyError:
		return False
	i = bisect_left(l, position)
	if i >= len(l):
		return False
	if l[i] != position:
		result = (i % 2) == 1
	else:
		result = (i % 2) == 0
	#print('in_annotated_region({},{})={}, i={}'.format(chromosome, position, result, i))
	return result

def main():
	parser = ArgumentParser(prog='compare-vcfs.py', description=__doc__)
	parser.add_argument('--only-genotype', default=None,
		help='Restrict analysis to given genotype from {"0/0", "0/1", "1/1", "."}')
	parser.add_argument('--names', default=None,
		help='Comma-separated list of data set names')
	parser.add_argument('--annotation', default=None,
		help='BED track with annotation.')
	parser.add_argument('--also-indels', default=False, action='store_true',
		help='Also work on indels.')
	parser.add_argument('--only-indels', default=False, action='store_true',
		help='Only work on indels.')
	parser.add_argument('sample', metavar='SAMPLE', help='sample to work on')
	parser.add_argument('vcf', nargs='+', metavar='VCF', help='VCF files to compare')
	args = parser.parse_args()

	if args.names is None:
		names = args.vcf
	else:
		names = args.names.split(',')
		assert len(names) == len(args.vcf)

	snps = not args.only_indels
	indels = args.only_indels or args.also_indels

	annotation = None
	if args.annotation is not None:
		annotation = read_annotation_track(args.annotation)
	#print(annotation)

	variant_tables = [read_variants(filename, args.sample, args.only_genotype, snps, indels) for filename in args.vcf]
	variant_sets = [set(vt.keys()) for vt in variant_tables]
	union = set()
	for s in variant_sets:
		union.update(s)
	intersection = set(union)
	for s in variant_sets:
		intersection.intersection_update(s)
	print('---------------------------------------------------------')
	print('Variants union: ', len(union))
	print('Variants intersection: ', len(intersection))
	if annotation is not None:
		annotated_variants = {v:'1/1' for v in union if in_annotated_region(annotation, v.chromosome, v.position)}
		variant_tables.append(annotated_variants)
		names.append('ANNOTATED')
	longest_name = max(len(n) for n in names)

	site_counts = defaultdict(int)
	for v in union:
		present = tuple((v in vt) for vt in variant_tables)
		site_counts[present] += 1
	
	print('---------------------------------------------------------')
	print('Site congruence (regardless of genotypes):')
	for present in sorted(site_counts.keys()):
		count = site_counts[present]
		present_names = [n for n,p in zip(names,present) if p]
		print('{%s}: %d' % (','.join(present_names), count))
		
	present_gts = frozenset(['0/1','1/1'])
	variant_presence_counts = defaultdict(int)
	for v in union:
		present = tuple((v in vt) and (vt[v] in present_gts) for vt in variant_tables)
		variant_presence_counts[present] += 1
	print('---------------------------------------------------------')
	print('Variants present (i.e. site in VCF and typed 0/1 or 1/1):')
	for present in sorted(variant_presence_counts.keys()):
		count = variant_presence_counts[present]
		present_names = [n for n,p in zip(names,present) if p]
		print('{%s}: %d' % (','.join(present_names), count))
	
	print('---------------------------------------------------------')
	print('Genotype congruence table: ')
	if annotation is not None:
		variant_tables = variant_tables[:-1]
		names = names[:-1]
	genotype_congruence_counts = defaultdict(int)
	for v in intersection:
		gts = tuple(vt[v] for vt in variant_tables)
		genotype_congruence_counts[gts] += 1
	gt_list = list(genotype_congruence_counts.keys())
	gt_list.sort()
	print(' '.join(n.rjust(longest_name) for n in names))
	for gts in gt_list:
		print(' '.join(gt.rjust(longest_name) for gt in gts) + '{:10d}'.format(genotype_congruence_counts[gts]))
	
	congruent = 0
	total = 0
	for gts, count in genotype_congruence_counts.items():
		gt_set = set(gts)
		if '.' in gt_set:
			continue
		total += count
		if len(gt_set) == 1:
			congruent += count
	print('Genotype congruence for sites successfully typed in all VCFs: {} / {} = {}%'.format(congruent, total, str(congruent/total*100) if total>0 else 'nan'))



if __name__ == '__main__':
	main()
