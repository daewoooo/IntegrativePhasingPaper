#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $usage = "Usage $0 <vcf_1, vcf_2...>\n";
die $usage if @ARGV < 1;

open OUT, ">", "indel_counts.txt" or die "Can't create file";
print OUT "filename\tindelCount\tSS_cells\tPB_cov\tIL_cov\ttenX_cov\n";

my @files = ();
foreach (@ARGV) { push @files, $_; }

foreach my $vcf (@files) {
	#print "$vcf\n";
	
	$vcf =~ /cells(\d+).pacbio(\d+|\w+).illumina(\d+|\w+).10x(\d+|\w+)/;

	my $SS = $1;
	my $PB = $2;
	my $IL = $3;
	my $tenX = $4;	

	open IN , "<", $vcf or die "Can't read file";

	my $indel_counter = 0;
	while(<IN>) {
		next if $_ =~ /\#+/;

		my ($chr, $pos, $rsID, $ref, $alt, $qual, $filter, $info, $format, $gen) = (split "\t", $_);

		my $del = 1 if length($ref) > length($alt);		
		my $ins = 1 if length($ref) < length($alt);

		next if $gen !~ /\|/;
		next unless $del or $ins;
		
		$indel_counter++;
		($del, $ins) = (0,0);
		#print "$_\n";
	}
	print OUT "$vcf\t$indel_counter\t$SS\t$PB\t$IL\t$tenX\n";
}
