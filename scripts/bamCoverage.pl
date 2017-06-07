#!/usr/bin/perl -W
use strict;

open OUT, '>', "coverage.txt" or die;

my @bam_files = ();
foreach (@ARGV) { push @bam_files, $_; }

#my $len = 2_864_785_223; #human_GRCh37 (hg19)
my $len = 225_280_621; #length for chromosome 1 (hg19) (without N's)

print OUT "Filename\tCoverage\tRange\n";

my $num = 0;
my $bases = 0;
foreach my $bam (@bam_files) { 
 
	open IN, "samtools view -hb -F256 $bam | samtools mpileup - |"; #flag 256 filters out secondary alignments

	while (<IN>) {
    		my @a = split /\t/;
    		$bases += $a[3];
    		$num++ if $a[3];
	}

	print "$bam\n";
	printf  "Mean coverage  : %.2f\n", ($num/$len)*100;
	#printf  "Depth of coverage :  %.4f\n", $bases/2_864_785_223;
	printf  "Depth of coverage :  %.4f\n", $bases/$len;

	print OUT"$bam\t";
	printf OUT "%.2f\t", ($num/$len)*100;
	printf OUT "%.4f\n", $bases/$len;

	($num, $bases)  = (0,0);

	close IN;
}
