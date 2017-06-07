#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $usage = "Usage $0 <bam_1, bam_2...>\n";
die $usage if @ARGV < 1;

open OUT, '>', "read_length.txt" or die "Can't create a file!!!\n";

my @files = ();

foreach (@ARGV) { push @files, $_; }

my @lengths = ();
my $sum = 0;
my $count = 0;
foreach my $bam (<@files>) {
	warn "Working on bamfile $bam\n";

	open (BAM, "samtools view $bam|") or die "Data processing failed";

	while(<BAM>) {
		chomp;
		my ($seq, $flag) = (split "\t", $_)[9,1];
		next if $flag & 256; #skip secondary alignments
		my $len = length($seq);
		$sum += $len;
		$count++;
	}
}

my $mean = $sum/$count;
print "\nMEAN read length: $mean\n";
print OUT "MEAN read length: $mean\n";

