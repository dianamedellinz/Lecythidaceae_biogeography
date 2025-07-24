#!/usr/bin/perl

use warnings;
use strict;

my %genelengths; # hash of hashes; 1st level keys are sample names, 2nd level keys are gene names, values are lengths in bp
my @genelist;
my %sampleNames; # keys are sample names, values are 1

my $geneListFile = shift @ARGV;
open(my $IN,"<",$geneListFile) or die "Couldn't open $geneListFile for reading.\n";

while (my $gene = <$IN>) {
	chomp $gene;
	push @genelist, $gene;
}

foreach my $g (@genelist) {
	my $fastaFile = "combined_a353." . $g . ".lecy.fa";
	open (my $FIN,"<",$fastaFile) or die "Couldn't open $fastaFile for reading.\n";
	my $sample;
	while (my $line = <$FIN>) {
		chomp $line;
		if ($line =~/^>/) {
			$sample = $line;
			$sampleNames{$sample} = 1;
		}
		else {
			$genelengths{$sample}{$g} = length($line);
		}
	}
	close $FIN;
}

print "Sample";
foreach my $gene (@genelist) {
	print "\t$gene";
}
print "\n";

foreach my $s (sort keys %sampleNames) {
	print $s;
	foreach my $g (@genelist) {
		if (exists $genelengths{$s}{$g}) {
			print "\t", $genelengths{$s}{$g};
		}
		else {
			print "\t0";
		}
	}
	print "\n";
}