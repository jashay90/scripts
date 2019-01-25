#!/usr/bin/perl
# Julie Shay
# October 31: The original midas snvs output format has changed with the new version of MIDAS, 
# so I am updated this script accordingly.
# March 20: This version will take the original midas snvs output, rather than the merge output.
# March 1, 2017
# This script will take a Mauve SNPs file and a snps_alt_alleles file from merge_midas output
# and for each sample in the snps_alt_alleles file, it will count the number of matching and non-matching
# SNV calls.
use strict;
use warnings;
use LWP::Simple;
use Getopt::Long;

my $mauvein = "mauve_salmonella_snvs_sorted";
my $midasin = "snps/output/Salmonella_enterica_58156.snps.gz";
my $mincoverage = 10;
my $notlcbin = "";	# optional to check that a site is in an LCB before labelling it as a false positive
my $help = 0;

GetOptions('mauve=s' => \$mauvein, 'midas=s' => \$midasin, 'mincoverage=i' => \$mincoverage, 'notlcb=s' => \$notlcbin, 'help' => \$help);
if ($help) {
	print "Options\:\n\-mauve\t\: file with list of true positive SNVs, sorted by position in reference sequence and with no header\n\-midas\t\: compressed MIDAS SNP output for a single species\n\-mincoverage\t\: minimum coverage to call a SNV. Default 10\n\-notlcb\t\: Optional file with regions of the reference sequence that were not aligned by Mauve.\n";
	exit;
}
my $maxref_freq = 0.5;	# not really sure what this is supposed to be.
my %nuclorder=("A", 4, "C", 5, "G", 6, "T", 7);
my $tp = ();
my $fp = ();
my $fn = ();
my $fnbcmissing = ();
my @mauvecol = ();
my @midascol = ();
my @notlcbstart = ();
my @notlcbend = ();
my $alt = 0;
if ($notlcbin){
	open(NONLCB, $notlcbin);
	while ($coords = <NONLCB>){
		chomp($coords);
		push(@notlcbstart, $coords);
		push(@notlcbend, $coords);
		$notlcbstart[$#notlcbstart] =~ s/\t.+//;
		$notlcbend[$#notlcbend] =~ s/.+\t//;
	}
	close(NONLCB);
}
open(MAUVE, $mauvein);
open(MIDAS, "gunzip -c $midasin |");
chomp($name[0]);
my @name = split(/\t/, $name[0]);
# <MAUVE>; don't need this if Mauve file is missing its header
while ($pos = <MAUVE>){
	if (!($pos =~ /^.N/)){
		chomp($pos);
		@mauvecol = split(/\t/, $pos);
		$mauvecol[1] = $mauvecol[0];
		$mauvecol[1] =~ s/^.//;
		$mauvecol[0] =~ s/.$//;
		if ($midascol[1] < $mauvecol[5]){
			while ($line = <MIDAS>){
				chomp($line);
				@midascol = split(/\t/, $line);
				if ($midascol[1] >= $mauvecol[5]){
					last;
				} elsif (($midascol[3] >= $mincoverage) && ($midascol[2] ne "N"))  {
					$alt = 0;
					for ($x = 4; $x <= 7; $x++){
						if (($x != $nuclorder{$midascol[2]}) && ($midascol[$x] >= ($midascol[3] / 2))){
							$alt = 1;
							last;
						}
					}
					if ($alt){
						if (!($notlcbin)){
							$fp += 1;
						} else {
							for ($x = 0; $x <= $#notlcbend; $x++){
								if ($notlcbend[$x] >= $midascol[1]){
									if ($notlcbstart[$x] > $midascol[1]){
										$fp += 1;
									}
									last;
								}
							}
						}
					}
				}
			}
		}
		# check whether $line contains the same SNV as $pos
		if ($midascol[3] >= $mincoverage){
				if ($midascol[2] ne $mauvecol[1]){
					if ($midascol[2] eq getcom($mauvecol[1])){
						$mauvecol[0] = getcom($mauvecol[0]);
					} else {
						print "Error at $pos";
					}
				}
				#test that the same alt allele is called here
				# if the same alt allele IS present, it's probably a tp
				if ($midascol[$nuclorder{$mauvecol[0]}] >= ($midascol[3] / 2)){
					$tp += 1
				} elsif ($midascol[$nuclorder{$midascol[2]}] >= ($midascol[3] / 2)){
					$fn += 1;
				} else {
					$alt = 0;
					for ($x = 4; $x <= 7; $x++){
						if ($midascol[$x] >= ($midascol[3] / 2)){
							$fp += 1;
							$alt = 1;
							last;
						} 
					}
					if (!($alt)){
						$fn += 1;
					}
				}
		} else {
			$fn += 1;
			$fnbcmissing += 1;
		}
	}
}
close(MAUVE);
close(MIDAS);

print "$midasin\ntrue positives\: $tp\nfalse positives\: $fp\nfalse negatives\: $fn\nfalse negatives because the site isnt sufficiently covered in the sample\: $fnbcmissing\n";

print "\n";

sub getcom {
	if ($_[0] eq "A"){
		return "T";
	} elsif ($_[0] eq "T"){
		return "A";
	} elsif ($_[0] eq "C"){
		return "G";
	} elsif ($_[0] eq "G"){
		return "C";
	} elsif ($_[0] eq "K"){
		return "M";
	} elsif ($_[0] eq "M"){
		return "K";
	} elsif ($_[0] eq "R"){
		return "Y";
	} elsif ($_[0] eq "Y"){
		return "R";
	} elsif ($_[0] eq "S"){
		return "S";
	} elsif ($_[0] eq "W"){
		return "W";
	} else {
		return "N";
	}
}
