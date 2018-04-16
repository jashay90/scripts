#!/usr/bin/perl
# Julie Shay
# December 7, 2017
# I made this by modifying ~/insilico/taxons.pl
# April 11, 2018: Cleaning this script up for Lee

use LWP::Simple;
use Getopt::Long;
use strict;
use warnings;

my $obsin = "species_profile.txt";
my $obsout = "kronainput.tab";
my $kronaout = "kronaout.html";
my $help = 0;
GetOptions('in=s' => \$obsin, 'tab=s' => \$obsout, => 'krona=s' => \$kronaout, 'h' => \$help, 'help' => \$help);

if ($help){
	print "Options\:\n\-in\t\: input file \(species_profile file from MIDAS\)\n\-tab\t\: tab-delimited output\n\-krona\t\: krona\'s html output\n";
	exit;
}

my @col;
my %observed;
my $sum = 0;
# grab the proportions and species names from the input file
open(OIN, "grep -vP \"\t0.0\t\" $obsin |");
<OIN>;
while (<OIN>){
	chomp($_);
	@col = split(/\t/, $_);
	$sum += $col[1];
	$col[2] = $col[0];
	$col[0] =~ s/\_.*//g;
	$col[2] =~ s/$col[0]_//;
	$col[2] =~ s/\_.*//;
	$col[3] = &gettaxon($col[0], $col[2]); # get the rest of the taxonomy from NCBI
	# If there's more than one strain of the same species in the input file,
	# just sum the proportions
	$observed{$col[3]} += $col[1];
}
close(OIN);

# print the proportions and taxonomies
my $shortened = "";
open(OBOUT, ">$obsout");
foreach (keys %observed){
	$shortened = $_;
	$shortened =~ s/\;\s/\t/g;
	print OBOUT "$observed{$_}\t$shortened\n";
}

#create the actual krona html output
`module load krona`;
`ktImportText -o $kronaout $obsout`;

# given a genus and species name, get the higher levels of the taxonomy.
sub gettaxon {
	my $ge = $_[0];
	my $sp = $_[1];
	my $base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";
	my $url = "$base/esearch.fcgi?db=taxonomy&term=";
	if ($sp eq "sp"){
		$url .= "genus%5BRank%5D%20AND%20%22$ge%22%5BOrganism%5D";
	} else {
		$url .= "species%5BRank%5D%20AND%20%22$ge+$sp%22%5BOrganism%5D";
	}
	my $result = get($url);
	if (($result =~ /No items found/) || (!($result =~ /\<Id\>/))){
		if ($sp eq "sp"){
			die "Cannot find this taxonomy.\n";
		} else {
			return &gettaxon($ge, "sp");
		}
	} else {
		$result =~ s/\n//g;
		$result =~ s/.+\<Id\>//g;
		$result =~ s/\<\/Id\>.+//g;
		$url = "$base/efetch.fcgi?db=taxonomy&id=$result";
		$result = get($url);
		$result =~ s/\n//g;
		$result =~ s/.+\<Lineage\>//g;
		$result =~ s/\<\/Lineage\>.+//g;
		if ($sp eq "sp"){
			$result .= "; $ge; unspecified";
		} else {
			$result .= "; $ge $sp";
		}
		$result =~ s/\;[^\;]*group[^\;]*//g;
		return $result;
	}
}
