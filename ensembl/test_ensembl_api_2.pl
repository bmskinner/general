#!/usr/bin/perl

# Where to find the Ensembl modules
use lib "/usr/local/ensembl/ensembl/BioPerl-1.6.924";
use lib "/usr/local/ensembl/ensembl/modules";
use lib "/usr/local/ensembl/ensembl-compara/modules";
use lib "/usr/local/ensembl/ensembl-funcgen/modules";
use lib "/usr/local/ensembl/ensembl-io/modules";
use lib "/usr/local/ensembl/ensembl-metadata/modules";
use lib "/usr/local/ensembl/ensembl-tools/modules";
use lib "/usr/local/ensembl/ensembl-variation/modules";

use strict; # strict language rules enforced
use warnings; # warn about poor code style

# Import relevant modules
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Slice;

#######################
# Connect to Ensembl
#######################

print "Loading registry (may take a while)\n";
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

print "Registry loaded\n";

#######################
# Setup
#######################
my $inFile = $ARGV[0]; # input file - first argument

#######################
# Fetch data
#######################

# Adaptors connect to a specific database
# Here we connect to raw DNA and the transcript databases
my $slice_adaptor = $registry->get_adaptor( 'Pig', 'Core', 'Slice' );
my $tr_adaptor    = $registry->get_adaptor( 'Pig', 'Core', 'Transcript' );

print "Got slice and transcript adaptors\n";

# Open a file handle to read the input file accessed by 'IN'
open(IN, "<", $inFile) or die "Can't open $inFile : $!"; 

print "Read input file $inFile\n";

while(<IN>) { # for every line in the file
	chomp; # remove trailing whitespace

	# assign new variables with each value in the tsv 
	my ($chr, $start, $stop) = split('\t',$_); 
	# Note these are locally scoped to the loop

	print $chr . "\t" . $start . "\t" . $stop . "\n";

	# Get the specified slice
	my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $start, $stop );
	my $sequence = $slice->seq(); # get the DNA sequence from the slice

	print $sequence . "\n";
}
close(IN); # close the file

print "Done\n";
