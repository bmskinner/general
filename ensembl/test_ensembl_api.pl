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
# Fetch data
#######################

# Adaptors connect to a specific database
# Here we connect to raw DNA and the transcript databases
# Note we don't actually use the transcript adaptor in this script;
# it's just here to show the syntax to connect
my $slice_adaptor = $registry->get_adaptor( 'Pig', 'Core', 'Slice' );
my $tr_adaptor    = $registry->get_adaptor( 'Pig', 'Core', 'Transcript' );

print "Got slice and transcript adaptors\n";

# Example using a fixed chromosome and coordinates (a 'slice' of the genome)

my $chr = "X";
my $start = 58063021;
my $stop = 58074244;

# Get the specified slice
my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $start, $stop );
my $sequence = $slice->seq(); # get the DNA sequence from the slice

print $sequence . "\n";

print "Done\n";




