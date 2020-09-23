#!/usr/bin/perl
use Bio::DB::Taxonomy;

my $idx_dir = $ARGV[0];
# taxon files can be downloaded from: ftp://ftp.ncbi.nih.gov/pub/taxonomy/
my $db = Bio::DB::Taxonomy->new(-source    => 'flatfile',
                                -nodesfile => $idx_dir . '/nodes.dmp',
                                -namesfile => $idx_dir . '/names.dmp',
                                -directory => $idx_dir);
# test
my $taxonid = 9606;
my $taxon = $db->get_taxon(-taxonid => $taxonid);
my $name = $taxon->scientific_name;

if ($name eq "Homo sapiens") {
  print "Index files for taxonomy database were successfully generated!\n";
} else {
  print "Something wrong happened while indexing taxonomy. Please try again!\n";
}

exit;
