#!/usr/bin/perl
use strict;
use warnings;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;
use Bio::TreeIO;
use Getopt::Std;
use Cwd 'abs_path';

sub usage {
    my $msg = shift;
    print "example: perl getSearchTaxa.pl -i genome_dir -b 0.00005 -h 0.00005 -r 10 -n mammalia -t taxonomy -o searchList.txt\n";
    print "-i\tFolder contains all search species (e.g. genome_dir)\n";
    die $msg."\n";
}

# global variables
our($opt_i,$opt_b,$opt_h,$opt_r,$opt_n,$opt_t,$opt_o);
getopts('i:b:h:r:n:t:o:');

# sanity checks;
my $genome_dir = ($opt_i) ? $opt_i : usage("ERROR: No input folder given\n");
my $eval_blast = ($opt_b) ? $opt_b : usage("ERROR: No eval_blast given\n");
my $eval_hmmer = ($opt_h) ? $opt_h : usage("ERROR: No eval_hmmer given\n");
my $eval_relaxfac = ($opt_r) ? $opt_r : usage("ERROR: No eval_relaxfac given\n");
my $group = ($opt_n) ? $opt_n : usage("ERROR: No group given\n");
my $idx_dir = ($opt_t) ? $opt_t : usage("ERROR: No taxonomy dir given\n");
my $output = ($opt_o) ? $opt_o : usage("ERROR: No output given\n");

open(OUT, ">$output") || die "Cannot create $output\n";
my $groupNode;
my %taxa;
my $db;

if($group ne "all") {
    $db = Bio::DB::Taxonomy->new(-source    => 'flatfile',
        -nodesfile => $idx_dir . 'nodes.dmp',
        -namesfile => $idx_dir . 'names.dmp',
        -directory => $idx_dir);
    checkGroup($group);
    # get tree
    %taxa = getTaxa($genome_dir);
    my $tree = getTree();
    my $final_eval_blast = $eval_blast*$eval_relaxfac;
    my $final_eval_hmmer = $eval_hmmer*$eval_relaxfac;
    if($groupNode) {
        foreach($tree->get_nodes()) {
            if($_->id == $groupNode->id) {
                $groupNode = $_;
            }
        }
        $tree->set_root_node($groupNode);
    }
    foreach (get_leaves($tree)) {
        my $tmp = @{$_->name('supplied')}[0];
        print OUT $tmp,"\n";
    }
} else {
    %taxa = getTaxa($genome_dir);
    foreach my $tax (keys %taxa) {
        print OUT $tax,"\n";
    }
}
exit;

sub checkGroup {
    my ($group) = $_[0];
	my $node = $db->get_taxon(-name => $group);
	if($node) {
		$groupNode = $node;
	} else {
		print "Your selected group " . $group . " was not found in the taxonomic tree... TERMINATING\n";
		exit;
	}
}

sub getTaxa {
    my ($genome_dir) = $_[0];
    ## removal of misplaced files in genome_dir
    if (-e "$genome_dir/query.sql"){
        unlink("$genome_dir/query.sql");
    }
    if (-e "$genome_dir/@@.fa"){
        unlink("$genome_dir/@@.fa");
    }
    my @taxonlist = `ls $genome_dir`;
    chomp @taxonlist;
    for (my $i = 0; $i < @taxonlist; $i++) {
        my ($taxon_name, $ncbi_id, $src_id) = split /@/, $taxonlist[$i];
        if (!$src_id) {
            $src_id = '';
        }
        $taxon_name = $taxonlist[$i];
        $taxa{$taxon_name} = $ncbi_id;
    }
	my $hashcount = keys(%taxa);
	return(%taxa);
}

sub getTree {
	# the full lineages of the species are merged into a single tree
	my $tree;
	foreach my $key (sort {lc $a cmp lc $b} keys %taxa) {
		my $node = $db->get_taxon(-taxonid => $taxa{$key});
		if (!defined $node){
			print "ISSUE in sub getTree. No correspodence found in taxonomy file for $key and taxid $taxa{$key}. Skipping...\n";
			next;
		}
		else {
			$node->name('supplied', $key);
			if($tree) {
				$tree->merge_lineage($node);
			}
			else {
				$tree = Bio::Tree::Tree->new(-verbose => $db->verbose, -node => $node);
			}
		}
	}
	return $tree;
}

sub get_leaves {
	my $tree = $_[0];
	my $delFlag = 0;
	if(defined($_[1])){
		$delFlag = $_[1];
	}

	my $node = $tree->get_root_node;
	my @leaves;
	my @children = ($node);
	for (@children) {
		push @children, $_->each_Descendent();
	}
	for (@children) {
		push @leaves, $_ if defined($_->name('supplied'));
	}
	# if the tree is set to be deleted
	if ($delFlag){
		@leaves = qw();
		return @leaves;
	}else{
		return @leaves;
	}
}
