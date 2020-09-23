#!/usr/bin/perl
use strict;
use File::Basename;
use lib dirname(__FILE__);
use Getopt::Long;
use Bio::Perl;
use File::Copy;

# PROGRAMNAME: translate.pl

# AUTHOR: INGO EBERSBERGER, ingo.ebersberger@univie.ac.at

# PROGRAM DESCRIPTION:

# DATE: Tue May 12 14:03:34 CEST 2009


# DATE LAST MODIFIED: 03.11.2010: Bug fix suggested by Todd Oakley.
# BUG -- BIOPERL GUESSES PROTEIN FILE FORMAT WHEN AMBIGUITY CODES ARE PRESENT
# CAUSING AN ERROR IN THE TRANLATE_6 FRAMES, WHICH INTERRUPTS ALL TRANSLATION -- THO

## Last modified: 10.01.2014
## added option -outpath
######################## start main #############################
my $help;
my @out;
my @estout;
my $infile;
my $trunc = 1;
my $outfile = "translate_tc.out";
my $outpath = '.';
my $limit = 20; ## this sets the maximum length for the sequence identifier. If sequence identifier are
## too long, then one can run into troubles with the parsing of the hmmsearch results.
#########
my $usage = "Name:\n\ttranslate.pl\n
Synopsis:\n\ttranslate_tc5.pl [-infile=FILE] [options] [-outfile=FILE]\n
Description:\n\tThis program takes a batch fasta-file with DNA
\tsequences as an input and translates the individual DNA sequences in
\tall six reading frames.
\t-infile: provide the relative or absolute path of the infile\n
\t-outfile: provide the relative or absolute path of the outfile
\tDefault is: translate_tc.out\n
\t-outpath: provide the path to the
\toutfile. Default is '.'\n
\ttrunc: set -trunc=0 to prevent truncation of the sequence header (see below).
\t-h: prints this help-message\n
NOTE: if the seq-id (everything up to the first [[:space:]]) contains a '|' everything between the '>' and the '|' will be taken as seq-id. Otherwise, the entire seq-id will be used. You can change this behavior by setting -trunc=0\n
NOTE: the script as an automated routine to check for unique sequence names in the input file. This may lead to cases where the $trunc value is overruled and additionally part of the sequence description may be included.";
##########

GetOptions (
    "h" => \$help,
    "infile=s" => \$infile,
    "outfile=s" => \$outfile,
    "outpath=s" => \$outpath,
    "trunc=s" => \$trunc);
if ($help) {
	print "$usage";
	exit;
}
if (-e "$outfile") {
    print LOG "an outfile $outfile already exists. Renaming to $outfile.old\n\n";
    my $newname = "$outfile.old";
    rename($outfile, $newname);
}
#my @seq_object = read_all_sequences($infile, 'fasta');

open (LOG, ">>$outpath/hamstrsearch.log") or warn "could not open logfile for writing\n";
print LOG "\n### TRANSLATE.PL: \n";

### changes suggested by Todd Oakley
my $tempseqio;
$tempseqio = Bio::SeqIO->new( '-file' => $infile, '-format' => 'fasta');
my @seq_object;

while( my $seq = $tempseqio->next_seq() ) {
     $seq->alphabet('dna');
     push(@seq_object,$seq);
}
### End changes Todd Oakley

## determine whether the seq-ids are unique given the chosen value for $trunc
my ($message, $cont, $check) = &checkIds();
if ($cont == 1) {
    ## the check for unique identifiers has failed and the programm is exiting
    print LOG "$message\n";
    close LOG;
	exit;
}
else {
    print LOG "All sequence identifier are unique!\n";
    if ($check == 2) {
	my $newname = "$infile.original";
	rename($infile, $newname);
	print LOG "Sequence description was needed to make seq-id unique. The original version of the infile was stored in $infile.original\n";
    }
    for (my $j = 0; $j < @seq_object; $j++) {
	my $finalid = $seq_object[$j]->{finalid};
	my $estseq = $seq_object[$j]->seq;
	my $inid = $seq_object[$j]->display_id;
	my @all_trans = Bio::SeqUtils->translate_6frames($seq_object[$j]);
	for (my $i = 0; $i < @all_trans; $i++) {
	    my $count = $i+1;
	    my $id = $all_trans[$i]->display_id;
	    my $seq = $all_trans[$i]->seq;
	    $id =~ s/$inid/$finalid/;
	    $id =~ s/-[0-9][RF]/_RF$count.0/;
	    push @out, ">$id\n$seq";
	}
	push @estout, ">$finalid\n$estseq";
	if ($j%100 == 0) {
	    print "$j Sequences processed\n";
	    open (OUT, ">>$outpath/$outfile") or die "failed to open outfile\n";
	    print OUT join "\n", @out;
	    print OUT "\n";
	    @out = qw();
	    close OUT;
	    if ($check == 2) {
		## part of the description was added to the seq-id
		open (OUT, ">>$infile");
		print OUT join "\n", @estout;
		print OUT "\n";
		@estout = qw();
	    }
	}
    }
    open (OUT, ">>$outpath/$outfile") or die "failed to open outfile\n";
    print OUT join "\n", @out;
    print OUT "\n";
    @out = qw();
    close OUT;
    if ($check == 2) {
	## part of the description was added to the seq-id
	open (OUT, ">>$infile");
	print OUT join "\n", @estout;
	print OUT "\n";
	close OUT;
	@estout = qw();
    }
}
close LOG;
exit;
########################## start sub ################
sub checkIds {
    my $message;
    my $check = 1;
    my $cont = 1;
    my $counter;
    ## Everything up to the first whitespace
    ## in the fasta header will be taken as sequence id by bioperl. If this
    ## id contains a '|' and $trunc is set to 1 (default), the ids may no longer
    ## be unique. This will be checked and if necessary the id will not be truncated
    ## for $check == 0, the truncated version of the id will be checked (only if $trunc == 1)
    ## for $check == 1, the complete id will be checked
    ## for $check == 2, the first 20 characters of the concatenated id and description
    ## will be checked
    if ($trunc == 1) {
	$check = 0;
    }

    while ($check < 3 and $cont == 1) {
	$cont = 0;
	for (my $i=0; $i < @seq_object; $i++) {
	    my $id = $seq_object[$i]->display_id;
	    $id =~ s/(.{0,$limit}).*/$1/;
	    if ($check == 0) {
		$id =~ s/|.*//;
	    }
	    elsif ($check == 2) {
		$id = $id . '_' . $seq_object[$i]->desc;
		$id =~ s/(.{0,$limit}).*/$1/;
	    }
	    if (defined $counter->{$id}) {
		if ($check == 0) {
		    $message = "trying next without truncating the id";
		}
		elsif ($check == 1) {
		    $message = 'trying next to include sequence description';
		}
		else {
		    $message = "Sequence identifier are not unique, using the first 20 characters. Aborting...";
		}
		print LOG "sequence ids are not unique in the file $infile, $message. The offending identfier is $id\n\n";
		$check ++;
		$cont = 1;
		$counter = undef;
		last;
	    }
	    else {
		$counter->{$id} = 1;
		$seq_object[$i]->{finalid} = $id;
	    }
	}
    }
    ## return the value of $cont. If this is 1, then the sequence id check has failed.
    return($message, $cont, $check);
}
