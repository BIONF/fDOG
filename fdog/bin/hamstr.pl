#!/usr/bin/perl
use strict;
use Getopt::Long;
use Parallel::ForkManager;
use Bio::SearchIO;
use Bio::Search::Hit::BlastHit;
use Bio::SeqIO;
use Bio::Align::ProteinStatistics;
use Bio::AlignIO;
use Term::Cap;
use POSIX;
use Cwd;
use Cwd 'abs_path';
use Statistics::R;
use File::Basename;
use lib dirname(__FILE__);
use run_genewise_hamstr;

# PROGRAMNAME: hamstr.pl

# Copyright (C) 2009 INGO EBERSBERGER, ingo.ebersberger@univie.ac.at
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 3 of the License
# or any later version.

# This program is distributed in the hope that it will be useful
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program; If not, see http://www.gnu.org/licenses

# PROGRAM DESCRIPTION: HaMStR is a program for targeted ortholog search in both EST/RNAseq
# and protein sequence data.

# DATE: Wed Dec 19 10:41:09 CEST 2007

# PROGRAM HISTORY
##23. 07. 2010: found a bug in the extraction of the
## hmm hit sequence from the sequnence_file. A end-of-line char was missing.

##09.08.2010: added the option to choose the new blastp program from ncbi. Just comment
##out line 45 in the script and uncomment line 46. Note, in order to make this work I have
##to slightly modify the blast output since otherwise it will not be parsed by the Bioperl
##Blast parser. Currently this is a pretty dirty $sedprog hack. It will also take care of removin
##the string lcl| that is added in some instances to the hit id in the blast output.

## I added the option that one can now provide a comma-separated string of phmm names as an
## argument for the option -hmm

## 08.03.2011:
## 1) BUG-FIX: Hamstr will now remove automatically newlines from the input sequence file
## 2) BUG-FIX: The sequence header remains now the same whether or not the flag -representative
## has been chosen.

## 10.04.2011
## 1) added some information to the log file.

## 20.05.2011
## 1) BUG-FIX: The grep for the EST sequence in the sub-routine predictORF received also a hit when
## the search pattern was only a substring of the EST sequence identifier. In some cases the wrong EST
## was then used to predict the ORF. This has been fixed.

## 30.05.2011
## 1) Extension: a command line option -longhead has been added. The user can now specify that the
## full sequence id including whitespaces will considered throughout the hamstr search. Note, the
## whitespaces will be replaced by the string specified in the variabel $idsep.
## 2) Modification from the bug fix from 20.05.2011. In the grep for the original EST it is no longer
## necessary that the search string and the EST id are identical over their entire length. Instead the
## search string may be a prefix of the EST id ending with a whitespace.

## 27.06.2011
## 1) Extension: I added the option to run a true reciprocal best hit search. Only the best hit from the
## hmmer search is used to check for reciprocity.

## 06.12.2011
## 1) Extension: I added the option -hit_limit to set the number of hmmsearch hits that HaMStR uses for
## the re-blast.

## 10.02.2012
## 1) Extension: I added checks for the appropriate hmmsearch version (HMMER 3) and for genewise and
## its environmental variable WISECONFIGDIR.
## 2) Bug fix in the -rbh option.

## 11.09.2012
## 1) Bug fix: -hitlimit, even if not set explicitely has been invoked resulting in a more stringent
## behaviour of HaMStR. This has been fixed resulting in longer run-times.
## 18.12.2012
## 1) Bug fix: There was a bug in the CDS extraction for reverse complemented
## sequences. A new line was moved to the beginning of the sequence
## leading to index errors.

## 18.12.2013
## 1) Bug fix: I have now adapted the script such that it no longer requires the default directory structure
## 2) Extension: Hamstr is now capable of identifying co-orthologs (sub routine IdentifyCoorthologsProt)
## 3) The re-blast for EST sequences is now a BlastX solving the problem of duplicated output for contigs with a pHMM
## hit in more than one frame.

## 08.01.2014
## Extension: check for co-orthology between ref-protein and best blast hit in cases both are not identical.
## Bug fix: option -rbh was disfunctional due to a missing function in new sub-routine parseHmmer3pm
## Bug fix: sortRef actually did not sort anything as it was accessing the unsorted file

## 09.01.2014
## Extension: Add the possibility to sort the hits by the hmmersearch Score rather than an alignment score
## This will make the best hmmersearch hit surviving the re-blast automatically the representative

## 10.01.2014
## Bug fix (minor): modified option -outpath to accept non-default value
## modification of the translate_tc call.

## 17.01.2014
## Bug fix (minor): added the option --anysymbol to the mafft-linsi command to avoid crash when protein sequences
## contain non-standard amino acids (U = selenocystein)

## 14.02.2014
## Extension: added the option to use ublast rather than blast for the reciprocity check

## 25.02.2014
## Extension: added the option -reuse. By default, old results will now be deleted unless flag -reuse has been set
## Extension: added syntax required for running fact implemented into Hamstr2.0. Option -fact does not occur in help
## as this works only in context with Hamstr2.0
## Extension: Hamstr now outputs a results summary at the end of the run.
## Extension: added the option -cleartmp to faciliate automatic removal of the tmp-dir prior to the hamstr run

## 05.03.2014
## Bug fix (minor): Variable $grepprog was not used throughout the script. In some routines 'grep' was hard coded.
## On MAC OS this could lead to unwanted call of BSD grep resulting in an error during re-blast.

## 16.03.2014
## Bug fix (major): There was a problem in translating ESTs in the correct frame. This has been fixed.
## Modification (minor): The alignment positions together with the score are no longer sorted externally.

## 02.04.2014
## Bug fix (minor): Flag $runublast was not changed to 1 when configuring hamstr
## with the ublast option.

## 05.08.2014
## Exentsion: Update to version 13.2.3. New features include the option to run hamstr in silent mode and the option
## to parallelize the hamstr search for individual core orthologs using the Parallel::ForkManager package

## 14.08.2014
## Bug fix (minor): corrected typo in sub routine call 'printOUT'

## 31.07.2015
## Extension: Update to version 13.2.4. New feaure provides the option to entirely remove intron sequences and incomplete codons
## from transcripts

## 03.07.2015
## Minor extension: Selected behavior with respect to introns in transcripts will be printed to hamstrsearch.log

## 05.07.2015
## Minor bug fix: A no-blast hit was not reported properly in the sub routine check4reciprocity resulting in rare
## cases in the acceptance of a spurious ortholog.

## 14.08.2015
## Change of output file naming. Upon selection of the strict option the reference species is no longer appended
## to the output file name.

## 30.01.2016
## Minor changes including the better integration of the onseq.pl script. Among others the blast files are no longer
## expected to have the '_prot' appendix.

## 12.02.2016
## Minor bug fix: In some instances the representative protein was not chosen correctly due to a bug in the subroutine
## sortRef. Analyses of transcript data are not affected at all.

## 19.12.2017
## Extension: HaMStR can now automatically determine the hit limit up to which candidates from the intial
## hmm search are evaluated as potential orthologs. Two options are available, either an hmm score driven
## cutoff determination, or alternatively, a lagPhase-based estimator.

## 02.02.2018
## Bug fix (solved): using grep within the checkcoorthologsref routine could cause an incomplete alignment of the reference gene,
## the candidate ortholog and the best blast hit. The resulting distance (kimura) calculation may caused an overoptimistic
## acceptence of co-ortholgy relations. The bug onyl occured while using the option checkCoOrthologsRef.
## HaMStR keeps original gene sets in FASTA format and *.fa.mod will link to the original FASTA file (no linebreaks within a sequence).

## 28.02.2018
## Minor Bug fix (solved): HaMStR is not longer asking infinite times for the replacement of already existing output files.
## Minor Bug fix (solved): Backward compatibility extended. Naming of reference fasta files: (*.fa and *_prot.fa)

## 20.07.2019
## fixed the issue of long proteins with best total hmm bit score but very poor domain scores. Allow now the option to sort
## hmmsearch output according to the best domain bit score. The current routine assumes that neither query nor
## hit has whitespaces in their names

## 14.04.2020 (Vinh)
## Bug fix (solved): existing symbolic link cannot recognized while checking the reference fasta file

## 10.07.2020 (v13.2.12 - vinh) solved problem when gene ID contains PIPE
## 13.07.2020 (v13.3.0 - vinh) solved problem when gene ID contains PIPE
## 22.07.2020 (v13.4.0 - vinh) moved tmp blast files to output folder and delete them when finished
## 01.12.2020 (v13.4.1 - vinh) add silent option to muscle for checkCoOrthologsRef
## 21.01.2021 (v13.4.2 - vinh) fiexed bug when refspec has "dot" in its name

######################## start main ###########################################
my $version = "HaMStR v.13.4.2";
######################## checking whether the configure script has been run ###
my $configure = 0;
if ($configure == 0){
	die "\n\n$version\n\nPLEASE RUN setup1s BEFORE USING HAMSTR\n\n";
}
########## EDIT THE FOLLOWING LINES TO CUSTOMIZE YOUR SCRIPT ##################
my $prog = 'hmmsearch'; #program for the hmm search
my $eval = 1; # default evalue cutoff for the hmm search
my $sedprog = 'sed';
my $grepprog = 'grep';
my $readlinkprog = 'readlink';
my $alignmentprog = 'clustalw';
my $alignmentprog_co = 'muscle';
########## EDIT THE FOLLOWING TWO LINES TO CHOOSE YOUR BLAST PROGRAM ##########
my $blast_prog = 'blastp';
my $filter = 'F'; # low complexity filter switch. Default 'on'. Set of 'F' to turn off permanently.
my $eval_blast = 10; # default evalue cutoff for the blast search
########## EDIT THE FOLLOWING LINES TO MODIFY DEFAULT PATHS ###################
my $path = abs_path(dirname(__FILE__));
$path =~ s/\/bin//;
my $hmmpath = "$path/core_orthologs"; #path where the hmms are located
my $blastpath = "$path/blast_dir"; #path to the blast-dbs
my $outpath = '.';
my $idsep = '__'; #character used to replace whitespaces in the sequence header with (flag -longhead)
my $hmm_dir = 'hmm_dir';
my $fa_dir  = 'fa_dir';
##############################
my $termios = new POSIX::Termios; $termios->getattr;
my $ospeed = $termios->getospeed;
my $t = Tgetent Term::Cap { TERM => undef, OSPEED => $ospeed };
my ($norm, $under, $bold) = map { $t->Tputs($_,1) } qw/me md us/;

############################## Variables ##############
my $fileobj;
## The main variable storing most of the results;
## $fileobj->{$taxon}->{prot}->[$hitcounter]
## $fileobj->{$taxon}->{ids}->[$hitcounter]
## $fileobj->{$taxon}->{cds}->[$hitcounter]
## $fileobj->{$taxon}->{hmmscore}->[$hitcounter]
#######################################################
my $pid = $$;
my $help;
my $debug;
my $seq2store_file='';
my $cds2store_file='';
my $hmm;
my @hmms;
my $fa;
my $fafile;
my @seqs2store;
my @cds2store;
my $dbpath;
my $dboutpath;
my $ep2eg;
my $dbfile_base;
my $aln;
my $idfile;
my $taxon_check = 0;
my $hmmset;
my $show_coreortholog_sets;
my $hmmsearch_dir;
my $dbfile; # the file hmmsearch is run against
my $dbfile_short;
my $taxon_file;
my $refspec_string;
my @refspec = qw();
my @primer_taxa;
my $refspec_name = '';
my $taxon_global;
my $fa_dir_neu = '';
my $gwrefprot;
my $seqtype;
my $align;
my $rep;
my $estflag;
my $proteinflag;
my $refseq;
my $strict;
my $relaxed;
my $refspec_final = '';
my $central;
my $concat;
my $seqs2store_file;
my $append;
my $longhead;
my $check = 1;
my @log = qw();
my $bhh;
my $hitlimit;
my $autoLimit;
my $scoreThreshold;
my $scoreCutoff = 10;
my $nonoverlappingCO;
my $algorithm = 'blastp';
my $frame;
my $checkCoRef;
my $sortalign;
my $check_genewise = 1;
my $outputfmt = 'blastxml';
my $fact;
my $runFACTparameter;
my $hmmcount;
my $reuse;
my $cleartmp;
my $ver;
my $silent;
my $cpu = 1;
my $force;
my $keepintron = 'k';
my $blastapp  = '';
my $blastdbend = '.pin';
######### ublast options #########
my $runublast = 1;
my $ublast = 0;
my $accel = 0.8;
#####determine the hostname#######
push @log, "VERSION:\t$version\n";
my $hostname = `hostname`;
chomp $hostname;
push @log, "HOSTNAME\t$hostname\n";
#################################
if (@ARGV==0) {
	$help = 1;
}
## help message
my $helpmessage = "
${bold}YOU ARE RUNNING $version on $hostname$norm

This program is freely distributed under a GPL.
Copyright (c) GRL limited: portions of the code are from separate copyrights

\n${bold}USAGE:${norm} hamstr -sequence_file=<> -hmmset=<> -taxon=<>  -refspec=<> [OPTIONS]

${bold}OPTIONS:$norm

${bold}REQUIRED$norm
-sequence_file=<>
		path and name of the file containing the sequences hmmer is run against.
-hmmset=<>
		specifies the name of the core-ortholog set.
		The program will look for the files in the default directory 'core-orthologs' unless you specify
		a different path via the option -hmmpath.
-refspec=<>
		sets the reference species. Note, it has to be a species that contributed sequences
		to the hmms you are using. NO DEFAULT IS SET! For a list of possible reference
		taxa you can have a look at the speclist.txt file in the default core-ortholog sets
		that come with this distribution. Please use the abreviations in this list. If you choose
		to use core-orthologs where not every taxon is represented in all core-orthologs, you
		can provide a comma-separated list with the preferred refspec first. The lower-ranking
		reference species will only be used if a certain gene is not present in the preferred
		refspecies due to alternative paths in the transitive closure to define the core-orthologs.
		CURRENTLY NO CHECK IS IMPLEMENTED!
		NOTE: A BLAST-DB FOR THE REFERENCE SPECIES IS REQUIRED!
-taxon
		You need to specify a default taxon name from which your ESTs or protein sequences are derived.
-est
		set this flag if you are searching in ESTs. Note, if neither the -est nor the -protein flag is set, HaMStR will
		guess the sequence type. If you select this flag, make sure to specify how to deal with introns retained in the
		ESTs. Check option -intron!
-protein
		set this flag if you are searching in protein sequences. Note, if neither the -est nor the -protein flag is set, HaMStR will
		guess the sequence type.

${bold}USING NON-DEFAULT PATHS$norm

-blastpath=<>
		Lets you specify the absolute or relative path to the blast databases. DEFAULT: $blastpath
-hmmpath=<>
		Lets you specify the absolute or relative path to the core ortholog set. DEFAULT: $hmmpath
-outpath=<>
		You can determine the path to the HaMStR output. Default: current directory.

${bold}ADDITIONAL OPTIONS$norm

-append
		set this flag if the output should be appended to the files *.out and *_cds.out. This becomes relevant when running
		hamstrsearch with individual hmms and you want to combine the results.
-central
		set this flag to store the modified infile in the same directory as the infile rather than in the output dir.
-checkCoorthologsRef
		If the re-blast does not identify the original reference protein sequence as best hit, HaMStR will check whether the best blast
		hit is likely a co-ortholog of the reference protein relative to the search taxon. NOTE: Setting this flag will substantially increase
		the sensitivity of HaMStR but most likely affect also the specificity, especially when the search taxon is evolutionarily only very
		distantly related to the reference taxon.
-cleartmp
		set this flag to remove existing tmp dir in the HaMStR output directory.
-concat
		set this flag if you want hamstr to concatenate sequences that align to non-overlapping parts of the reference protein.
		If you choose this flag, no co-orthologs will be predicted.
-cpu
		You can specify the number of parallel jobs in the HaMStR search. HaMStR uses the Parallel::ForkManager module for this purpose.
-eval_blast=<>
		This option allows to set the e-value cut-off for the Blast search. Default: 10
-eval_hmmer=<>
		This options allows to set the e-value cut-off for the HMM search.Default: 1
-filter=<T|F>
		Set this flag to F if the re-blast should be performed without low-complexity filtering. Default is T.
-force
		Setting this flag forces hamstr to overwrite existing output files (files ending with .out) without further asking.
-hit_limit=<>
		By default, HaMStR will re-blast all hmmsearch hits against the reference proteome. Reduce the number
		of hits for reblast with this option.
-autoLimit
		Setting this flag will invoke a lagPhase analysis on the score distribution from the hmmer search. This will determine automatically
		a hit_limit for each query.
-scoreThreshold
		Instead of setting an automatic hit limit, you can specify with this flag that only candidates with an hmm score no less
		than x percent of the hmm score of the best hit are further evaluated. Default is x = 10.
		You can change this cutoff with the option -scoreCutoff. Note, when setting this lag, it will be effective for
		both the core ortholog compilation and the final ortholog search.
-scoreCutoff=<>
		In combination with -scoreThreshold you can define the percent range of the hmms core of the best hit up to which a
		candidate of the hmmsearch will be subjected for further evaluation. Default: 10%.
-hmm
		Option to provide only a single hmm to be used for the search.
		Note, this file has to end with .hmm
-intron=<${bold}k${norm}eep|${bold}m${norm}ask|${bold}r${norm}emove>
		Specify how to deal with introns that may occur in transcript sequences. Default: keep - Introns will be retained in the transcript
		but will be identified by lower case letters.
-longhead
		Set this flag in the case your sequence identifier contain whitespaces and you whish to keep
		the entire sequence identifier throughout your analysis. HaMStR will then replace the whitespaces with
		a '__'. If this flag is not set, HaMStR will truncate the sequence
		Identifier at the first whitespace, however if and only if the sequence identifier then remain unique.
		NOTE: too long sequence headers (~ > 30 chars) will cause trouble in the hmmsearch as the program will truncate
		the output!
-nonoverlapping_cos
		If you set this flag, non-overlapping co-orthologs will be reported as well. NOTE: this flag is still experimental
-rbh
		set this flag if you want to use a reciprocal best hit criterion. Only the highest scoring
		hit from the hmmer search will be used for re-blast.
-relaxed
		set this flag if the reciprocity criterion is fulfilled when the re-blast against
		any of the primer taxa was successfull. Note that setting this flag will substantially decrease the
		stringency of the ortholog assignment with the consequence of an increased number of false positives.
-representative
		From all sequences that fulfill the reciprocity criterion the one showing the highest similarity to the
		core ortholog sequence in the reference species is identified and selected as representative.
-reuse
		Set this flag if you want to prevent HaMStR from overwriting previous results.
-show_hmmsets
		setting this flag will list all available core ortholog sets in the specified path. Can be combined with -hmmpath.
-silent
		Supresses (almost) all print statements to the screen.
-debug
		Get some additional meta information as print out to the screen.
-sort_global_align
		Setting this flag will tell hamstr to sort ortholog candidates according to their global alignment score to the reference
		sequence rather than according to the score they have achieved in the hmmer search (local). NOTE: In the case of searching
		EST data this flag is automatically set.
-strict
		Set this flag if the reciprocity criterion is only fulfilled when the re-blast against
		all primer taxa was successfull
-aligner
		Choose between muscle or mafft-linsi for the alignment of multiple sequences. DEFAULT: muscle
		\n\n";

GetOptions (
	"append" => \$append,
	"autoLimit"	=> \$autoLimit,
	"aligner=s"	=> \$alignmentprog_co,
	"blastpath=s" => \$blastpath,
	"checkCoorthologsRef" => \$checkCoRef,
	"concat" => \$concat,
	"cpu=s"  => \$cpu,
	"central" => \$central,
	"debug" => \$debug,
	"est" => \$estflag,
	"eval_blast=s" => \$eval_blast,
	"eval_hmmer=s" => \$eval,
	"fasta_file=s" => \$fafile,
	"filter=s" => \$filter,
	"force"    => \$force,
	"h"        => \$help,
	"hit_limit=s" => \$hitlimit,
	"hmm=s"    => \$hmm,
	"hmmset=s" => \$hmmset,
	"hmmpath=s" => \$hmmpath,
	"intron=s" => \$keepintron,
	"longhead" => \$longhead,
	"nonoverlapping_cos" => \$nonoverlappingCO,
	"outpath=s" => \$outpath,
	"protein"=> \$proteinflag,
	"rbh" => \$bhh,
	"refspec=s" => \$refspec_string,
	"relaxed" => \$relaxed,
	"representative" => \$rep,
	"reuse"	=> \$reuse,
	"sequence_file=s" => \$dbfile,
	"scoreCutoff=s" => \$scoreCutoff,
	"scoreThreshold" => \$scoreThreshold,
	"show_hmmsets" => \$show_coreortholog_sets,
	"silent" => \$silent,
	"sort_global_align" => \$sortalign,
	"strict" => \$strict,
	"taxon_file=s" => \$taxon_file,
	"taxon=s"  => \$taxon_global,
	"ublast"	=> \$ublast,
	"v"	=> \$ver,
	"accel=s" => \$accel,
	"fact" => \$fact,
	"cleartmp" => \$cleartmp
);

if ($help) {
	print $helpmessage;
	exit;
}
elsif($ver){
	print "$version\n";
	exit;
}

## 1) check if all information is available to run HaMStR
($check, @log) = &checkInput();
if ($check == 0) {
	print "\n\n${bold}There was an error running $version$norm\n\n";
	print join "\n", @log;
	exit;
}
else {
	open (OUT, ">$outpath/hamstrsearch.log") or die "could not open logfile\n";
	print OUT join "\n", @log;
	close OUT;
}
my $tmpdir = "$outpath/tmp";

### read in of the core-ortholog sequences
my $co_seqs = parseSeqfile("$fafile");

## initialize the forking procedure
my $pm = new Parallel::ForkManager($cpu);

## collect all the entries of the final output file
#my ($spid, $exit_code, $ident, $exit_signal, $core_dump, $data);
#$pm->run_on_finish(sub {
#	($spid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
#	$core_dump = undef;
#	if ($seqderef){
#		push @seqs2store, @$seqderef;
#		if ($estflag) {
#			my $estderef = $data->[1];
#			push @cds2store, @$estderef;
#		}
#	}
#});

## 2) loop through the hmms
## process each hmm file separately
$hmmcount = scalar(@hmms);

for (my $i = 0; $i < @hmms; $i++) {
	my $pid = $pm->start and next;
	my $localid = $$;
	$frame = undef;
	$fileobj = undef;
	my @seqs = qw();
	my @newseqs = qw();## var to contain the sequences to be added to the orthologous cluster
	my @newcds = qw();
	my $hmm = $hmms[$i];
	printOUT("Processing $hmm\n");
	my $hmmout = $hmm;
	$hmmout =~ s/\.hmm/\.out/;
	## 3) run the hmm search
	if (!(-e "$hmmsearch_dir/$hmmout")) {
		printOUT("\n\nnow running $prog using $hmm\n");
		my $hmmOutFile = "$hmmsearch_dir/$hmmout";
		my $hmmModel = "$hmm_dir/$hmm";
		my $hmmInfile = "$dboutpath/$dbfile";
		`$prog --noali --tblout \"$hmmOutFile\" -E $eval \"$hmmModel\" \"$hmmInfile\"` or die "Problem running hmmsearch as $prog --noali --tblout \"$hmmOutFile\" -E $eval \"$hmmModel\" \"$hmmInfile\". No output $hmmsearch_dir/$hmmout\n";
	}
	else {
		printOUT("an hmmresult $hmmout already exists. Using this one!\n");
	}

	## 4) process the hmm search result
	my $hitcount = 0;
	## 4a) loop through the individual results
	## now the modified version for hmmer3 comes
	my $hitlimit_local = $hitlimit;
	my ($query_name, $results, $hitlimit_local, $criticalValue) = parseHmmer4pm($hmmout, $hmmsearch_dir);
	if (! $results) {
		printOUT("no hit found for $query_name\n");
		$pm->finish;
		next;
	}
	## Automatic hit limit information
	if (defined $autoLimit) {
		printDebug("Automatic cutoff estimation via a lag Phase analysis was selected. Estimated lag point is $criticalValue. Limiting the number of hits for the evaluation from " . scalar(@$results) . " to $hitlimit_local");
	}
	elsif (defined $scoreThreshold) {
		printDebug("Automatic cutoff estimation via a minimal score was selected. Cutoff: $scoreCutoff percent of the best hmm score. Hits with an hmm score below $criticalValue are not considered. Limiting the number of hits for the evaluation from " . scalar(@$results) . " to $hitlimit_local");
	}
	##
	printOUT("Results for $query_name\n");
	my ($check, $refspec_final) = &determineRefspecFinal($query_name, @refspec);
	if ($check == 0) {
		die "error in retrieving refspec data\n";
	}
	if (!defined $hitlimit_local or $hitlimit_local > scalar(@$results)) {
		$hitlimit_local = scalar(@$results);
	}
	for (my $k = 0; $k < $hitlimit_local; $k++) {
		my $hitname = $results->[$k]->{id};
		my $hithmmscore = $results->[$k]->{hmmscore};
		printOUT("$hitname\n");
		my $keep = 0;
		my $hitseq = '';
		$refseq = '';
		## 4b) test for the reciprocity criterion fulfilled
		($keep, $hitseq, $frame)  = &check4reciprocity($localid, $query_name, $hitname, $refspec_final, @refspec);
		if ($keep == 1) {
			## blast search with the hmm hit identifies the core-ortholog sequence of the reference species
			my $taxon = $taxon_global;
			## put the info about the hits into an object for later post-processing
			### HERE COMES THE NEW STUFF THAT DEALS WITH THE DIFFERENT POSSIBILITIES: STRICT, RELAXED OR WHATEVER...
			$fileobj = &determineReferences ($localid, $fileobj, $taxon, $refspec_final, $hitname, $hithmmscore, $hitseq, $hitcount);
			$hitcount++;
		}
		else {
			printOUT("Reciprocity not fulfilled!\n\n");
		}
	}
	## 5) do the rest only if at least one hit was obtained
	if (defined $fileobj) {
		## 5a) if the hits are derived from ESTs, get the best ORF
		if ($estflag) {
			$fileobj =  &predictORF($frame);
		}
		&processHits($localid, $fileobj);
		if (!$rep and !$concat) {
			## identify co-orothologs only for protein sequences. This adds a key 'coorthologs' to the $fileobj->{$taxon} that
			## holds the index values for the sequences in the $fileobj->{$taxon}->{ids} and the corresponding {prot} array ref
			## that made it into the co-ortholog field
			&identifyCoorthologsProt($localid, $taxon_global);
		}
		## 6) prepare the output
		my @taxa = keys(%$fileobj);
		for (my $i = 0; $i< @taxa; $i++) {
			push @newseqs, ">$query_name|$fileobj->{$taxa[$i]}->{refspec_final}|$taxa[$i]|$fileobj->{$taxa[$i]}->{refid}|1";
			push @newseqs, $fileobj->{$taxa[$i]}->{refprot};
			if ($estflag) {
				push @newcds, ">$query_name|$fileobj->{$taxa[$i]}->{refspec_final}|$taxa[$i]|$fileobj->{$taxa[$i]}->{refid}|1";
				push @newcds, $fileobj->{$taxa[$i]}->{refcds};
			}
			if (!$rep and !$concat){
				## print the remaining sequences only when the -representative option has not been chosen.
				my $coorthologsobj = $fileobj->{$taxa[$i]}->{coorthologs};
				my $idobj = $fileobj->{$taxa[$i]}->{ids};
				my $protobj = $fileobj->{$taxa[$i]}->{prot};
				my $cdsobj  = $fileobj->{$taxa[$i]}->{cds};
				my $refspecobj = $fileobj->{$taxa[$i]}->{refspec};
				for (my $j = 0; $j < @$coorthologsobj; $j++) {
					my $index = $coorthologsobj->[$j];
					push @newseqs, ">$query_name|$refspecobj->[$index]|$taxa[$i]|$idobj->[$index]|0";
					push @newseqs, $protobj->[$index];
					if ($estflag) {
						push @newcds, ">$query_name|$refspecobj->[$index]|$taxa[$i]|$idobj->[$index]|0";
						push @newcds, $cdsobj->[$index];
					}
				}
			}
			my $refs = $co_seqs->{$query_name};
			for (keys %$refs) {
				my $line = ">$query_name|$_|" . $refs->{$_}->{seqid} . "\n" . $refs->{$_}->{seq};
				push @seqs, $line;
			}
			chomp @seqs;
			printOUT("\n");
			@seqs = (@seqs, @newseqs);
			open (OUT, ">$fa_dir_neu/$query_name.fa");
			print OUT join "\n", @seqs;
			print OUT "\n";
			close OUT;
			if ($estflag) {
				open (OUT, ">$fa_dir_neu/$query_name.cds.fa");
				print OUT join "\n", @newcds;
				close OUT;
			}
			open (OUT, ">>$seqs2store_file") or die "failed to open output file\n";
			if ($estflag){
				open (OUT2, ">>$cds2store_file") or die "failed to open output file for cds\n";
			}
			for (my $i = 0; $i < @newseqs; $i+= 2) {
				my $line = $newseqs[$i] . "|" . $newseqs[$i+1];
				$line =~ s/>//;

				print OUT $line;
				print OUT "\n";
				#			push @seqs2store, $line;
				if ($estflag) {
					my $cdsline = $newcds[$i] . "|" . $newcds[$i+1];
					$cdsline =~ s/>//;
					print OUT2 $cdsline;
					print OUT2 "\n";
					push @cds2store, $cdsline;
				}
			}
			close OUT;
			close OUT2;
		}
	}
	if (@seqs2store > 0) {
		my $seqref = \@seqs2store;
		my $estref = \@cds2store;
		$pm->finish;
	}
	else {
		$pm->finish;
	}
}

$pm->wait_all_children;

### The following bit of code has been out-commented as shared memory between forked child
### processes does not exist. The handing back of return values from the child to the parent
### does work, however leads to memory problems.
### all HaMStR searches have been completed and all children have finished. Do the output

#if (@seqs2store > 0) {
#  if ($append) {
#   open (OUT, ">>$seqs2store_file") or die "failed to open output file\n";
# }
# else {
#   open (OUT, ">$seqs2store_file") or die "failed to open output file\n";
# }
# print OUT join "\n", @seqs2store;
# print OUT "\n";
# close OUT;
# if ($estflag) {
#   if ($append) {
#     open (OUT, ">>$cds2store_file") or die "failed to open output file\n";
#   }
#   else {
#   open (OUT, ">$cds2store_file") or die "failed to open output file\n";
#   }
#   print OUT join "\n", @cds2store;
#   print OUT "\n";
#   close OUT;
# }
#}
###########################################################################################

my $orthologs = 0;
if (-e $seqs2store_file) {
	$orthologs = `less $seqs2store_file |wc -l`;
	if ($fact){
		## starting funFACT.pl <user_defined_taxonName> <comma-separated string_of_all_used_reference_species>
		system("perl runFact.pl $runFACTparameter $outpath $blastpath $taxon_global $refspec_string");
	}
}
else {
	printOUT("no hits found\n\n");
}
### WRAP UP #####
my $fa_dir_neu_tmp = $fa_dir_neu; $fa_dir_neu_tmp =~ s/\|/\\\|/g;
my $ortholog_groups = `ls $fa_dir_neu_tmp |$grepprog -v 'cds.fa' |wc -l`;
my $hmmsearch_dir_tmp = $hmmsearch_dir; $hmmsearch_dir_tmp =~ s/\|/\\\|/g;
my $hmmsearched = `ls $hmmsearch_dir_tmp |wc -l`;
chomp ($ortholog_groups, $hmmsearched, $orthologs);

if (!defined $silent) {
	print "\n\n
####HaMStR completed!#########
Results of HaMStR search in $taxon_global
Number of core_orthologs searched: $hmmcount
Number of core_orthologs with hmmsearch output: $hmmsearched
Number of ortholog_groups extended: $ortholog_groups
Number of orthologous sequences: $orthologs
##############################\n\n";
} else {
	# print "$taxon_global done\n";
}
exit;


##################### start sub ###############

####### checkInput performs a number of checks whether sufficient information
### and all data are available to run HaMStR
sub checkInput {
	######### check a number of flags that only serve for providing the user with some information
	if (defined $show_coreortholog_sets) {
		## Do nothing but just list all available core ortholog sets in $hmmpath
		my @coresets = (`ls $hmmpath`);
		chomp @coresets;
		if (scalar(@coresets > 0)){
			print "\n${bold}THE FOLLOWING CORE ORTHOLOG SETS ARE AVAILABLE IN $hmmpath:${norm}\n\n";
			for (my $i = 0; $i < @coresets; $i++){
				my @available = qw();
				my @unavailable = qw();
				print "\n${bold}$coresets[$i]${norm}\n\n";
				my @refspec = `head -n 20 $hmmpath/$coresets[$i]/$coresets[$i].fa |$grepprog '>' |cut -d '|' -f 2 |sort |uniq`;
				chomp @refspec;
				for (my $j = 0; $j < @refspec; $j++){
					if (-e "$blastpath/$refspec[$j]"){
						push @available, "\t$refspec[$j]";
					}
					else {
						push @unavailable, "\t$refspec[$j]";
					}
				}
				print "\tAvailable reference taxa:\n";
				print join "\n", @available;
				if (@unavailable > 0){
					print "\n\n\tUnvailable reference taxa (no Blast db at $blastpath)\n";
					print join "\n", @unavailable;
				}
			}
		}
		else {
			print "\n${bold}NO CORE ORTHOLOG SETS ARE AVAILABLE! CHECK $hmmpath!${norm}\n\n";
		}
		print "\n\n";
		exit;
	}
	######### push all user defined variables into the log file ################
	push @log, "\nUSER DEFINED PARAMTERS (inc. default values)\n";
	my %parameters = (append => $append,
	blastpath => $blastpath,
	checkCoorthologsRef => $checkCoRef,
	cleartmp => $cleartmp,
	concat => $concat,
	est => $estflag,
	eval_blast => $eval_blast,
	eval_hmmer => $eval,
	filter => $filter,
	hit_limit => $hitlimit,
	hmm => $hmm,
	hmmset => $hmmset,
	hmmpath => $hmmpath,
	intron => $keepintron,
	longhead => $longhead,
	nonoverlapping_cos => $nonoverlappingCO,
	outpath => $outpath,
	protein => $proteinflag,
	rbh => $bhh,
	refspec => $refspec_string,
	relaxed => $relaxed,
	representative => $rep,
	reuse => $reuse,
	sequence_file => $dbfile,
	show_hmmsets => $show_coreortholog_sets,
	sort_global_align => $sortalign,
	strict => $strict,
	taxon => $taxon_global,
	ublast => $ublast);

	foreach ( sort keys %parameters) {
		if (defined $parameters{$_}) {
			push @log, "\t -$_:\t$parameters{$_}";
		}
		else {
			push @log, "\t -$_:\tnot set";
		}
	}

	#############################################################################
	my $check = 1;

	if (!defined $dbfile) {
		push @log, "You need to specify a valid infile with the option -sequence_file\n\n";
		$check = 0;
		return($check, @log);
	}
	### for FACT use the unmodified value of $dbfile
	$runFACTparameter = $dbfile;
	## extract the path from the dbpath if available and prune of trailing '/'
	if ($dbfile =~ /(.*\/)/) {
		$dbpath = $1;
		$dbpath =~ s/\/$//;
	}
	else {
		$dbpath = '.';

	}
	$dbfile =~ s/.*\///;
	# $dbfile_short = $dbfile;
	# $dbfile_short =~ s/\..*//;
	my @dbfileTMP = split(/\./, $dbfile); pop @dbfileTMP;
	$dbfile_short = join(".", @dbfileTMP);
	if ($central) {
		$dboutpath = $dbpath;
		# print "setting dboutpath to $dboutpath";
	}
	##
	## 0) Check for presence of the file with the sequences that should be hamstered
	if (-e "$dbpath/$dbfile") {
		push  @log, "\t$dbfile ready";
	}
	else {
		#the provided infile does not exist:
		push @log, "${bold}FATAL:${norm} The specified infile $dbpath/$dbfile does not exist. PLEASE PROVIDE A VALID INFILE!\n";
		$check = 0;
		return ($check, @log);
	}
	## 1) check for filetype
	printOUT("checking for sequence type:\n");
	if (!defined $estflag and !defined $proteinflag) {
		push @log, "\nCHECKING SEQUENCE TYPE\n";
		push @log, "\tNo file sequence type was determined. HaMStR will guess whether EST or protein sequences are analyzed";
		my $seq = `head -n 2 $dboutpath/$dbfile |tail -n 1`;
		my $orilength = length($seq);
		$seq =~ s/[AGCTN]//ig;
		if (length($seq) / $orilength >0.1) {
			$proteinflag = 1;
			printOUT("Guessing sequence type: Protein\n");
			push @log, "\tMore than 10% of the first sequence in the file are non-AGCTN. Guessing sequence type: Protein";
		}
		else {
			$estflag = 1;
			printOUT("Guessing sequence type: DNA\n");
			push @log, "\tLess than 10% of the first sequence in the file are non-AGCTN. Guessing sequence type: DNA\n";
		}
		$check = 1;
	}
	if ($estflag and !$check_genewise) {
		push @log, "\n\nHaMStR has been configured with the flag --protein_only and will not accept DNA sequences as input. I am stopping tests here! If you really want to analyse DNA sequence data please reconfigure.\n";
		$check = 0;
		return ($check, @log);
	}
	## $dbfile_base hat den originalen file namen, egal ob est oder protein
	$dbfile_base = $dbfile;

	if ($ublast){
		if ($runublast){
			$blast_prog = 'usearch';
			$algorithm = 'ublast';
			$outputfmt = 'blasttable';
			$blastdbend = '.udb';
		}
		else {
			push @log, "\n\nHaMStR has been configured with the --noublast option. Either re-start without the -ublast flag or reconfigure\n";
			$check = 0;
			return($check, @log);
		}
	}
	if ($estflag) {
		$dbfile = "$dbfile.tc";
		$algorithm = 'blastx';
		if ($blast_prog eq 'blastp'){
			$blast_prog = 'blastx';
		}
		$sortalign = 1;
		push @log, "HaMStR will run on the ESTs in $dbfile_base";
		push @log, "\nTRANSLATING ESTs\n";
		if (!(-e "$dboutpath/$dbfile")) {
			printOUT("translating $dbfile_base, this may take a while\n");
			`$path/bin/translate.pl -infile=$dboutpath/$dbfile_base -outfile=$dbfile -outpath=$dboutpath`;
			open (LOG, "$outpath/hamstrsearch.log");
			my @info = <LOG>;
			@log = (@log, @info);
			close LOG;
		}
		else {
			push @log, "Translated file already exists, using this one";
		}
		if (! -e "$dboutpath/$dbfile") {
			push @log, "${bold}FATAL:${norm} The translation of $dbfile_base failed. Check the script translate.pl";
			print "failed\n";
			$check = 0;
		}
		else {
			## file type is protein
			printOUT("succeeded\n");
		}
	}
	## 2) Check for presence of the blast program
	push @log, "\nCHECKING FOR PROGRAMS\n";
	printOUT("checking for the blast program:\t");
	if (`which $blast_prog` =~ / no /) {
		push @log, "${bold}FATAL:${norm} could not execute $blast_prog. Please check if this program is installed and executable";
		print "failed\n";
		$check = 0;
	}
	else {
		push @log, "\tcheck for $blast_prog succeeded";
		unless ($silent) {
			print "succeeded\n";
		}
	}
	## 3) Check for presence of hmmsearch
	printOUT("checking for hmmsearch:\t");
	my $hmmcheck = `$prog -h |$grepprog -c 'HMMER 3'`;
	if (! `$prog -h`) {
		push @log, "${bold}FATAL:${norm} could not execute $prog. Please check if this program is installed and executable";
		print "failed: $prog is not installed or not executable\n";
		$check = 0;
	}
	elsif ($hmmcheck != 1) {
		push @log, "${bold}FATAL:${norm} It seems that $prog is not from the HMMER 3 package. Please check!";
		print "failed: $prog is not from the HMMER 3 package\n";
		$check = 0;
	}
	else {
		push @log, "\tcheck for $prog succeeded";
		printOUT("succeeded\n");
	}
	## 3b) Check for genewise
	if ($check_genewise) {
		printOUT("checking for genewise:\t");
		if (! `genewise -help`) {
			push @log, "${bold}FATAL:${norm} Could not execute genewise. Please check if this program is installed and executable";
			print "failed: genewise is not executable\n";
			$check = 0;
		}
		else {
			my $gwcheck = `echo \$WISECONFIGDIR`;
			if (length($gwcheck) < 1) {
				push @log, "${bold}FATAL:${norm} The environmental variable WISECONFIGDIR has not been set. I am expecting troubles when invoking genewise.
				Please consult the installation manual for genewise and set this variable";
				print "failed: the environmental variable WISECONFIGDIR has not been set.\n";
				$check = 0;
			}
			else {
				printOUT("\tsucceeded\n");
			}
		}
	}
	else {
		push @log, "${bold}GENEWISE-CHECK skipped:${norm} The hamstr-script has been configured with the option --protein_only. To override this setting set reconfigure the script or set the variable $check_genewise to 1";
	}
	## 4) Check for presence of the directory structure

	push @log, "\nCHECKING FOR HMMs\n";
	printOUT("checking for presence of the hmm files:\t");
	if ( ! defined $hmmset or ! -e "$hmmpath/$hmmset") {
		push @log, "${bold}FATAL:${norm} You need to specify a valid core ortholog set. Make also sure that you provide the path to this set if it is not in the default location $hmmpath. You can check available core ortholog sets using the option -show_hmmsets.";
		print "failed\n";
		$check = 0;
	}
	else {
		$hmmpath = "$hmmpath/$hmmset";
		$fafile = "$hmmpath/$hmmset" . '.fa';
		$hmm_dir = "$hmmpath/$hmm_dir";
		$hmmsearch_dir = $outpath .'/hmm_search_' . $dbfile_short . '_' . $hmmset;

		## 4b) check for the presence of the hmm-files and the fasta-file
		if (!(-e "$hmm_dir")) {
			push @log, "${bold}FATAL:${norm} Could not find $hmm_dir";
			print "failed\n";
			$check = 0;
		} else {
			if (defined $hmm) {
				@hmms = split ',', $hmm;
				chomp @hmms;
				### check for the presence of all hmms
				for (my $k = 0; $k < @hmms; $k++) {
					if (! -e "$hmm_dir/$hmms[$k]") {
						push @log, "${bold}FATAL:${norm} $hmms[$k] has been defined but could not be found in $hmm_dir/$hmms[$k]";
						$check = 0;
						last;
					} else {
						push @log, "\t$hmms[$k] has been found";
					}
				}
			} else {
				push @log, "\trunning HaMStR with all hmms in $hmm_dir";
				my $hmm_dir_tmp = $hmm_dir; $hmm_dir_tmp =~ s/\|/\\\|/g;
				@hmms = `ls $hmm_dir_tmp`;
			}
			chomp @hmms;
			printOUT("\tsucceeded\n");
		}
	}
	## 6) Test for presence of the fasta file containing the sequences of the core-ortholog cluster
	printOUT("checking for presence of the core-ortholog file:\t");
	if (defined $fafile) {
		if (! -e "$fafile") {
			push @log, "Fatal: Could not find the file $fafile";
			print "failed\n";
			$check = 0;
		}
		else {
			push @log, "\tcheck for $fafile succeeded";
			printOUT("\tsucceeded\n");
		}
	}
	else {
		push @log, "${bold}FATAL:${norm} Please provide path and name of fasta file containing the core-ortholog sequences";
		$check = 0;
		print "failed\n";
	}
	## 7) Checks for the taxon_file
	push @log, "\nCHECKING TAXON NAME\n";
	printOUT("testing whether the taxon has been determined:\t");
	if (defined $taxon_global) {
		push @log, "\tusing default taxon $taxon_global for all sequences";
		printOUT("succeeded\n");
		$taxon_check = 2;
	}
	else {
		push @log, "${bold}FATAL:${norm} No taxon_file found. Please provide a global taxon name using the option -taxon";
		print "failed\n";
		$check = 0;
	}
	## 8) Check for reference taxon
	push @log, "\nCHECKING FOR REFERENCE TAXON\n";
	printOUT("checking for reference species and blast-dbs:\t");
	if (!(defined $refspec_string) and (! defined $strict and ! defined $relaxed)) {
		push @log, "${bold}FATAL:${norm} Please provide a reference species for the reblast!";
		print "failed\n";
		$check = 0;
	}
	elsif (defined $strict or defined $relaxed) {
		if (! defined $refspec_string) {
			## The user has not provided a string of reference taxa. Chose all from the fasta file containing
			## the core orthologs.
			@refspec = `$grepprog '>'  $fafile |cut -d '|' -f 2 |sort |uniq`;
			chomp @refspec;
			$refspec_string = join ',', @refspec;
		}
		else {
			@refspec = split (/,/, $refspec_string);
		}
		if ($strict) {
			push @log, "\tStrict flag has been set. Reference species for the reblast: All of $refspec_string";
		}
		else {
			push @log, "\tRelaxed flag has been set. Reference species for the reblast: Any of $refspec_string";
		}
		if (@refspec == 0) {
			print "failed\n";
			$check = 0;
		}
		else {
			printOUT("succeeded\n");
		}
	}
	else {
		push @log, "\t Reference species for the re-blast: $refspec_string";
		@refspec = split(/,/, $refspec_string);
		$refspec_name = $refspec[0];
		printOUT("succeeded\n");
	}
	## 9) Check for presence of the required blast dbs
	printOUT("checking for blast-dbs:\t");
	push @log, "\nCHECKING FOR BLAST DATABASES\n";
	for (my $i = 0; $i < @refspec; $i++) {
		my $blastpathtmp = "$blastpath/$refspec[$i]/$refspec[$i]";
		if (-e $blastpathtmp . $blastdbend) {
			push @log, "\tcheck for $blastpathtmp succeeded";
			printOUT("succeeded\n");
		}
		elsif (-e $blastpathtmp . '_prot' . $blastdbend){
			## the check for the file naming '_prot' is only to maintain backward compatibility
			$blastapp = '_prot';
			$blastpathtmp = $blastpathtmp . $blastapp;
			push @log, "\tcheck for $blastpathtmp succeeded";
			printOUT("succeeded\n");
		}
		else {
			push @log, "${bold}FATAL:${norm} please edit the blastpath. Could not find $blastpathtmp or blast database blastpathtmp.pin does not exist.";
			print "$blastpathtmp failed\n";
			$check = 0;
		}
	}
	## 9.1) Check for presence of the required FASTA file of reference species
	printOUT("checking for reference fasta files:\t");
	push @log, "\nCHECKING FOR REFERENCE FASTA FILES\n";
	for (my $i = 0; $i < @refspec; $i++) {
		my $referencedb = "$blastpath/$refspec[$i]/$refspec[$i]".".fa";
		my $referencedb_prot = "$blastpath/$refspec[$i]/$refspec[$i]"."_prot.fa"; # backward compatibility
		my $ref_dir = "$blastpath/$refspec[$i]";
		my $link = $referencedb;
		unless (-e $referencedb) {
			$link = `$readlinkprog $referencedb`;
			unless ($link =~ /^\./ || $link =~ /^\//) {
				my $cwd = cwd();
				die "Linked source for $referencedb not found in $cwd!";
			}
		}
		# my $ref_location = $referencedb; # not used anywhere else
		chomp($link);
		if (-e $referencedb || -e $link) {
			push  @log, "\tinfile ready";
		} elsif (-e "$referencedb_prot"){
			push  @log, "\tinfile ready";
		} else {
			#the provided reference fasta file does not exist or link to file does not exist:
			push @log, "${bold}FATAL:${norm} FASTA file for the specified reference $refspec[$i] does not exist. PLEASE PROVIDE A VALID REFERENCE SPECIES!\n";
			$check = 0;
			return ($check, @log);
		}
	}

	## 10) Set the file where the matched seqs are found
	my $strictstring = '';
	if (defined $strict) {
		$strictstring = '.strict';
	}
	$seqs2store_file = $outpath . '/hamstrsearch_' . $dbfile_short . '_' . $hmmset . $strictstring . '.out';
	$cds2store_file = $outpath . '/hamstrsearch_' . $dbfile_short . '_' . $hmmset . '_cds' . $strictstring . '.out';

	if (! $append){
		if (-e "$seqs2store_file") {
			my $answer = 'Y';
			my $breaker = 0;
			if (!$force){
				print "A HaMStR outfile $seqs2store_file already exists and option -force has not been chosen! Shall I overwrite this file [Y|N]: ";
				$answer = <STDIN>;
				chomp $answer;
				while ($answer !~ /[YN]/i and ($breaker < 4)) {
					$breaker ++;
					print "Please answer with 'Y' or 'N':\t";
					$answer = <STDIN>;
					chomp $answer;
					if (($breaker > 3) and ($answer !~ /[YN]/i)){
						print "No proper answer is given: exiting.\nPlease re-start HaMStR with the -append option, or alternatively remove the file manually, or force the replacement of exsiting files with option -force.\n";
						exit;
					}
				}
			}
			if ($answer =~ /Y/i) {
				open (OUT, ">$seqs2store_file") or die "failed to open $seqs2store_file\n";
				print OUT '';
				close OUT;
				if ($estflag){
					open (OUT, ">$cds2store_file") or die "failed to open $cds2store_file\n";
					print OUT '';
					close OUT;
				}
			}
			else {
				print "You chose to not overwrite the existing output files. Please re-start HaMStR with the -append option, or alternatively remove the file manually.\n";
				exit;
			}
		}
	}
	## 11) apply the evalue-cut-off to the hmmsearch program
	push @log, "\nPROGRAM OPTIONS\n";
	push @log, "\thmmsearch will run with an e-value limit of $eval";

	## 11b) hit limit for the re-blast
	if ($hitlimit) {
		push @log, "\tre-blast hit_limit: $hitlimit";
	}
	else {
		push @log, "\tre-blast hit_limit: none applied";
	}
	## 11c) The blast evalue limit
	push @log, "\tBlast will run with an evalue limit of $eval_blast\n";

	## 12) check for filter setting for BLAST
	printOUT("checking for low complexity filter setting:\t");
	$filter =~ tr/ft/FT/;
	if ($filter ne 'T' and $filter ne 'F') {
		push @log, "${bold}FATAL:${norm} Filter is set to $filter. Please set the low complexity filter either to F or T.";
		print "low complexity filter check failed\n";
		$check = 0;
	}
	else {
		push @log, "\tcheck for low complexity filter setting succeeded. Chosen value is $filter";
		if ($blast_prog ne 'blastall'){
			$filter = 'yes' if $filter eq 'T';
			$filter = 'no' if $filter eq 'F';
		}
		printOUT("succeeded\n");
	}

	## 13) setting up the directories where the output files will be put into.
	$fa_dir_neu = $outpath . '/fa_dir_' . $dbfile_short . '_' . $hmmset . '_' . $refspec[0];
	$tmpdir = $outpath . '/tmp'; # . $tmpdir;
	if (!$strict) {
		$fa_dir_neu = $outpath . '/fa_dir_' . $dbfile_short . '_' . $hmmset . '_' . $refspec[0];
	}
	if ($strict) {
		$fa_dir_neu = $outpath . '/fa_dir_' . $dbfile_short . '_' . $hmmset;
		$fa_dir_neu .= '_strict';
	}

	if ($relaxed) {
		$fa_dir_neu .= '_relaxed';
	}
	if ($check == 1) {
		if (!(-e "$hmmsearch_dir")) {
			`mkdir "$hmmsearch_dir"`;
		}
		elsif (-e "$hmmsearch_dir" and ! $reuse) {
			`rm -rf "$hmmsearch_dir"`;
			`mkdir "$hmmsearch_dir"`;
		}
		if (!(-e "$fa_dir_neu")) {
			`mkdir "$fa_dir_neu"`;
		}
		elsif (-e "$fa_dir_neu" and ! $reuse) {
			`rm -rf "$fa_dir_neu"`;
			`mkdir "$fa_dir_neu"`;
		}
		if (!(-d "$tmpdir")) {
			`mkdir "$tmpdir"`;
		}
		elsif (-d "$tmpdir" and $cleartmp) {
			`rm -rf "$tmpdir"`;
			`mkdir "$tmpdir"`;
		}
	}
	## 14) determin whether or not the -representative flag has been set
	if (defined $rep) {
		push @log, "\tHaMStR will run with the -representative option";
	}
	else {
		push @log, "\tHaMStR was called without the -representative option. More than one ortholog may be identified per core-ortholog group!";
	}

	## check further options
	if (defined $nonoverlappingCO){
		push @log, "\tThe flag -nonoverlapping_cos has been set. HaMStR will output co-orthologs even when they align to non-overlapping parts of the reference sequence";
	}
	if (defined $checkCoRef){
		push @log, "\tThe flag -CheckCoorthologsRef has been set.";
	}
	if (defined $bhh){
		push @log, "\tThe flag -rbh has been set. HaMStR will run with the reciprocal best hit option.";
	}
	if ($sortalign){
		push @log, "\tThe flag -sort_global_align has been set. HaMStR will sort hits according to the global alignment score against the reference sequence. (Default for EST data)."
	}

	## check how hamstr should deal with possible introns in transcripts:
	if ($estflag) {
		my $breaker = 0;
		while ($keepintron !~ /^[kmr]/i and ($breaker < 4)){
			$breaker ++;
			print "option intron was set to $keepintron: Please answer either with 'k(eep)', 'm(ask)', or 'r(emove)':\t";
			$keepintron = <STDIN>;
			chomp $keepintron;
			if (($breaker > 3) and ($keepintron !~ /^[kmr]/i)){
				print "No proper answer is given: exiting.\nPlease re-start HaMStR with the option -intron=[kmr].\nOptions are 'k(eep)', 'm(ask)', or 'r(emove)'. Default is 'k(eep)' introns.\n";
				exit;
			}
		}
		if ($keepintron =~ /^k/i) {
			push @log, "\tKeep introns (Default) has been chosen. HaMStR will keep any introns in lower case in the reported CDS. Thus, CDS cannot be directly translated into the aa sequence.";
		}
		elsif ($keepintron =~ /^m/i) {
			push @log, "\tMask introns has been chosen. HaMStR will keep any introns but masks them as 'N' in the reported CDS. Thus, CDS cannot be directly translated into the aa sequence."
		}
		elsif ($keepintron =~ /^r/i) {
			push @log, "\tRemove introns has been chosen. HaMStR will remove any position that genewise could not align to the reference protein rendering the CDS consistent with the amino acid sequence";
		}

	}

	return ($check, @log);
}
#################
## check4reciprocity is the second major part of the program. It checks
## whether the protein sequence that has been identified by the hmmsearch
## identifies in turn the protein from the reference taxon that was used to
## build the hmm.
sub check4reciprocity {
	my $frame;
	my ($localid, $query_name, $hitname, $refspec_final, @refspec) = @_;
	my $searchdb;
	my $strict_suc = -1; # keeps track of success for all taxa
	my $relaxed_suc = 0; # keeps track of success for at least one taxon
	## get the sequence that was identified as hit in the pHMM search from the db_file
	my $hitseq = `$grepprog -m 1 -A 1 ">$hitname\$" $dboutpath/$dbfile_base | tail -n 1`;
	if (!defined $hitseq) {
		print "could not retrieve a sequence for $hitname. Skipping...\n";
		return(0, '', '', '');
	}
	## continue with the blast
	chomp $hitseq;
	## now run the blast
	open (OUT, ">$tmpdir/$$.fa") or die "could not open out for writing\n";
	print OUT ">$hitname\n$hitseq";
	close OUT;
	## now comes the new part that does one to many blast searches. We need to iterate through all
	## entries in the file $refspec_final and perform the Blast against each reftaxon. Note, unless
	## $strict or $relaxed flags are set, there will be only a single reftaxon. If $relaxed is chosen
	## then we can stop the blast searches as soon as the reciprocity is fulfilled.
	for (my $k = 0; $k < @$refspec_final; $k++) {
		my $orthocount = $refspec_final->[$k]->{orthocount};
		## 1) Perform the blast search with the k-th reftaxon
		printOUT("Reftaxon: $refspec_final->[$k]->{refspec}\n");
		$tmpdir =~ s/\|/\\\|/g;
		if ($blast_prog =~ /blast[px]/) {
			!`$blast_prog -db $refspec_final->[$k]->{searchdb} -seg '$filter' -max_target_seqs 10 -evalue $eval_blast -outfmt 5 -query $tmpdir/$$.fa  -out $tmpdir/$$.blast` or die "Problem running $blast_prog\n";
			### postprocess the outfile
		}
		elsif ($blast_prog =~ /blastall/) {
			!`blastall -p $algorithm -d $refspec_final->[$k]->{searchdb} -F $filter -e $eval_blast -m7 -i $tmpdir/$$.fa -o $tmpdir/$$.blast` or die "Problem running $blast_prog\n"
		}
		else {
			if ($estflag){
				`$blast_prog -ublast $tmpdir/$$.fa -db $refspec_final->[$k]->{searchdb}.udb -strand both -accel $accel -evalue $eval_blast -blast6out $tmpdir/$$.blast` or die "Problem running $blast_prog\n;"
			}
			else {
				`$blast_prog -ublast $tmpdir/$$.fa -db $refspec_final->[$k]->{searchdb}.udb -accel $accel -evalue $eval_blast -blast6out $tmpdir/$$.blast` or die "Problem running $blast_prog\n;"
			}
			## sort the output as ublast does not do it (at least not for ESTs)
			`sort -n -r -k 12 $tmpdir/$$.blast >$tmpdir/blastsort.tmp`;
			`mv $tmpdir/blastsort.tmp $tmpdir/$$.blast`;
			####################
		}
		## 2) now parse the best blast hit
		my $hits = &getBestBlasthit("$tmpdir/$$.blast");
		if (defined $hits and @$hits > 0) {
			## at least one blast hit
			$frame = $hits->[0]->{frame};
			my $idsref = $refspec_final->[$k]->{refid};
			my @original_ids = @$idsref;
			my $suc = 0; # keeps track of success for a single taxon
			if ($checkCoRef == 0) {
				## the user does not want to check further in case that id of best blast hit and of reference species differ
				printOUT("core_orthologs: ", join "\t", @original_ids , "\n");
				## now loop through the best hits with the same score and check whether
				## among these I find the same seq as in $original
				my $i = 0;
				while ($suc == 0 and $i <@$hits) {
					printOUT("blast-hit: $hits->[$i]->{name}");
					## now loop through all the refspec-sequences in the hmm file; this is the case when co-orthologs have been determine in the core-ortholog
					my $j = 0;
					while ($suc == 0 and $j < @original_ids) {
						if ($original_ids[$j] eq $hits->[$i]->{name}) {
							printOUT("\thitting\n");
							$refspec_final->[$k]->{hit} = $j;
							$suc = 1;
							$relaxed_suc = 1;
						}
						else {
							printOUT("\nnot hitting $original_ids[$j]\n");
							$j ++;
						}
						if ($suc == 1) {
							$relaxed_suc = 1;
							if ($strict_suc == -1) {
								$strict_suc = 1;
							}
						}
					}
					$i++;
				}
				if ($suc == 0) {
					# none of the blast hits matched against the the reftaxon seq
					$strict_suc = 0;
				}
			}

			else {
				## The user has chosen to search more sensitive, asking whether the best blast hit might be a co-ortholog to the reference sequence
				my $qhdistance;
				my $rhdistance;
				printOUT("core_orthologs: $original_ids[0]\n");
				## we will check only the best blast hit and impose a distance criterion
				## in case of an EST, we will have to predict the reading frame and translate it...
				my $bestid = $hits->[0]->{name};
				my $refid = $original_ids[0];
				## get the sequences from the blast db. Currently, I'm using a simple grep
				my $bestseq = `$grepprog -m 1 -A 1 ">$bestid" $refspec_final->[$k]->{searchdb}.fa |tail -n 1` or die "Could not retrieve original sequence for besthit\n";
				my $refseq = `$grepprog -m 1 -A 1 ">$refid" $refspec_final->[$k]->{searchdb}.fa |tail -n 1` or die "Could not retrieve original sequence for refseq\n";
				chomp ($bestseq, $refseq);
				printOUT("blast-hit: $bestid");
				my $queryseq = $hitseq;
				if ($bestid eq $refid) {
					printOUT("\thitting\n");
					$refspec_final->[$k]->{hit} = 0;
					$suc = 1;
					$relaxed_suc = 1;
				}
				else {
					printOUT("\nBest hit $bestid differs from reference sequence $refid! Doing further checks\n");
					if ($estflag){
						printOUT("Frame is $hits->[0]->{frame} or $frame\n");
						my ($hitseqtr) = &findORF($hitseq, $bestseq, $frame);
						($suc, $qhdistance, $rhdistance) = &checkCoorthologRef($localid, $hitseqtr, $bestseq, $refseq);
					}
					else {
						($suc, $qhdistance, $rhdistance) = &checkCoorthologRef($localid, $hitseq, $bestseq, $refseq);
					}
					## print distances (debug mode)
					if ($debug){
						my $distDebugFile = $path . "/output/" . $taxon_global . ".debug.dist";
						unless (-e $distDebugFile){
							open (my $DISTDEBUG, ">>$distDebugFile") or die "Error, could not create file: ". "$distDebugFile";
							print $DISTDEBUG "hmmset\trefid\tbestid\tqueryid\tqhdist\trhdist\n";
							close $DISTDEBUG;
						}
						if (-e $distDebugFile){
							open (my $DISTDEBUG, ">>$distDebugFile") or die "Error, could not create file: ". "$distDebugFile";
							print $DISTDEBUG "$query_name\t$refid\t$bestid\t$hitname\t$qhdistance\t$rhdistance\n";
							close $DISTDEBUG;
						}
					}

					if ($suc == 1) {
						printOUT("\t Distance query - blast hit: $qhdistance, Distance blast hit - reference: $rhdistance\tAccepting\n");
						$refspec_final->[$k]->{hit} = 0;
					}
					else {
						printOUT("\t Distance query - blast hit: $qhdistance; Distance blast hit - reference: $rhdistance Rejecting\n");
					}
				}
				if ($suc == 1){
					$relaxed_suc = 1;
					if ($strict_suc == -1) {
						$strict_suc = 1;
					}
				}
				else {
					$strict_suc = 0;
				}
			}
		}
		else {
			printOUT("no hit obtained\n");
			$strict_suc = 0;
		}
		## when the user has chosen the strict flag, there is no reason to continue when $suc
		## has remained 0 (reciprocity criterion not fulfilled). Thus, return to main.
		if ($strict and $strict_suc == 0) {
			return (0, $hitseq);
		}
	}
	if ($relaxed_suc == 1) {
		if ($estflag and $frame eq '-') {
			## reverse sequence
			$hitseq = &revComp($hitseq);
		}
		return (1, $hitseq, $frame);
	}
	else {
		return (0, $hitseq);
	}
}

#############
sub getBestBlasthit {
	my $hits;
	my $count = 0;
	my ($file) = @_;
	$file =~ s/\\//g;
	my $searchio = Bio::SearchIO->new(
	-file        => "$file",
	-format      => $outputfmt)
	or die "parse failed";
	while(my $result = $searchio->next_result){
		my $sig;
		my $sig_old;
		while( my $hit = $result->next_hit) {
			my $frameval = $hit->strand('query');
			if ($frameval >0){
				$frame = '+';
			}
			elsif ($frameval <0 ) {
				$frame = '-';
			}
			elsif (!defined $frameval and $estflag) {
				die "error in obtaining frame in sub getBestBlasthit\n";
			}
			else {
				$frame = 'na';
			}

			## now I enter all top hits having the same evalue into the result
			$sig = $hit->score;
			if (!defined $sig_old) {
				$sig_old = $sig;
			}
			if ($sig == $sig_old) {
				if ($estflag){
					printOUT("frame is $frame\n");
					$hits->[$count]->{frame} = $frame;
				}
				$hits->[$count]->{name} = $hit->name;
				$count ++;
			}
			else {
				## there is no lower ranking hit with the same score as the best hit. End the loop.
				last;
			}
		}
	}
	return($hits);
}
##################
sub getTaxon {
	my ($hitname) = @_;
	if ($hitname =~ /\D/) {
		$hitname =~ s/_.*//;
	}
	my $taxon = `$grepprog -m 1 "^$hitname," $taxon_file | $sedprog -e 's/^.*,//'`;
	chomp $taxon;
	$taxon =~ s/^[0-9]+,//;
	$taxon =~ s/\s*$//;
	$taxon =~ s/\s/_/g;
	if ($taxon) {
		return ($taxon);
	}
	else {
		return();
	}
}
###############
sub determineReferences {
	my ($localid, $fileobj, $taxon, $refspec_final, $hitname, $hithmmscore, $hitseq, $hitcount) = @_;
	my $refseq = '';
	my $refspec;
	## now we have to distinguish between three cases:
	## 1) hamstr is running in normal mode and one refspec has been determined. In this case, $refspec_final
	## contains data only from a single species.
	## 2) hamstr is running in normal mode and alternative refspecs have been determined by the user.
	## $refspec_final may contain results from more than one species, but we need to consider only the first
	## entry.
	## 3) hamstr is running in the strict mode. In this case $refspec_final contains data from several taxa and we need
	## to select the taxon and sequence that is most similar to the hamstered sequence.
	## 4) hamstr is running in the relaxed mode. In this case $refspec_final may contain data from several taxa and
	## we need to select the taxon and the sequence that is most similar to the hamstered sequence.
	if (defined $strict or defined $relaxed) {
		## more than one refspec. Now find the one that fits best
		my $max_score = 0;
		for (my $i = 0; $i < @$refspec_final; $i++) {
			## first, check whether the reciprocity criterion has been fulfilled
			if (defined $refspec_final->[$i]->{hit}) {
				my $rcn = $refspec_final->[$i]->{hit};
				my $refseq_cand = $refspec_final->[$i]->{sequence}->[$rcn];
				my $refspec_cand_id = $refspec_final->[$i]->{refid}->[$rcn];
				my $refspec_cand = $refspec_final->[$i]->{refspec};
				my $score = &getAlignmentScore($localid, $refseq_cand, $hitseq);
				if ($score > $max_score) {
					$refspec = $refspec_cand;
					$refseq = $refseq_cand;
					$max_score = $score;
				}
			}
		}
	}
	else { ## no choice, just one refspec
	my $rcn = $refspec_final->[0]->{hit};
	$refseq = $refspec_final->[0]->{sequence}->[$rcn];
	$refspec = $refspec_final->[0]->{refspec};
}
$fileobj->{$taxon}->{prot}->[$hitcount] = $hitseq;
$fileobj->{$taxon}->{ids}->[$hitcount] = $hitname;
$fileobj->{$taxon}->{hmmscore}->[$hitcount] = $hithmmscore;
$fileobj->{$taxon}->{refseq}->[$hitcount]= $refseq;
$fileobj->{$taxon}->{refspec}->[$hitcount] = $refspec;
return($fileobj);
}
###############
sub processHits {
	my ($localid, $fileobj) = @_;
	## 1) align all hit sequences for a taxon against the reference species
	my @taxa = keys(%$fileobj);
	for (my $i = 0; $i < @taxa; $i++) {
		&orfRanking($localid, $taxa[$i]);
	}
}


################
sub predictORF {
	my $fileobj_new;
	my @taxa = keys(%$fileobj);
	for (my $i = 0; $i < @taxa; $i++) {
		my $protobj = $fileobj->{$taxa[$i]}->{prot};
		my $idobj = $fileobj->{$taxa[$i]}->{ids};
		my $refseqobj = $fileobj->{$taxa[$i]}->{refseq};
		my $refspecobj = $fileobj->{$taxa[$i]}->{refspec};
		my @ids = @$idobj;
		for (my $j = 0; $j < @ids; $j++) {
			my $refseq = $refseqobj->[$j];
			my $refspec = $refspecobj->[$j];
			my $est = $protobj->[$j];
			if (! $est) {
				die "error in retrieval of est sequence for $ids[$j] in subroutine processHits\n";
			}
			### debuggin IUB code
			if ($est =~ /[^AGCT]/i) {
				$est =~ s/[^AGCTagct]/n/g;
			}
			printOUT("running genewise using frame $frame\n");
			my $gw = run_genewise_hamstr->new($est, $refseq, $tmpdir, $keepintron);
			my $translation = $gw->translation;
			my $cds = $gw->codons;
			$translation =~ s/[-!]//g;
			$fileobj_new->{$taxa[$i]}->{ids}->[$j] = $ids[$j];
			$fileobj_new->{$taxa[$i]}->{prot}->[$j] = $translation;
			$fileobj_new->{$taxa[$i]}->{cds}->[$j] = $cds;
			$fileobj_new->{$taxa[$i]}->{refseq}->[$j] = $refseq;
			$fileobj_new->{$taxa[$i]}->{refspec}->[$j] = $refspec;
		}
	}
	return($fileobj_new);
}
############################
sub orfRanking {
	my ($localid, $spec) = @_;
	my $result;
	my $refprot;
	my $refcds;
	my @toalign;
	my $protobj = $fileobj->{$spec}->{prot};
	my $idobj = $fileobj->{$spec}->{ids};
	my $refcluster; ## variables to take the cluster and its id for later analysis
	my $refid;
	if (@$protobj == 1) {
		## nothing to chose from
		$refprot = $protobj->[0];
		$refcds = $fileobj->{$spec}->{cds}->[0];
		my $length = length($refprot);
		$refid = $idobj->[0];
	}
	else {
		## more than one cluster
		## note, I set the refseq fix to the first entry. This is to avoid that in this routine
		## sequences from different taxa are used.
		push @toalign, ">$fileobj->{$spec}->{refspec}->[0]";
		push @toalign, $fileobj->{$spec}->{refseq}->[0];
		## now walk through all the contigs
		for (my $i = 0; $i < @$protobj; $i++) {
			my @testseq = (">$idobj->[$i]", $protobj->[$i]);
			@testseq = (@testseq, @toalign);
			open (OUT, ">$tmpdir/$localid.ref.fa") or die "could not open file for writing refseqs\n";
			print OUT join "\n", @testseq;
			close OUT;
			## run clustalw
			!(`$alignmentprog -infile=$tmpdir/$localid.ref.fa -output=fasta -outfile=$tmpdir/$localid.ref.aln 2>&1 >$tmpdir/$localid.ref.log`) or die "error running clustalw\n";
			## get the alignment score
			$result->[$i]->{score} =  `$grepprog "Alignment Score" $tmpdir/$localid.ref.log |$sedprog -e 's/[^0-9]//g'`;
			if (!$result->[$i]->{score}) {
				die "error in determining alignment score\n";
			}
			chomp $result->[$i]->{score};
			## get the aligned sequence
			open (ALN, "$tmpdir/$localid.ref.aln") or die "failed to open alignment file\n";
			my @aln = <ALN>;
			close ALN;
			my $aseq = extractSeq($idobj->[$i], @aln);
			## remove the terminal gaps
			$aseq =~ s/-*$//;
			$result->[$i]->{aend} = length $aseq;
			my ($head) = $aseq =~ /^(-*).*/;
			($result->[$i]->{astart}) = length($head)+1;
			## add the hmmscore to $result
			$result->[$i]->{hmmscore} = $fileobj->{$spec}->{hmmscore}->[$i];
		}
		### the results for all seqs has been gathered, now order them according to alignment start in the refseq
		$result = &sortRef($result);
		($refprot, $refcds, $refid) = &determineRef($result,$spec);
	}
	$fileobj->{$spec}->{refprot} = $refprot;
	$fileobj->{$spec}->{refcds}  = $refcds;
	$fileobj->{$spec}->{refid}   = $refid;
	$fileobj->{$spec}->{refspec_final} = $fileobj->{$spec}->{refspec}->[0];
	return();
}
###########################
sub sortRef {
	my $result = shift;
	my @sortref;
	for (my $i = 0; $i < @$result; $i++) {
		$sortref[$i]->{index} = $i;
		$sortref[$i]->{astart} = $result->[$i]->{astart};
		$sortref[$i]->{aend} = $result->[$i]->{aend};
		$sortref[$i]->{score} = $result->[$i]->{score};
		$sortref[$i]->{hmmscore} = $result->[$i]->{hmmscore};
	}
	@sortref = sort { $a->{astart} <=> $b->{astart} } @sortref;
	for (my $i = 0; $i < @sortref; $i++) {
		($result->[$i]->{id}, $result->[$i]->{start}, $result->[$i]->{end}, $result->[$i]->{score}, $result->[$i]->{hmmscore}) = ($sortref[$i]->{index}, $sortref[$i]->{astart}, $sortref[$i]->{aend}, $sortref[$i]->{score}, $sortref[$i]->{hmmscore});
	}
	return($result);
}
########################
sub determineRef {
	my ($result, $spec) = @_;
	my $lastend = 0;
	my $lastscore = 0;
	my $final;
	my $count = 0;
	my $id = '';
	my $scorekey = 'hmmscore';
	if ($sortalign){
		$scorekey = 'score';
	}
	for (my $i = 0; $i < @$result; $i++) {
		if ($result->[$i]->{start} < $lastend or $lastend == 0) {
			if ($result->[$i]->{$scorekey} > $lastscore) {
				$lastend = $result->[$i]->{end};
				$lastscore = $result->[$i]->{$scorekey};
				$id = $result->[$i]->{id};
				printOUT("ref is $id with score $lastscore\n");
			}
		}
		elsif ($result->[$i]->{start} > $lastend) {
			## a new part of the alignment is covered. Fix the results obtained so far
			$final->[$count]->{id} = $id;
			$lastend = $result->[$i]->{end};
			$id = $result->[$i]->{id};
			$count++;
		}
	}
	$final->[$count]->{id} = $id;
	## now concatenate the results
	my $refprot = '';
	my $refid = '';
	my $refcds = '';

	## now comes a dirty hack. The user has the chance to maximize phylogentic information by concatenating
	## orthologous sequences that do no align to the same part of the reference protein (option -concat). If so,
	## the co-ortholog-detection at a later step will not work and will be disabled.
	my $looplimit = 1;
	if ($concat) {
		$looplimit = scalar(@$final);
	}
	for (my $i = 0; $i < $looplimit; $i++) {
		my $seq = $fileobj->{$spec}->{prot}->[$final->[$i]->{id}];
		my $cdsseq = $fileobj->{$spec}->{cds}->[$final->[$i]->{id}];
		my $length = length($seq);
		if ($concat){
			$refid .= "$fileobj->{$spec}->{ids}->[$final->[$i]->{id}]-$length" . "PP";
		}
		else {
			$refid .= "$fileobj->{$spec}->{ids}->[$final->[$i]->{id}]";
		}
		$refprot .= $seq;
		if ($estflag) {
			$refcds .= $cdsseq;
		}
	}
	$refid =~ s/PP$//;
	return($refprot, $refcds, $refid);
}
#############################
sub extractSeq {
	my ($id, @aln) = @_;
	my $seq = '';
	my $start = 0;
	for (my $i = 0; $i < @aln; $i++) {
		if ($aln[$i] =~ $id) {
			$start = 1;
		}
		elsif ($aln[$i] =~ />/ and $start == 1) {
			last;
		}
		elsif ($start == 1) {
			$seq .= $aln[$i];
		}
	}
	$seq =~ s/\s//g;
	return ($seq);
}
##############################
sub revComp {
	my ($seq) = @_;
	chomp($seq);
	$seq =~ tr/AGCTYRKMWSagct/TCGARYMKWSTCGA/;
	$seq = reverse($seq);
	return($seq);
}
##############################
sub parseHmmer3pm {
	my ($file, $path) = @_;
	my $hits;
	my $query;
	my %tmphash;
	if (!defined $path){
		$path = '.';
	}
	$file = $path . '/' . $file;
	my $in = Bio::SearchIO->new(
	-format => 'hmmer',
	-file   => $file
	);
	while( my $result = $in->next_result ) {
		# this is a Bio::Search::Result::HMMERResult object
		if (!defined $query){
			$query = $result->query_name();
			printOUT("query is $query\n");
		}
		my $hitcount = 0;
		while( my $hit = $result->next_hit ) {
			my $tmp = $hit->name();
			my $tmpscore = $hit->score();
			$tmp =~ s/_RF.*//;
			if (!defined $tmphash{$tmp}){
				$hits->[$hitcount]->{id} = $tmp;
				$hits->[$hitcount]->{hmmscore} = $tmpscore;
				$hitcount++;
				$tmphash{$tmp}=1;
				if (defined $bhh){
					last;
				}
			}
		}

		if (defined $hits->[0]) {
			####### a quick hack to obtain the lagPhase value
			my $criticalValue; # takes the value used for candidate discrimination
			my $hitLimitLoc = $hitlimit;
			if (defined $autoLimit) {
				printDebug("Entering getLag Routine\n");
				## the user has invoked the autmated inference of a hit limit
				($hitLimitLoc, $criticalValue)  = getLag($hits, $hitcount);
				if (!defined $criticalValue) {
					## there was a problem in the computatation of the lagPhase
					print "Computation of lagPhase did not succeed, switching to score threshold using a default cutoff of $scoreCutoff\n";
					($hitLimitLoc, $criticalValue) = getHitLimit($hits, $hitcount);
				}
			}
			elsif (defined $scoreThreshold) {
				printDebug("entering the scoreThreshold routine");
				($hitLimitLoc, $criticalValue) = getHitLimit($hits, $hitcount);
				printDebug("hitlimitloc is now $hitLimitLoc");
			}

			return ($query, $hits, $hitLimitLoc, $criticalValue);
		}
		else {
			return ($query);
		}
	}
}
##############################
sub parseHmmer4pm {
	my ($file, $path) = @_;
	my $hmmhits;
	my $hits;
	my $query;
	my @rest;
	my %tmphash;
	my $hitcount = 0;
	if (!defined $path){
		$path = '.';
	}
	$file = $path . '/' . $file;

	$file =~ s/\|/\\\|/g;
	my @hmmout = `$grepprog -v '#' $file |sort -rnk 9 |sed -e 's/ /@/g'`;
	for (my $i = 0; $i < @hmmout; $i++) {
		($hmmhits->[$i]->{target_name}, $hmmhits->[$i]->{target_accession}, $hmmhits->[$i]->{query_name}, $hmmhits->[$i]->{query_accession},  $hmmhits->[$i]->{total_evalue},  $hmmhits->[$i]->{total_score},  $hmmhits->[$i]->{total_bias},  $hmmhits->[$i]->{domain_evalue}, $hmmhits->[$i]->{domain_score},  $hmmhits->[$i]->{domain_bias}, @rest) = split(/@+/, $hmmout[$i]);

		if (!defined $query){
			$query = $hmmhits->[$i]->{query_name};
			printOUT("query is $query\n");
		}
		my $tmp = $hmmhits->[$i]->{target_name};
		my $tmpscore = $hmmhits->[$i]->{domain_score};
		$tmp =~ s/_RF.*//;
		if (!defined $tmphash{$tmp}){
			$hits->[$hitcount]->{id} = $tmp;
			$hits->[$hitcount]->{hmmscore} = $tmpscore;
			$hitcount++;
			$tmphash{$tmp}=1;
			if (defined $bhh){
				last;
			}
		}

	}
	if (defined $hits->[0]) {
		####### limit the list of hmm hits
		my $criticalValue; # takes the value used for candidate discrimination
		my $hitLimitLoc = $hitlimit;
		if (defined $scoreThreshold) {
			printDebug("entering the scoreThreshold routine");
			($hitLimitLoc, $criticalValue) = getHitLimit($hits, $hitcount);
			printDebug("hitlimitloc is now $hitLimitLoc");
		}

		return ($query, $hits, $hitLimitLoc, $criticalValue);
	}
	else {
		return ($query);
	}

}
##############################
sub parseSeqfile {
	my $seqref;
	my $id;
	my $spec;
	my $seqid;
	my $seq;
	my $file = shift;
	open (IN, "$file") or die "failed to open $file\n";
	my @seqs = <IN>;
	close IN;
	chomp @seqs;
	for (my $i = 0; $i < @seqs; $i++) {
		if ($seqs[$i] =~ />/) {
			$seqs[$i] =~ s/>//;
			if (defined $id and defined $seq) {
				$seqref->{$id}->{$spec}->{seqid} = $seqid;
				$seqref->{$id}->{$spec}->{seq} = $seq;
				$seq = undef;
			}
			($id, $spec, $seqid) = split (/\|/, $seqs[$i]);
		}
		else {
			$seq .= $seqs[$i];
		}
	}
	if (defined  $id and defined $seq) {
		$seqref->{$id}->{$spec}->{seqid} = $seqid;
		$seqref->{$id}->{$spec}->{seq} = $seq;
		$seq = undef;
	}
	return ($seqref);
}
##################
sub getAlignmentScore{
	my ($localid, $refseq_cand, $hitseq) = @_;
	my @testseq = ('>hitseq', $hitseq, '>refseq', $refseq_cand);
	open (OUT, ">$tmpdir/$localid.ref.fa") or die "could not open file for writing refseqs\n";
	print OUT join "\n", @testseq;
	close OUT;
	## run clustalw
	!(`$alignmentprog -infile=$tmpdir/$localid.ref.fa -output=fasta -outfile=$tmpdir/$localid.ref.aln 2>&1 >$tmpdir/$localid.ref.log`) or die "error running clustalw\n";
	## get the alignment score
	my $score =  `$grepprog "Alignment Score" $tmpdir/$localid.ref.log |$sedprog -e 's/[^0-9]//g'`;
	if (!$score) {
		die "error in determining alignment score! Problem with ClustalW\n";
	}
	chomp $score;
	return ($score);
}
######################3
sub determineRefspecFinal {
	my ($query_name, @refspec) = @_;
	my $refspec_final;
	## now get the id and the sequence used for building the hmm. Note, the latter will be
	## needed at a later step to determine the best hit
	my @original;
	my $ac = 0;
	for (my $i = 0; $i < @refspec; $i++) {
		$fafile =~ s/\|/\\\|/g;
		@original = `$grepprog -A 1 "^>$query_name|$refspec[$i]" $fafile |$sedprog -e "s/.*$refspec[$i]\|//"`;
		chomp @original;

		if (@original > 0) {
			$refspec_final->[$ac]->{refspec} = $refspec[$i];
			$refspec_final->[$ac]->{searchdb} = "$blastpath/$refspec[$i]/$refspec[$i]" . $blastapp;
			## now allow for more than one sequence per core-ortholog cluster and species
			$refspec_final->[$ac]->{orthocount} = 0;
			for (my $j = 0; $j < @original; $j+= 2) {
				$refspec_final->[$ac]->{refid}->[$refspec_final->[$ac]->{orthocount}] = $original[$j];
				$refspec_final->[$ac]->{sequence}->[$refspec_final->[$ac]->{orthocount}] = $original[$j+1];
				$refspec_final->[$ac]->{orthocount} += 1;
			}
			$ac++;
			@original = qw();
			if (!defined $strict and !defined $relaxed) {
				## one reftaxon is enough
				last;
			}
		}
		else {
			printOUT("original sequence not be found with grepping for ^>$query_name|$refspec[$i]. Proceeding with next refspec\n");
		}
	}
	if (! defined $refspec_final->[0]->{refid}) {
		print "original sequence not found\n";
		return (0, $refspec_final);
	}
	## now print some wordy information...
	if (!defined $strict and !defined $relaxed) {
		printOUT("REFSPEC is $refspec_final->[0]->{refspec}\n");
	}
	return(1, $refspec_final);
}

############## co-ortholog prediction using a alignment score criterion as in InParanoid.
sub identifyCoorthologsProt{
	my ($localid, $spec) = @_;
	my @get;
	my @infofile;
	my $protobject = $fileobj->{$spec}->{prot}; #this is an array ref
	my $idobject = $fileobj->{$spec}->{ids};
	my @genes2check = @$idobject;
	my $refseq = $fileobj->{$spec}->{refprot};
	my $refid = $fileobj->{$spec}->{refid};
	if ($estflag) {
		$refid =~ s/-[0-9]+$//;
	}
	my $refspec_final = $fileobj->{$spec}->{refspec_final};
	## initialize the array with the sequences to be aligned with the reference sequence
	my @out = qw();
	my @hitids = qw();
	push @out, ">$refspec_final";
	push @out, $fileobj->{$spec}->{refseq}->[0];
	for (my $i = 0; $i < @genes2check; $i++) {
		my $seq = $protobject->[$i];
		chomp $seq;
		push @out, ">" . $spec .'|' . $genes2check[$i];
		push @out, $seq;
	}
	## writing sequences to file
	my $tmpdirTmp = $tmpdir; $tmpdirTmp =~ s/\\//g;
	open (OUT, ">$tmpdirTmp/$localid.orth.fa") or die "failed to open $localid.orth.fa\n";
	print OUT join "\n", @out;
	close OUT;
	## aligning sequences
	if ($alignmentprog_co eq 'mafft-linsi'){
		`mafft --maxiterate 1000 --localpair --anysymbol --quiet $tmpdir/$localid.orth.fa > "$tmpdirTmp/$localid.orth.aln"`;
	}
	elsif ($alignmentprog_co eq 'muscle') {
		`$alignmentprog_co -quiet -in $tmpdir/$localid.orth.fa -out "$tmpdirTmp/$localid.orth.aln"`;
	}
	else {
		die "$alignmentprog_co is neither mafft-linsi nor muscle\n";
	}
	if (! -e "$tmpdirTmp/$localid.orth.aln") {
		die "something wrong running $alignmentprog_co\n";
	}
	## do the matrix caluclation
	my $in = Bio::AlignIO->new(-format => 'fasta',
	-file   => "$tmpdirTmp/$localid.orth.aln");
	my $aln = $in->next_aln;
	my $pepstats = Bio::Align::ProteinStatistics->new();
	my $kimura = $pepstats->distance(-align => $aln,
	-method => 'Kimura');
	## do the evaluation
	### get the pairwise distances to the yeast sequences
	#### get the represenative id
	my $smallestdist = $kimura->get_entry("$refspec_final", "$spec|$refid");
	push @get, $spec.'|'.$refid;
	push @infofile, $spec.'|'.$refid.'|'.$smallestdist.'|'.1;
	printOUT("smalles dist is $smallestdist, besthit is $refid from $spec\n");
	## now get any other hit protein that is closer to besthit than the representative seq
	## is to the refspec
	my $count = 0; #this counter keeps track of the entries in the coorthologs field of $fileobj->{$spec}
	for (my $i = 0; $i < @genes2check; $i++) {
		if ($genes2check[$i] ne $refid) {
			my $dist = $kimura->get_entry("$spec|$refid", "$spec|$genes2check[$i]");
			if ($dist <= $smallestdist or ($dist =~ /NaN/  and defined $nonoverlappingCO)) {
				printOUT("co-ortholog detected: $genes2check[$i] with distance $dist compared to $smallestdist of $refid\n");
				$fileobj->{$spec}->{coorthologs}->[$count] = $i;
				$count++;
				push @infofile, $spec.'|'.$genes2check[$i].'|'.$dist.'|'.0;
			}
			else {
				printOUT("co-ortholog rejected: $genes2check[$i] with distance $dist compared to $smallestdist of $refid\n");
			}
		}
	}
	my $counter = 0;
	while (defined $fileobj->{$spec}->{coorthologs}->[$counter]){
		my $index = $fileobj->{$spec}->{coorthologs}->[$counter];
		$counter ++;
	}
	printOUT(join "\n", @infofile);
}
######## sub checkCoorthologRef
sub checkCoorthologRef {
	## relevant steps are
	## 1) get query sequence and query id,
	## 2) get refseq
	## 3) get seq for best blast hit
	## compute the distance query - best blast hit and best blast hit - reference seq
	## return '1' if d(q,b)>d(r,b), else return '0';
	my ($localid, $query, $best, $ref) = @_;
	open (OUT, ">$tmpdir/$localid.co.fa") or die "failed to open $localid.co.fa\n";
	print OUT ">query\n$query\n>best\n$best\n>ref\n$ref\n";
	close OUT;
	## aligning sequences
	if ($alignmentprog_co eq 'mafft-linsi'){
		`mafft --maxiterate 1000 --localpair --anysymbol --quiet $tmpdir/$localid.co.fa > "$tmpdir/$localid.co.aln"`;
	}
	elsif ($alignmentprog_co eq 'muscle') {
		`$alignmentprog_co -quiet -in $tmpdir/$localid.co.fa -out "$tmpdir/$localid.co.aln"`;
	}
	else {
		die "$alignmentprog_co is neither mafft-linsi nor muscle\n";
	}
	if (! -e "$tmpdir/$localid.co.aln") {
		die "something wrong running $alignmentprog_co in routine checkCoorthologRef\n";
	}
	## do the matrix caluclation
	my $in = Bio::AlignIO->new(-format => 'fasta',
	-file   => "$tmpdir/$localid.co.aln");
	my $aln = $in->next_aln;
	my $pepstats = Bio::Align::ProteinStatistics->new();
	my $kimura = $pepstats->distance(-align => $aln,
	-method => 'Kimura');
	## do the evaluation
	### get the pairwise distances to the yeast sequences
	#### get the represenative id
	my $querydist = $kimura->get_entry('query', 'best');
	my $refdist = $kimura->get_entry('best','ref');
	if (($querydist > $refdist) or ($querydist == 0 and $refdist == 0)){
		return(1, $querydist, $refdist);
	}
	else {
		return(0, $querydist, $refdist);
	}
}
####### sub findORF
sub findORF{
	my ($est, $prot, $frame) = @_;
	if ($frame eq '-') {
		$est = revComp($est);
	}
	### debuggin IUB code
	if ($est =~ /[^AGCT]/i) {
		$est =~ s/[^AGCTagct]/n/g;
	}
	printOUT("\trunning genewise using frame $frame\n");
	my $gw = run_genewise_hamstr->new($est, $prot, "$tmpdir");
	my $translation = $gw->translation;
	return ($translation, $est);
}
####### sub printOUT
sub printOUT {
	my $message = shift;
	if (!defined $silent) {
		print $message;
	}
	return();
}
###### sub getLag
sub getLag {
	print "\nInside getlag\n";
	my ($hits, $hitcount) = @_;
	my $minScore = $hits->[$hitcount-1]->{hmmscore};
	my $maxScore = $hits->[0]->{hmmscore};
	if ($minScore == $maxScore) {
		## there is nothing to do, since there is either only one hit, or all hits have the same
		## hmmscore. Return the value of $hitcount.
		return($hitcount, 1);
	}
	## debug
	else {
		print "hitcount is $hitcount, max is $maxScore, min is $minScore\n";
		my @yData = qw();
		my @xData = qw();
		my @xDataLog = qw();
		## now we generate a reversed list of the normalized bitscores
		for (my $i = 0; $i < $hitcount; $i++) {
			push(@yData, 1 - ($hits->[$i]->{hmmscore} - $minScore)/($maxScore - $minScore));
			push(@xDataLog, log(0.1*($i+1)));
			push(@xData, (0.1*($i+1)));
		}
		## The module requires a sufficient amount of trailing 1 to measure the lag point,
		## so we just append them
		for (my $i = $hitcount; $i < ($hitcount+20); $i++) {
			push(@yData, 1);
			push(@xData, 0.1*($i));
			push(@xDataLog, 0.1*($i));
		}
		### calculate end point of lag phase
		my $R = Statistics::R->new();
		# set variables for R
		my $lagPoint = computeLagPoint($R, \@xDataLog, \@yData);
		if ($lagPoint eq 'NA'){
			print "Least square fit to data failed! Trying log-transformed data.\n";
			my $lagPoint = computeLagPoint($R, \@xDataLog, \@yData);
		}
		### compute the cutoff
		if ($lagPoint eq 'NA') {
			return();
		}
		else {
			my $hitLimitGetLag;
			print "limit is $lagPoint. Abs is " . abs($lagPoint) . "\n";
			for (my $i = 0; $i < @xData; $i++) {
				if ($xData[$i] > abs($lagPoint)) {
					$hitLimitGetLag = $i + 1;
					print "Setting hl to $hitLimitGetLag\n";
					last;
				}
			}
			print "hitlimit in getLag is $hitLimitGetLag\n";
			return ($hitLimitGetLag, $lagPoint);
		}
	}
}
##########################
sub getHitLimit {
	my ($hits, $hitcount) = @_;
	my $hitLimitLoc = 0;
	my $maxScore = $hits->[0]->{hmmscore};
	my $limit = $maxScore / 100 * (100 - $scoreCutoff);
	for (my $i = 0; $i < $hitcount; $i++) {
		if ($hits->[$i]->{hmmscore} >= $limit) {
			$hitLimitLoc++;
		}
		else {
			last;
		}
	}
	return ($hitLimitLoc, $limit);
}
## debug
##########################
sub printDebug{
	my @message = @_;
	if ($debug){
		print join "\n", @message;
		print "\n";
	}
}
##########################
sub computeLagPoint {
	my ($R, $xdata, $ydata) = @_;
	$R->set( 'x', \@$xdata);
	$R->set( 'y', \@$ydata );
	# define function
	$R->run( q`func = function(t,params){ 1/(1 + exp(4 * params[1] * (params[2] - x) + 2)) }`);
	# do Nonlinear Least Squares Fitting
	$R->run(q`try <- try(nls(y ~ 1/(1 + exp(4 * mean * (lamda - x) + 2)),
	start = list(mean=1.4, lamda=0.5),
	control = list(maxiter=500)), TRUE)`);
	$R->run(q`if(class(try) != "try-error"){
		f = nls(y ~ 1/(1 + exp(4 * mean * (lamda - x) + 2)),
		start = list(mean=1.4, lamda=0.5),
		control = list(maxiter=500))
		p = coef(f)
		lagPoint = p[2]
	} else {
		lagPoint = "NA"
	}`);


	### return lag point
	my $lagPoint = $R->get('lagPoint');
	return($lagPoint);
}
