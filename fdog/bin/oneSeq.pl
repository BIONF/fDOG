#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;
use File::Copy qw(move);
use File::Basename;
use File::Path;
use File::Path qw/make_path/;
use File::Path 'rmtree';
use File::Which;
use lib dirname(__FILE__);
use Parallel::ForkManager;
use IO::Handle;
use Getopt::Long;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;
use Bio::TreeIO;
use Bio::Tools::Run::StandAloneBlast;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Term::Cap;
use POSIX;

use Capture::Tiny qw/capture/;
use IPC::Run qw( run timeout );
use Time::HiRes;
use List::Util qw(shuffle);
use Cwd;
use Cwd 'abs_path';
use Array::Utils qw(:all);
use Try::Tiny;

my $startTime = gettime();

# Copyright (C) 2009 INGO EBERSBERGER, ebersberger@bio.uni-frankfurt.de
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

# PROGRAM DESCRIPTION: oneSeq.pl is a program for targeted ortholog search in protein sequence data.

# PROGRAM HISTORY
## This script is based on a perl script authored by Peter Schmitzberger in the course
## of his Master's project at the CIBIV, MFPL, Vienna, Austria

## MODIFIED: 13. Aug. 2015 - solved path issues. Script will now work together with
## HaMStR 13.2.5

## Modified 14. Aug. 2015 - added the options -outpath and -hmmpath and a more refined
## clean up after the search, in cases where a custom outpath has been chosen.

## Modified: 01. Feb. 2016: restructured major parts

## Modified: 04. Feb. 2016: - Additions - feature architecture similarity (fas) score support
##                                      - alternations in the program flow
##                                      - global/local option for alignments
##                                      - additional options
##                                      - autocleanup
##                                      - ENV SWATDIR for alignment support (local copy required or use SWATDIR=/home/holger/appz/phredphrap)
##                                      - if you run oneSeq.pl in DB mode, please adapt /bin/run-query.sh to your username and passwort
##                                      -

## Modified 07. Aug. 2017: - Changes:   - change of alignment program, swat replaced by
##                                        ssearch (local:local) and glsearch (global:local) and ggsearch (global:global)
##                                      - selection of best fitting ortholog candidate modified
##                                      - coreFilter: strict, relaxed and none
##

## Modified 19. Jan. 2018: - Additions 	- added option to prioritize closer taxon if two taxa have a similar score
##										- after a taxon has been choosen every taxa closer gets discarded in the next cycles
##										- added commandline parameter to choose the deviation allowed for two taxa to be considered similar
##

## Modified 09. Feb. 2018: - Changes 	- Now the HaMStR candidate search climbs the tree and evaluates only one taxon at a time
##										- The FAS score for a candidate will now only be calculated, if the alignment score is high enough,
##										  to top the current best candidate
##										- If a candidate reaches the maximum score the search stops and a new round starts
##										- If a candidate is within deviation range of the maximum score only the taxa, which are on the same tree branch,
##										  will get evaluated and then the search gets canceled and a new round starts

## Modified 24. Nov. 2018: Release      - release oneSeq v1.3.1
##                                      - Not included feature/feature-updated-fas-util

## Modified 19. July 2019: - Changes	- added option to run muscle instead of mafft

## Modified 22. July 2019: - invoked priority mode for the fas score computation if t = 30

## Modified 2. Dec. 2019
## Bug Fix: Check for taxa with invalid NCBI Taxonomy Id runs now properly and crashes are avoided
## Implemented cleanup of the core ortholog directory to avoid accumulation of feature annotations

## Modified 05. Feb. 2020 (Vinh):   - added option to set number of CPUs for FAS annotation
##									- input faste file must not be present in data folder or working directory
##									- output files will be stored either in user defined directory set via -outpath option, or in working directory by default

## Bug fix 14. April 2020 (Ingo):	- fixed bug that inactivated the -append option

## Modified 14. April 2020 (Vinh): - added option for using user-defined blast_dir, genome_dir and weight_dir
##									- reference species (and taxa for core-set compilation) specified from blast_dir

## Modified 16. Juni 2020 (Vinh): major change in FAS score calculation (v1.7.0)
##									- no need for profile_prog, architecture_prog and visualsPath
##									- final FAS score calculation is done using hamstrFAS

## Modified 16. Juni 2020 v1.7.1 (Vinh)	- replace greedyFAS by calcFAS
## Modified 07. July 2020 v1.7.2 (Vinh)	- check if FAS executable
## Modified 10. July 2020 v1.7.3 (Vinh)	- solved problem when gene ID contains PIPE
## Modified 13. July 2020 v1.8.0 (Vinh)	- added initial check, no longer use .mod files
## Modified 22. July 2020 v1.9.0 (Vinh)	- moved tmp blast files to output folder and delete them when finished
## Modified 27. Aug 2020 v2.1.0 (Vinh)	- option to input newick tree for search taxa
## Modified 07. Sep 2020 v2.2.0 (Vinh)	- append seed sequence to output extended.fa if no ortholog was found in refspec
## Modified 22. Sep 2020 v2.2.1 (Vinh)	- make sure that seed sequence always at the beginning of extended.fa output
## Modified 23. Sep 2020 v2.2.3 (Vinh)	- use full taxonomy name instead of abbr taxon name for LOG
## Modified 01. Dec 2020 v2.2.4 (Vinh)	- fixed bug while creating final extended.fa (and replaced grep and sed by bioperl)

############ General settings
my $version = 'oneSeq v.2.2.4';
##### configure for checking if the setup.sh script already run
my $configure = 0;
if ($configure == 0){
	die "\n\n$version\n\nPLEASE RUN fdog.setup BEFORE USING fdog\n\n";
}
##### hostname
my $hostname = `hostname`;
chomp $hostname;
#############
my $termios = new POSIX::Termios; $termios->getattr;
my $ospeed = $termios->getospeed;
my $t = Tgetent Term::Cap { TERM => undef, OSPEED => $ospeed };
my ($norm, $under, $bold) = map { $t->Tputs($_,1) } qw/me md us/;
#### Paths
my $path = abs_path(dirname(__FILE__));
$path =~ s/\/bin//;
$path =~ s/\/$//;
printDebug("Path is $path");

#### Programs and output
my $sedprog = 'sed';
my $grepprog = 'grep';
my $readlinkprog = 'readlink';

my $globalaligner = 'ggsearch36';
my $glocalaligner = 'glsearch36';
my $localaligner = 'ssearch36';
my $fasta36Path = which('fasta36');
if ( !(defined $fasta36Path) || $fasta36Path eq "") {
	$globalaligner = $path.'/bin/aligner/bin/'.'ggsearch36';
	$glocalaligner = $path.'/bin/aligner/bin/'.'glsearch36';
	$localaligner = $path.'/bin/aligner/bin/'.'ssearch36';
	unless (-e $globalaligner) {
		print "fasta36 not found! Please install it before using fdog!\n";
		exit();
	}
}

my $algorithm = "blastp";
my $blast_prog = 'blastp';
my $outputfmt = 'blastxml';
my $eval_blast_query = 0.0001;
my $filter = 'T';
my $annotation_prog = "annoFAS";
my $fas_prog = "calcFAS";
my $fdogFAS_prog = "fdogFAS";

##### ublast Baustelle: not implemented yet
my $runublast = 0;
my $ublast = 0;
my $accel = 0.8;

############ database connection details
my $dbname="";
my $username="";
my $pw="";
my $database = "DBI:mysql:database=dbdmpng;host=$dbname";
my $getThemAll = 0;
my $updateBlast_dir = 0;

############ directory paths
my $currDir = getcwd;
my $coreOrthologsPath = "$path/core_orthologs/";
my $outputPath = $currDir; #"$path/output"; ## DEFAULT OUTPUT PATH
my $hamstrPath = "$path/bin/hamstr";
my $homeDir = $path;
my $alignmentscoreMatrix = "BP62"; ## opt given by ssearch and glsearch [codaa.mat idnaa.mat P250 P120 BL50 MD40 MD20 MD10 BL62 BL80 BP62 VT160 OPT5]
my $genome_dir = "$path/genome_dir";
my $taxaPath = "$genome_dir/";
my $blastPath = "$path/blast_dir/";
my $idx_dir = "$path/taxonomy/";
my $dataDir = $path . '/data';
my $weightPath = "$path/weight_dir/";

my @defaultRanks = (
	'superkingdom', 'kingdom',
	'superphylum', 'phylum', 'subphylum',
	'superclass', 'class', 'subclass', 'infraclass',
	'superorder', 'order', 'suborder', 'parvorder', 'infraorder',
	'superfamily', 'family', 'subfamily',
	'tribe', 'subtribe',
	'genus', 'subgenus',
	'species group', 'species subgroup', 'species'
);

################## some variables
my $finalOutput;
my $dbHandle;
my $core_hitlimit = 3; # number of hmm hits to consider for reblast during core set generation
# number of hmm hits to consider for reblast during final ortholog search.
# Note, this limits the number of co-orthologs that can be found.
my $hitlimit = 10;
## lagPhase test. Setting the autolimit option to decide from the score distribution how many hits to evaluate.
my $autoLimit;
my $scoreThreshold = 1; # evaluate only hmmsearch hits whose score is within the 10% margin of the best hmmsearch hit
my $scoreCutoff = 10; #value in percent of the hmmscore of the best hit
# Setup for FAS score support (FAS support is used by default)
# Note, fas_t is set to 0.75 by default. Changes will influence sensitivity and selectivity
my $fas_support = 1;
my $countercheck = 0;
my $fasoff      = 0;
my $fasstrict   = 0;
my $fas_T       = 0.75;
my $priThreshold = '-t 30';
my %profile     = ();
my %fas_score_keeper = ();
my $eval_filter = 0.001;
my $inst_eval_filter = 0.01;

my $help;
my @profile = qw();
# my $showTaxa;
my $refSpec;
my $seqFile = '';
my $seqId= '';
my $seqName;
my $minDist;
my $maxDist;
my $minCoreOrthologs;
my $coreTaxa;
my $strict;
my $force = 0;
my $group;
my $groupNode;
my $blast;
my $batch;
my $blastNode;
my $representative;
my $core_rep;
my $debug;
my $corestrict;
my $inputSeq = "";
my $rbh;
my $append = 0;
# Note, the evalue defaults ($eval_blast, $eval_hmmer) will be relaxed for final ortholog run by $eval_relaxfac
my $eval_blast = 0.00001; #1E-5
my $eval_hmmer = 0.00001; #1E-5
my $eval_relaxfac = 10; #checked in checkInput
my $coreOnly;
my $cpu = 1;    #sets number of forks for final ortholog search (can be set via option -cpu=<>)
my $corecpu = 1;    #sets number of forks for core-ortholog assembly (MUST BE 1, due to directed search process through the tree)
my $hyperthread;
my $silent;
my $checkcoorthologsref;
my $cccr;
my $tree;
my $wholeTree;
my $treeDelFlag;
my $currentNoRankDistNode;
my $currentChildsToIgnoreNode;
my $currentDistNode;
my @logOUT = qw();
### Details about the alignment strategy
### Note, the alignment strategy can be local, glocal, or global
### Default: local
my $local;
my $global;
my $glocal;
my $core_filter_mode;
my $dbmode = 0;         ## default run in dbmode. consider setting this in the configure step
my $vlevel = 2;         ## verbosity level
my @taxonlist = qw();
my @refTaxonlist = qw();
my $seqio_object;
my %taxa;
my %refTaxa;
my $autoclean;
my $getversion;
my $coreex; ## flag to set when the core set already exists
my $addenv;
my $ignoreDistance = 0; 	## flag to normalise the score by the distance in the tree
my $distDeviation = 0.05; 	## Span in which a score is consideren similar
my $breakAfter = 5; 		## Number of Significantly bad candidates after which the current run cancels
my %hashTree;
my $aln = 'muscle';
my $searchTaxa;
################# Command line options
GetOptions (
	"h"                 => \$help,
	"append"	=> \$append,
	# "showTaxa"          => \$showTaxa,
	"refSpec=s"         => \$refSpec,
	"db"                => \$dbmode,
	"filter=s"          => \$filter,
	"seqFile=s"   => \$seqFile,
	"seqId=s"           => \$seqId,
	"seqName=s"         => \$seqName,
	"silent"            => \$silent,
	"minDist=s"         => \$minDist,
	"maxDist=s"         => \$maxDist,
	"coreOrth=i"        => \$minCoreOrthologs,
	"coreTaxa=s"        => \$coreTaxa,
	"strict"            => \$strict,
	"rbh"               => \$rbh,
	"evalBlast=s"      => \$eval_blast,
	"evalHmmer=s"      => \$eval_hmmer,
	"evalRelaxfac=s"       => \$eval_relaxfac,
	"checkCoorthologsRef"       => \$checkcoorthologsref,
	"coreCheckCoorthologsRef"   => \$cccr,
	"hitlimitHamstr=s"          => \$hitlimit,
	"coreHitlimitHamstr=s"      => \$core_hitlimit,
	"autoLimitHamstr"	=> \$autoLimit,
	"scoreCutoff=s" => \$scoreCutoff,
	"scoreThreshold" => \$scoreThreshold,
	"coreRep"           => \$core_rep,
	"coreStrict"        => \$corestrict,
	"coreOnly"          => \$coreOnly,
	"group=s"           => \$group,
	"blast"             => \$blast,
	"batch=s"           => \$batch,
	"fas"               => \$fas_support,
	"countercheck"      => \$countercheck,
	"fasoff"            => \$fasoff,
	"coreFilter=s"     => \$core_filter_mode,
	"minScore=s"        => \$fas_T,
	"local"             => \$local,
	"global"            => \$global,
	"glocal"            => \$glocal,
	"rep"               => \$representative,
	"cpu=s"             => \$cpu,
	"outpath=s"         => \$outputPath,
	"hmmpath=s"         => \$coreOrthologsPath,
	"blastpath=s"         => \$blastPath,
	"searchpath=s"         => \$genome_dir,
	"weightpath=s"         => \$weightPath,
	"debug"             => \$debug,
	"coreHitlimit=s"   => \$core_hitlimit,
	"hitlimit=s"        => \$hitlimit,
	"force"             => \$force,
	"cleanup"           => \$autoclean,
	"addenv=s"          => \$addenv,
	"version"           => \$getversion,
	"reuseCore"        => \$coreex,
	"ignoreDistance"	=> \$ignoreDistance,
	"distDeviation=s"	=> \$distDeviation,
	"aligner=s"	=> \$aln,
	"hyperthread" => \$hyperthread,
	"searchTaxa=s" => \$searchTaxa
);

$outputPath = abs_path($outputPath);
unless (-d $coreOrthologsPath) {
	make_path($coreOrthologsPath);
}
$coreOrthologsPath = abs_path($coreOrthologsPath)."/";
$blastPath = abs_path($blastPath)."/";
$weightPath = abs_path($weightPath)."/";
$genome_dir = abs_path($genome_dir)."/";
$taxaPath = $genome_dir;

############# do initial check
if (!defined $help && !defined $getversion) { #} && !defined $showTaxa) {
	print "Validity checking....\n";
	my $checkStTime = gettime();
	initialCheck($seqFile, $seqName, $blastPath, $taxaPath, $weightPath, $fasoff);
	print "Check finished in " . roundtime(gettime() - $checkStTime). " sec!\n";

	if (!defined $coreex) {
		if (!grep(/$minDist/, @defaultRanks)) {
			die "ERROR: minDist $minDist invalid!\n";
		}

		if (!grep(/$maxDist/, @defaultRanks)) {
			die "ERROR: maxDist $maxDist invalid!\n";
		}

		if (!defined $minCoreOrthologs) {
			die "ERROR: coreOrth not defined (must be integer)!";
		}
	}
}

############# show version
if ($getversion){
	print "You are running $version\n";
	print "This version supports FAS comparison.\n";
	exit;
}

############# show help
if($help) {
	my $helpmessage = helpMessage();
	print $helpmessage;
	exit;
}

############# connect to the database
if ($dbmode) {
	$dbHandle = DBI->connect($database, $username, $pw)
	or die "Can not open the database!";
}

############# show all taxa
# if ($showTaxa) {
# 	#get all taxa from database
# 	#hash example: sacce_2336 -> NCBI ID for sacce_2336
# 	printTaxa();
# 	exit;
# }

#switched from online version to flatfile because it is much faster
#taxon files can be downloaded from: ftp://ftp.ncbi.nih.gov/pub/taxonomy/
my $indexStart = gettime();
print "Please wait while the taxonomy database is indexing...\n";
my $db = Bio::DB::Taxonomy->new(-source    => 'flatfile',
	-nodesfile => $idx_dir . 'nodes.dmp',
	-namesfile => $idx_dir . 'names.dmp',
	-directory => $idx_dir);
my $indexTime = gettime() - $indexStart;
my $db_bkp = $db;
print "Indexing done in ",roundtime($indexTime)," sec!\n";

%taxa = getTaxa();
%refTaxa = getRefTaxa();
## debugging message
my $taxcount = keys(%taxa);
printDebug("receiving hash of taxa with $taxcount elements from sub getTaxa");
###
for (keys %taxa){
	printDebug("value of $_ is $taxa{$_}");
}

my $outputFa =  $coreOrthologsPath . $seqName . "/" . $seqName . ".fa";
my $outputAln = $coreOrthologsPath . $seqName . "/" . $seqName . ".aln";
my $tmpdir = $outputPath . '/' . $seqName . '/tmp';
make_path($tmpdir);
checkOptions();
createFoldersAndFiles($outputFa, $seqName, $inputSeq, $refSpec);

my $curCoreOrthologs = 0;
my $hamstrSpecies = $refSpec;
my $addedTaxon = $refSpec;
my $noMoreOrthologs = 0;
my $coremode;
my %finalcontent;
my %candicontent;
my $maxAlnScore = 0;

# create weight_dir in oneseq's home dir (used for annotations,weighting,feature extraction)
# get annotations for seed sequence if fas support is on
if ($fas_support){
	if (!$weightPath) {
		createWeightFolder();
	}
	getAnnotation($outputFa);
}

my $coreStTime = gettime(); #time;
#core-ortholog search
if (!$coreex) {
	print "\nCore compiling...\n";
	$coremode = 1;
	$taxaPath = $blastPath;
	#### moved from above
	my $taxBuildSt = gettime();
	unless ($silent) {
		print "Building up the taxonomy tree...\n";
	}
	push @logOUT, "Building up the taxonomy tree...\n";
	$tree = getRefTree();
	$treeDelFlag = 0;
	if($group) {
		foreach($tree->get_nodes()) {
			if($_->id == $groupNode->id) {
				$groupNode = $_;
			}
		}
		$tree->set_root_node($groupNode);
	}
	unless ($silent) {
		print "Finished building the taxonomy tree in ". roundtime(gettime() - $taxBuildSt) ." sec\n";
	}
	push @logOUT, "Finished building the taxonomy tree in ". roundtime(gettime() - $taxBuildSt) ." sec\n";
	## Tree without deletions
	$wholeTree = getRefTree();
	if($group) {
		foreach($wholeTree->get_nodes()) {
			if($_->id == $groupNode->id) {
				$groupNode = $_;
			}
		}
		$wholeTree->set_root_node($groupNode);
	}
	## initialise control nodes
	$currentDistNode = $wholeTree->find_node(-ncbi_taxid => $refTaxa{$refSpec});
	$currentNoRankDistNode = $currentDistNode->ancestor; ## the node from which the distance to other species will be calculated
	$currentChildsToIgnoreNode = $currentDistNode; ## the node containing all child species which will not be included in the candidates file

	%hashTree = buildHashTree();
	removeMaxDist();
	printDebug("Subroutine call removeMinDist\nRefspec is $refSpec\nTaxon is $refTaxa{$refSpec}\n");
	$treeDelFlag = removeMinDist($refTaxa{$refSpec});
	#### end moved from above

	if ($ignoreDistance){
		$distDeviation = 0;
		$breakAfter = -1;
	}

	## some variables used later
	my $firstRun = 1;

	while (get_leaves($tree, $treeDelFlag) > 0 && $curCoreOrthologs < $minCoreOrthologs && $noMoreOrthologs == 0) {

		# checking the tree which determines the taxa that are going to be searched for hits
		# printDebug("Subroutine call from core-ortholog compilation\nNumber of leaves is ".get_leaves($tree)."\nCurrent core-orthologs: $curCoreOrthologs\nVar \$noMoreOrthologs is set to $noMoreOrthologs\n");
		if ($debug){
			print "\nTaxonomic Tree as text:\n";
			my $tree_as_string = $tree->as_text("tabtree");
			print $tree_as_string;
			print "\n";
		}

		#generate new aln
		if($curCoreOrthologs > 0) {
			createAlnMsf();
		}

		unless ($silent) {
			print "In round $curCoreOrthologs running hmmbuild on $outputAln\n";
		}
		hmmbuild($coreOrthologsPath.$seqName."/hmm_dir/".$seqName.".hmm", $outputAln);

		## get the max alignment score per ortholog
		printDebug("Discovering maximum alignmentscore");

		## Align every current core ortholog against all curretn core orthologs
		## the maximum found in this alignment is the maximun any other sequence can reach
		copy($outputFa, $outputFa . ".extended") or die "Error, could not copy to file: ". "$outputFa" . ".extended\n";

		## get the max alnscore
		my %maxAlnScores = getCumulativeAlnScores();
		foreach my $score (values %maxAlnScores){
			if ($score > $maxAlnScore){
				$maxAlnScore = $score;
			}
		}
		printDebug("The maximum alignmentscore is: $maxAlnScore");
		clearTmpFiles();

		my $addedTaxon = getBestOrtholog();
		my $addedTaxonName = getTaxonName($addedTaxon);
		print "Added TAXON: $addedTaxon\_$addedTaxonName\n";
		#if a new core ortholog was found
		if($addedTaxon ne "") {
			$hamstrSpecies = $hamstrSpecies . "," . $addedTaxon;

			clearTmpFiles();

			++$curCoreOrthologs;
			printDebug("Subroutine call from core-ortholog compilation\nTaxon is $addedTaxon\nNCBI Id is $refTaxa{$addedTaxon}\n");
			$treeDelFlag = removeMinDist($refTaxa{$addedTaxon});
		}
		else {
			#there are no more core orthologs
			$noMoreOrthologs = 1;
			print "\nThe desired number of core orthologs could not be reached.\n";
		}
	}

	## This is now the final round of alignment and profile hidden Markov model building
	## It concludes the core ortholog set compilation
	if ($curCoreOrthologs < $minCoreOrthologs ){
		print "\nWARNING: The desired number of core orthologs could not be reached. Training with only $curCoreOrthologs sequences\n";
	}
	createAlnMsf();
	hmmbuild($coreOrthologsPath.$seqName."/hmm_dir/".$seqName.".hmm", $outputAln);
}
print "==> Core set compilation finished in " . roundtime(gettime() - $coreStTime). " sec!\n";
push @logOUT, "Core set compilation finished in " . roundtime(gettime() - $coreStTime). " sec!";

#after having calculated the core orthologous set,
#start fdog to find all orthologs
# my $finalOutput = $outputPath . '/' . $seqName . '.extended.fa';
my $orthoStTime = gettime();
if (!$coreOnly) {
	$coremode = 0;
	push @logOUT, "Performing the final ortholog search...";
	print "\nPerforming the final ortholog search...\n";
	my $startTmp = gettime();
	#using $eval_relaxfac to relax the evalues for final search
	my $final_eval_blast = $eval_blast*$eval_relaxfac;
	my $final_eval_hmmer = $eval_hmmer*$eval_relaxfac;

	$taxaPath = $genome_dir;
	my @searchTaxa;
	unless($groupNode) {
		@searchTaxa = keys %taxa;
	} else {
		unless ($searchTaxa) {
			# %taxa = getTaxa();
			# print "GET TAXA TIME: ", roundtime(gettime() - $startTmp),"\n";
			my $tree = getTree();
			# print "GET TREE TIME: ", roundtime(gettime() - $startTmp),"\n";
			if($groupNode) {
				foreach($tree->get_nodes()) {
					if($_->id == $groupNode->id) {
						$groupNode = $_;
					}
				}
				$tree->set_root_node($groupNode);
			}
			foreach (get_leaves($tree)) {
				push(@searchTaxa, @{$_->name('supplied')}[0]);
			}
		} else {
			open(SEARCH, $searchTaxa) || die "Cannot open $searchTaxa file!\n";
			@searchTaxa = <SEARCH>;
			close (SEARCH);
		}
	}
	# print "PREPARE TIME: ", roundtime(gettime() - $startTmp),"\n";

	my $pm = new Parallel::ForkManager($cpu);
	if ($hyperthread) {
		$pm = new Parallel::ForkManager($cpu*2);
	}

	foreach (sort @searchTaxa) {
		chomp(my $searchTaxon = $_);
		my $pid = $pm->start and next;
		my $searchTaxonName = getTaxonName($searchTaxon);
		if (defined($searchTaxonName)) {
			unless ($silent) {
				print $searchTaxon, "\t", $searchTaxonName, "\n";
			} else {
				print $searchTaxonName, "\n";
			}
		}
		runHamstr($searchTaxon, $seqName, $finalOutput, $refSpec, $hitlimit, $representative, $strict, $coremode, $final_eval_blast, $final_eval_hmmer, $aln);
		$pm->finish;
	}
	$pm->wait_all_children;
}
push @logOUT, "Ortholog search completed in ". roundtime(gettime() - $orthoStTime) ." sec!";
print "==> Ortholog search completed in ". roundtime(gettime() - $orthoStTime) ." sec!\n";

## Evaluation of all orthologs that are predicted by the final run
if(!$coreOnly){
	my $fasStTime = gettime();
	my $processID = $$;

	# check if final extended.fa exists
	unless (-e $finalOutput) {
		die "ERROR: Could not find $finalOutput\n";
	}
	# check and add seed to final extended.fa if needed
	addSeedSeq($seqId, $seqName, $coreOrthologsPath, $refSpec, $finalOutput); # BLABLABLABLA

	# calculate FAS scores for final extended.fa
	if ($fas_support) {
		print "Starting the feature architecture similarity score computation...\n";
		my $fdogFAScmd = "$fdogFAS_prog -i $finalOutput -w $weightPath -t $tmpdir -o $outputPath --cores $cpu";
		unless ($countercheck) {
			$fdogFAScmd .= " --bidirectional"
		}
		system($fdogFAScmd)
		# print $fdogFAScmd,"\n";
	} else {
		fasta2profile($finalOutput, $seqName)
	}
	push @logOUT, "FAS calculation completed in " . roundtime(gettime() - $fasStTime). " sec!\n";
	print "==> FAS calculation completed in " . roundtime(gettime() - $fasStTime). " sec!\n";
	if($autoclean){
		print "Cleaning up...\n";
		runAutoCleanUp($processID);
	}
}

## Delete tmp folder
unless ($debug) {
	my $delTmp = "rm -rf $tmpdir";
	system ($delTmp) == 0 or die "Error deleting tmp files in $tmpdir\n";
	my $delcommandTmp = "rm -rf $outputPath/tmp";
	system ($delcommandTmp) == 0 or die "Error deleting tmp files in $outputPath/tmp\n";
}
print "==> fdog finished after " . roundtime(gettime() - $startTime) . " sec!\n";
push @logOUT, "fdog finished after " . roundtime(gettime() - $startTime) . " sec!\n";

#### writing the log
open (LOGOUT, ">$outputPath/fdog.log") or warn "Failed to open fdog.log for writing";
print LOGOUT join "\n", @logOUT;
close LOGOUT;
exit;


######################## SUBROUTINES ########################

#################################
## Clears Temporary files
sub clearTmpFiles {
	#clear temporary result file
	if(-e $outputFa.".extended") {
		unlink($outputFa.".extended");
	}

	#clear all alignment files
	my @files = glob("*.scorefile");
	foreach my $file (@files) {
		unlink($file);
	}
}

sub getCandicontent{
	my %candicontent;
	my $candidatesFile = $outputFa . ".extended";
	if (-e $candidatesFile) {

		########################
		## step: 2
		## setup
		## candidates to hash
		## %candicontent keeps info about all candidates (header and sequence)
		open (CANDI, "<".$candidatesFile) or die "Error: Could not find $candidatesFile\n";
		my $head;
		%candicontent = ();
		while(<CANDI>){
			my $line = $_;
			chomp($line);
			if ($line =~ m/^>/){
				$line =~ s/>//; # clip '>' character
				$head = $line;
			}else{
				$candicontent{$head} = $line;
			}
		}
		close (CANDI);
	}
	return %candicontent;
}

#################################
## Get the alinment score for the current candidate file
## only works for files holding only one candidate
sub getCumulativeAlnScores{
	chdir($coreOrthologsPath . $seqName);
	my $candidatesFile = $outputFa . ".extended";
	my $scorefile = $$ . ".scorefile";
	my %scores;
	########################
	## step: 1
	## setup
	## set alignment command (glocal, local, or global)
	#local      local:local    ssearch36   Smith-Waterman
	#glocal     global:local   glsearch36  Needleman-Wunsch
	#global     global:global  ggsearch36  Needleman-Wunsch
	my $loclocCommand = "$localaligner \"" . $outputFa . "\" \"" . $candidatesFile . "\" -s " . $alignmentscoreMatrix . " -m 9 -d 0 -z -1 -E 100" . " > " . $scorefile;
	my $globlocCommand = "$glocalaligner \"" . $outputFa . "\" \"" . $candidatesFile . "\" -s " . $alignmentscoreMatrix . " -m 9 -d 0 -z -1 -E 100" . " > " . $scorefile;
	my $globglobCommand = "$globalaligner \"" . $outputFa . "\" \"" . $candidatesFile . "\" -s " . $alignmentscoreMatrix . " -m 9 -d 0 -z -1 -E 100" . " > " . $scorefile;
	########################
	## step: 2
	## setup
	## candidates to hash
	## %candicontent keeps info about all candidates (header and sequence)
	my %candicontent = getCandicontent();

	########################
	## step: 3
	## get alignment scores
	chdir($coreOrthologsPath . $seqName);
	if ($glocal){
		system($globlocCommand);
	}elsif ($global){
		system($globglobCommand);
	}elsif ($local){
		system($loclocCommand);
	}
	########################
	## step: 4
	## collect alignment score
	## keep track about min and max for each query/coreortholog vs candidate set
	my $max = -10000000;
	my $min = 10000000;

	%scores = cumulativeAlnScore($scorefile, \%candicontent);
	return %scores;
}

#################################
## Get the alinment scores for the current candidate file
sub getAlnScores{
	chdir($coreOrthologsPath . $seqName);
	my $candidatesFile = $outputFa . ".extended";
	my $scorefile = $$ . ".scorefile";
	my %scores;

	########################
	## step: 1
	## setup
	## set alignment command (glocal, local, or global)
	#local      local:local    ssearch36   Smith-Waterman
	#glocal     global:local   glsearch36  Needleman-Wunsch
	#global     global:global  ggsearch36  Needleman-Wunsch
	my $loclocCommand = "$localaligner " . $outputFa . " " . $candidatesFile . " -s " . $alignmentscoreMatrix . " -m 9 -d 0 -z -1 -E 100" . " > " . $scorefile;
	my $globlocCommand = "$glocalaligner " . $outputFa . " " . $candidatesFile . " -s " . $alignmentscoreMatrix . " -m 9 -d 0 -z -1 -E 100" . " > " . $scorefile;
	my $globglobCommand = "$globalaligner " . $outputFa . " " . $candidatesFile . " -s " . $alignmentscoreMatrix . " -m 9 -d 0 -z -1 -E 100" . " > " . $scorefile;

	########################
	## step: 2
	## setup
	## candidates to hash
	## %candicontent keeps info about all candidates (header and sequence)
	my %candicontent = getCandicontent();

	########################
	## step: 3
	## get alignment scores
	chdir($coreOrthologsPath . $seqName);
	if ($glocal){
		system($globlocCommand);
	}elsif ($global){
		system($globglobCommand);
	}elsif ($local){
		system($loclocCommand);
	}

	########################
	## step: 4
	## collect alignment score
	## keep track about min and max for each query/coreortholog vs candidate set
	my $max = -10000000;
	my $min = 10000000;

	%scores = cumulativeAlnScore($scorefile, \%candicontent);

	## Normalize Alignment scores (unity-based)
	printDebug("Normalize alignment scores:\n");
	foreach my $key (keys %scores){
		my $score = $scores{$key};
		unless ($silent) {
			print "Cumulative alignmentscore is: $score\n";
		}
		$scores{$key} = $scores{$key} / $maxAlnScore;
		$score = $scores{$key};
		unless ($silent) {
			print "Normalised alignmentscore is: $score\n";
		}
	}
	return %scores;
}

#################################
## Get the fas scores for the current candidate file
sub getFasScore{
	printDebug("Changing to $coreOrthologsPath$seqName", "Candidate file is $outputFa".'.extended');
	chdir($coreOrthologsPath . $seqName);
	my %fas_box;
	my $scorefile = $$ . ".scorefile";
	my $rankscore;

	########################
	## step: 1
	## setup
	## candidates to hash
	## %candicontent keeps info about all candidates (header and sequence)
	my %candicontent = getCandicontent();

	########################
	## step: 2
	## get FAS score
	## fas support: on/off
	if ($fas_support){
		my @candidateIds = keys(%candicontent);
		my ($name,$gene_set,$gene_id,$rep_id) = split(/\|/, $candidateIds[0]);
		unless (-e "$weightPath/$gene_set.json") {
			print "ERROR: $weightPath/$gene_set.json not found! FAS Score will be set as zero.\n";
			$fas_box{$candidateIds[0]} = 0.0;
		} else {
			my $lnCmd = "ln -fs $weightPath/$gene_set.json \"$coreOrthologsPath$seqName/fas_dir/annotation_dir/\"";
			system($lnCmd);
			my $fasOutTmp = `$fas_prog -s \"$coreOrthologsPath$seqName/$seqName.fa\" -q $blastPath/$gene_set/$gene_set.fa --query_id \"$gene_id\" -a \"$coreOrthologsPath$seqName/fas_dir/annotation_dir/\" -o \"$coreOrthologsPath$seqName/fas_dir/annotation_dir/\" --raw --tsv --domain --cpus 1 | grep "#" | cut -f 3,4`;
			my @fasOutTmp = split(/\t/,$fasOutTmp);
			$fas_box{$candidateIds[0]} = $fasOutTmp[1];
		}
	}
	return %fas_box;
}

#################################
## Add fas and alignment score together while using specified filters (coreFilter option)
sub getFilteredRankScore{
	my $alnScore = $_[0];
	my $fasScore = $_[1];
	my $rankscore = 0;
	# $rankscore: keeps alignment and fas score, decider about $bestTaxon
	if ($core_filter_mode){
		if ($core_filter_mode eq "strict"){
			# case 1: filter
			if ($fasScore < $fas_T){
				#eliminate candidate $key
				print "Deleting candidate from list due to insufficient FAS score.\n";
				$rankscore = 0;
			}else{
				#keep
				if ($alnScore){
					$rankscore = $fasScore + $alnScore;
				}else{
					$rankscore = $fasScore;
				}
			}
		}elsif ($core_filter_mode eq "relaxed"){
			# case 2: disadvantage
			if ($fasScore < $fas_T){
				# ignore FAS score for rankscore
				printDebug("Candidate will be disadvantaged.\n");
				if ($alnScore){
					$rankscore = $alnScore;
				}else{
					$rankscore = 0;
				}
			}
			else{
				#keep
				if ($alnScore){
					$rankscore = $fasScore + $alnScore;
				}else{
					$rankscore = $fasScore;
				}
			}
		}
	}else{
		# case 3: no filter
		if($fasScore) {
			if ($alnScore){
				$rankscore = $fasScore + $alnScore;
			}else{
				$rankscore = $fasScore;
			}
		}
	}
	return $rankscore;
}

sub getHeaderSeq{
	my $bestTaxon = $_[0];
	open (EXTFA, $outputFa.".extended");
	my $sequenceLine = 0;
	my $bestSequence = "";

	########################
	## step: 7
	## get best sequence from candidate file
	## (will be added to the model)
	while(<EXTFA>) {
		my $line = $_;
		chomp($line);
		if($sequenceLine == 1) {
			$bestSequence = $line;
			$sequenceLine = -1;
		}

		if($line eq $bestTaxon) {
			$sequenceLine = 1;
		}
	}
	close EXTFA;
	my @best = split '\|', $bestTaxon;
	my $header = $best[0] . '|' . $best[1] . '|' . $best[2];
	return ($header, $bestSequence);
}

## create profile from extended.fa
sub fasta2profile{
	my ($file, $out) = ($_[0], $_[1]);
	my ($fO_base, $fO_path, $fO_suffix) = fileparse( $file, qr/\.[^.]*/ );
	my $outFile = $fO_path.$out.".phyloprofile";
	open(FA, $file);
	open(PPOUT, ">$outFile");
	print PPOUT "geneID\tncbiID\torthoID\n";
	foreach my $line(<FA>) {
		if ($line =~ /^>/) {
			chomp($line);	# test|ANOGA@7165@1|Q7Q3C2|1
			$line =~ s/>//;
			my @lineTMP = split(/\|/, $line);
			my $geneID = $lineTMP[0];
			my @orthoTMP = split(/@/, $lineTMP[1]);
			my $ncbiID = "ncbi".$orthoTMP[1];
			print PPOUT $geneID, "\t", $ncbiID, "\t", $line,"\n";
		}
	}
	close(FA);
	close(PPOUT);
}

## auto clean up can be invoked via the "-cleanup" option
# $processID: given process ID
sub runAutoCleanUp {
	my $processID = $_[0];
	unless ($silent) {
		print "Deleting $outputPath/tmp\n";
	}
	my $delCommandTmp = "rm -rf \"$outputPath/tmp\"";
	system ($delCommandTmp) == 0 or die "Error deleting result files\n";
	my $seedName = $seqName . '_seed';
	my $annopath = $coreOrthologsPath.$seqName."/fas_dir/annotation_dir";
	if( -l "$currDir/$seqFile" ) {
		my $delLnSeedFile = "rm $currDir/$seqFile";
		system ($delLnSeedFile);
	}
	unless ($silent) {
		print "Deleting $annopath\n";
	}
	if (!$fasoff) {
		opendir(ANNODIR, $annopath) or warn "Could not open $annopath in sub runAutoCleanup\n";
		my @annodirs = grep (!/$seedName/, readdir(ANNODIR));
		unless ($silent) {
			print scalar(@annodirs) . " content of $annopath\n";
		}
		for (my $anno = 0; $anno < @annodirs; $anno++){
			if ($annodirs[$anno] ne '..' and $annodirs[$anno] ne '.' and $annodirs[$anno] ne $seqName.".json") {
				unless ($silent) {
					print "Deleting $annopath/$annodirs[$anno]\n";
				}
				rmtree ($annopath."/".$annodirs[$anno]);
			}
		}
		closedir (ANNODIR);
	}
}

## starting annotation_prog for given seed sequence file
# $seedseqFile: fasta file with seed sequence
sub getAnnotation {
	my ($seedseqFile) = ($_[0]);
	my $inputAnno = $seedseqFile;
	$inputAnno =~ s/\|/\\\|/g;
	my $outputAnno = $coreOrthologsPath . $seqName . "/fas_dir/annotation_dir";
	$outputAnno =~ s/\|/\\\|/g;
	my $annotationCommand = "$annotation_prog" . " -i $inputAnno" . " -o $outputAnno --cpus 1" . " --name \"$seqName\""; #" --name " . $seqName . "_seed" . " --cpus 1";
	system($annotationCommand);
}

## determines the reference species and/or the sequence id of the input sequence.
sub determineRef {
	my ($infile, @refspec) = @_;
	#run blast for all available proteomes if the given sequence is not in the database
	unless ($silent) {
		print "One moment please!\nLooking for the most similar sequence to your input sequence.\n";
	}
	my $bestHit->{score} = 1;
	$bestHit->{evalue} = 10;
	my $outname = $$;
	## Baustelle: Currently, we need to loop through all possible taxa to id the best matching one
	for (my $i = 0; $i < scalar(@refspec); $i++) {
		my $curTaxon = $refspec[$i];
		## run the blast search
		printDebug("running blast on $curTaxon");
		my $resultFile = runBlast($seqFile, $dataDir, $outname, $tmpdir, "$blastPath/$curTaxon/$curTaxon");
		my $hits = &getBestBlasthit($tmpdir, $resultFile);
		if (defined $hits and @$hits > 0) {
			#only use the best hit with the index [0]. Note, $hits is an array ref of hashrefs.
			if($hits->[0]->{score} > $bestHit->{score}) {
				$bestHit->{name} = $hits->[0]->{name};
				$bestHit->{score} = $hits->[0]->{score};
				$bestHit->{evalue} = $hits->[0]->{evalue};
				$bestHit->{species} = $curTaxon;
			}
		}
	}
	return($bestHit);
}

sub checkGroup {
	my $group = shift;
	my $node = $db->get_taxon(-name => $group);
	if($node) {
		$groupNode = $node;
	}
	else {
		print "Your selected group " . $group . " was not found in the taxonomic tree... TERMINATING\n";
		exit;
	}
}

#################################
sub checkOptions {
	if($eval_relaxfac < 1){
		# rethink
		if($eval_relaxfac <= 0){
			printOut("\nThe specified factor for evalue relaxation is <= 0. Please see the help text for option -eval_relaxfac. We recommend a factor > 1. Default is 10.\n",1);
			my $answer = '';
			my $breaker = 0;
			while (($answer !~ /[0-9]/i) and ($breaker < 4)) {
				$breaker++;
				$answer = getInput("Please choose a new factor (Integer) for evalue relaxation. [1,100]");
				if (($breaker > 3) and ($answer !~ /[0-9]/i)){
					print "No proper factor given ... exiting.\n";
					exit;
				}
			}
			if ($answer =~ /[0-9]/i) {
				$eval_relaxfac = $answer;
			}
		}
	}
	### check for colision of force and append. Change in favor of append
	if ($force == 1 and $append ==1) {
		$force = 0;
	}
	### check the presence of the pre-computed core set
	if ($coreex) {
		if (! -e "$coreOrthologsPath/$seqName/$seqName.fa") {
			print "You selected the option -reuseCore, but the core ortholog group $coreOrthologsPath/$seqName/hmm_dir/$seqName.hmm does not exist\n";
			exit;
		}
	}
	### begin move up
	### checking reference species
	my $optbreaker = 0;
	while ((!$refSpec or !$refTaxa{$refSpec})  && !$blast) {
		if ($optbreaker >= 3){
			print "No proper refspec given ... exiting.\n";
			exit;
		}
		my $output = '';
		for (my $i = 0; $i < @refTaxonlist; $i++) {
			$output = $output . "[$i]" . "\t" . $refTaxonlist[$i] . "\n";
		}
		### for debug?
		# for (keys %taxa){
		# 	print "value of $_ is \'$taxa{$_}\'";
		# }
		# printDebug("taxa contains $taxa{$refSpec}"); # cannot print this if $taxa{$refSpec} not exists!
		my $refSpecIdx = getInput("\n" . $output . "\n" . "You did not provide a valid reference species ($refSpec). Please choose the number for the reference species your input sequence is derived from", 1);
		$optbreaker++;
		$refSpec = $refTaxonlist[$refSpecIdx];
		checkBlastDb($refSpec, $refSpec);
	}
	### end move up
	### adding new routine to generate the input sequence if -reuseCore has been set
	if ($coreex) {
		my @refseq=`$grepprog -A 1 ">$seqName|$refSpec" $coreOrthologsPath/$seqName/$seqName.fa`;
		chomp @refseq;
		unless ($silent) {
			print "$refseq[0]\n";
		}
		(my $tmp1, my $tmp2, $seqId) = split '\|', $refseq[0];
		if (length($seqId) == 0){
			die "error in retrieving sequence while using -reuseCore\n";
		}
		print "overruling the provided seed sequence since you used the option -reuseCore. Setting seed id to $seqId\n";
		open OUT, (">$currDir/$seqName.fa") or die "could not open $currDir/$seqFile for writing in retrieve refseq\n";
		print OUT join "\n", @refseq;
		close OUT;
		$seqFile = "$seqName.fa";
	}
	### end mod
	### check input file
	$optbreaker = 0;
	while ((length $seqFile == 0) or ((! -e "$currDir/$seqFile") and (! -e "$dataDir/$seqFile"))) {
		if ($optbreaker >= 3){
			print "No proper file given ... exiting.\n";
			exit;
		}
		if (length $seqFile > 0){
			if (-e $seqFile) {
				my @seqFileTMP = split(/\//, $seqFile);
				unless (-e "$currDir/$seqFileTMP[@seqFileTMP-1]") {
					system("ln -fs \"$seqFile\" \"$currDir/$seqFileTMP[@seqFileTMP-1]\"");
				}
				$seqFile = $seqFileTMP[@seqFileTMP-1];
			} else {
				printOut("\nThe specified file $seqFile does not exist!\n",1);
			}
		}
	}
	if (-e "$currDir/$seqFile"){
		$dataDir = $currDir;
		printDebug("Setting datadir to $currDir in sub checkOptions");
	}

	### checking the number of core orthologs. Omit this check if the option -reuseCore has been selected
	$optbreaker = 0;
	while(!$minCoreOrthologs and !$coreex) {
		if ($optbreaker >= 3){
			print "No proper number given ... exiting.\n";
			exit;
		}
		$minCoreOrthologs = getInput("Please specify the desired number of core orthologs!", 1);
		$minCoreOrthologs = checkInt($minCoreOrthologs);
		$optbreaker++;
	}
	## check for blast filter
	if ($blast_prog ne 'blastall'){
		$filter = 'yes' if $filter eq 'T';
		$filter = 'no' if $filter eq 'F';
	}

	$inputSeq = fetchSequence($seqFile, $dataDir);

	## the user has not provided a sequence id, however, the refspec is determined.
	if($seqId eq '') {
		my $besthit;
		if (!$blast){
			## a refspec has been determined
			#use blast to search in the proteome of the specified reference species for the input sequence
			#in order to obtain a valid sequence id
			$besthit = determineRef($seqFile, ($refSpec));
		}
		else {
			$besthit = determineRef($seqFile, @refTaxonlist);
		}
		$seqId = $besthit->{name};
		$refSpec = $besthit->{species};
		my $details = "Evalue: " . $besthit->{evalue};
		printOut("Seq id has been determined as $seqId in $refSpec with $details", 2);
		if($seqId eq '') {
			print "There was no significant hit for your sequence in " . $refSpec . ".\nPlease specify a sequence id on your own.\n";
			exit;
		}
	}

	if($coreTaxa) {
		if(! -e $coreTaxa) {
			print "Please specify a valid file with taxa for the core orthologs search\n";
			exit;
		}
		my @userTaxa = parseTaxaFile();
		my %newTaxa = ();
		foreach (@userTaxa) {
			$newTaxa{$_} = $taxa{$_};
		}
		$newTaxa{$refSpec} = $refTaxa{$refSpec};
		%taxa = %newTaxa;
	}

	if($group) {
		checkGroup($group);
	}

	if(!$seqName) {
		my $i = 0;
		while($i < 7) {
			my $j = chr(int(rand(127)));
			if($j =~ /[a-zA-Z]/) {
				$seqName .=$j;
				$i++;
			}
		}
		print "Your sequence was named: " . $seqName . "\n\n";
	}
	$outputPath = $outputPath . "/$seqName";
	if (! -d "$outputPath"){
		mkdir "$outputPath", 0777  or die "could not create the output directory $outputPath";
	}
	## check whether a result file already exists:
	$finalOutput = $outputPath . '/' . $seqName . '.extended.fa';
	if ($outputPath && -e "$finalOutput"){
		## an ouput file is already existing
		if (!$force && !$append){
			## The user was not aware of an existing output file. Let's ask him
			my $input = '';
			my $breaker = 0;

			while (($input !~ /^[aor]/i) and ($breaker < 4)) {
				$breaker++;
				die "\nAn outputfile $finalOutput already exists. Please consider option -force for overwriting it or option -append for appending to it.\n"
				# $input = getInput("\nAn outputfile $finalOutput already exists. Shall I overwrite it [o], or rename it [r], or [a] append to it?", 1);
				# if (($breaker > 3) and ($input !~ /[aor]/i)){
				# 	print "Please consider option -force or option -append.\n";
				# 	die "No proper answer is given: Quitting\n";
				# }
			}
			# if ($input =~ /[aA]/) {
			# 	$append = 1;
			# 	$force = 0;
			# }
			# elsif ($input =~ /[oO]/){
			# 	$append = 0;
			# 	$force = 1;
			# }
			# else {
			# 	$append = 0;
			# 	$force = 0;
			# }
		}
		if ($force){
			## the user wants to overwrite
			printOut("Removing existing output directory $outputPath", 1);
			rmtree ([ "$outputPath" ]) or die "could not remove existing output directory $outputPath\n";
			mkdir $outputPath or die "could not re-create the output directory $outputPath\n";
		}
		elsif ($append) {
			printOut("Appending output to $finalOutput\n", 1);
			if (-e "$outputPath/$seqName.extended.profile") {
				## read in the content for latter appending
				printOut("Appending output to $outputPath/$seqName.extended.profile", 1);
				open (IN, "<$outputPath/$seqName.extended.profile") or die "failed to open $outputPath/$seqName.extended.profile after selection of option -append\n";
				while (<IN>) {
					chomp $_;
					my @keys = split '\|', $_;
					$profile{$keys[1]} = 1;
				}
			}
			elsif ($fasoff) {
				## no extended.profile file exists but not necessary, because user switched off FAS support -> do nothing
			}
			else {
				printOut("Option -append was selected, but the existing output was incomplete. Please restart with the -force option to overwrite the output");
				exit;
			}
		}
		else {
			printOut("Renaming existing output file to $finalOutput.old", 2);
			my $bu_dir = $outputPath.'_bkp';
			!`mv $outputPath $bu_dir` or die "Could not rename existing output file $outputPath to $bu_dir\n";
			mkdir $outputPath or die "could not recreate $outputPath after renaming the old output\n"
		}
	}

	#### checking for the min and max distance for the core set compilation
	#### omit this check, if the option reuseCore has been selected (added 2019-02-04)
	$optbreaker = 0;
	if (!$coreex) {
		my $node;
		$node = $db->get_taxon(-taxonid => $refTaxa{$refSpec});
		$node->name('supplied', $refSpec);
		if (lc($maxDist) eq "root"){
			$maxDist = 'no rank';
		}
		while (!$maxDist or (checkRank($maxDist, $node) == 0)) {
			if ($optbreaker >= 3){
				print "No proper maxDist given ... exiting.\n";
				exit;
			}
			print "You have not defined a valid maximum distance rank!\n";
			printTaxonomy($node);
			my $in = getInput('Please choose a rank by giving the number in square brackets', 1);
			$optbreaker++;
			$maxDist = parseInput($node, $in);
			print "You selected ". $maxDist . " as maximum rank\n\n";
		}
		while (!$minDist or (checkRank($minDist, $node) == 0)) {
			if ($optbreaker >= 3){
				print "No proper minDist given ... exiting.\n";
				exit;
			}
			print "You have not defined a minimum distant rank!\n";
			printTaxonomy($node);
			my $in = getInput('Please choose a rank by giving the number in square brackets', 1);
			$optbreaker++;
			$minDist = parseInput($node, $in);
			print "You selected " . $minDist . " as minimum rank\n\n";
		}
	}
	$optbreaker = 0;

	#### checking in fas options
	if($fasoff){
		print "You have turned FAS support off. Candidate orthologs will not be evaluated by their FAS.\n";
		# turn FAS support off
		$fas_support = 0;
	}
	## check if user defined fas_T is off limits
	if ($fas_T < 0 or $fas_T > 1){
		print "You chose an odd FAS score filter (-minScore), default is 0.75.\n";
		my $answer = '';
		$optbreaker = 0;
		while ($answer < 0 or $answer > 1) {
			if ($optbreaker >= 3){
				print "No proper fas filter given ... exiting.\n";
				exit;
			}
			$answer = getInput("Please choose a FAS score filter [0,1] between 0 (relaxed) and 1 (stringent):");
			$optbreaker++;
		}
		if ($answer > 0 and $answer < 1) {
			$fas_T = $answer;
		}
	}
	### rather strict fas filter for core orthologs: OFF
	if(!$core_filter_mode){
		unless ($silent) {
			print "No FAS filter for core-orthologs set.\n";
		}
	}elsif($core_filter_mode eq "relaxed"){
		#core ortholog candidates with a FAS score below the threshold will be disadvantaged
	}elsif($core_filter_mode eq "strict"){
		#core ortholog candidates with a FAS score below the threshold will not be considered any more
	}else{
		print "No known filter mode for core-orthologs specified. Continuing with default settings\n";
		$core_filter_mode = 0;
	}

	### check alignment strategy
	if (($local && $global) or ($local && $glocal) or ($global && $glocal)){
		print "Please specify only one alignment strategy!\n";
		print "Possible options are: -glocal, -local, or -global\n";
		print "... exiting.\n";
		exit;
	}elsif(!$local && !$global && !$glocal){
		unless ($silent) {
			print "No specific alignment strategy set. Continuing with local alignments (Smith-Waterman-Algorithm).\n";
		}
		$local = 1;
	}
}

####################### sub check the systematic rank
sub checkRank {
	my $rank = $_[0];
	my $node = $_[1];
	my $rankExists = 0;
	while($node->ancestor && $rankExists == 0) {
		if($node->rank eq $rank) {
			$rankExists = 1;
		}
		$node = $node->ancestor;
	}

	if($node->rank eq $rank) {
		$rankExists = 1;
	}

	return $rankExists;
}

############
## modified by Ingo - Added Option to run Muscle
sub createAlnMsf {
	my $linsiCommand = '';
	if (!defined $aln or $aln eq 'mafft-linsi') {
		my $linsiCommand = "mafft --maxiterate 1000 --localpair --anysymbol --quiet \"" . $outputFa . "\" > \"" . $outputAln . "\"";
	}
	elsif ($aln eq 'muscle') {
		$linsiCommand = "muscle -quiet -in \"" . $outputFa . "\" -out \"" .$outputAln. "\"";
	}
	else {
		die "issues with the msa. You need to select either mafft or muscle\n";
	}
	system($linsiCommand) == 0 or die "Could not run alignment\n$linsiCommand\n";
}

################ creating folders for fas support usage
sub createWeightFolder{
	#create weight_dir in hamstr1seq home dir
	my $weightdir = $path."/"."weight_dir";
	mkdir "$weightdir", 0777 unless -d "$weightdir";
}

################
sub createFoldersAndFiles {
	my ($outputFa, $seqName, $inputSeq, $refSpec) = (@_);
	#create core orthologs directory
	my $dir = $coreOrthologsPath . $seqName;
	if (!$coreex){
		mkdir "$dir", 0755 unless -d "$dir";
		my $header = $seqName . "|" . $refSpec . "|" . $seqId;

		#create FA file
		open (OUTPUT,  ">$outputFa") or die "Error creating fa file $outputFa\n";
		print OUTPUT  ">" . $header . "\n";
		print OUTPUT $inputSeq;
		close OUTPUT;

		#create the Aln file initially only with a single species in there
		open (OUTPUT,  ">$outputAln") or die "Error creating fa file $outputAln\n";
		print OUTPUT  ">" . $header . "\n";
		print OUTPUT $inputSeq;
		close OUTPUT;

		#create the folder for the hmm output
		my $hmmdir = $dir . "/hmm_dir";
		mkdir "$hmmdir", 0755 unless -d "$hmmdir";
	}
	#create the fas_dir for core orthologs if fas support is ON
	if ($fas_support){
		my $fasdir = $dir. "/fas_dir";
		mkdir "$fasdir", 0777 unless -d "$fasdir";

		my $annodir = $fasdir."/annotation_dir";
		mkdir "$annodir", 0777 unless -d "$annodir";
	}
}
#################
sub fetchSequence {
	my ($file, $filepath) = @_;
	if (! defined $filepath){
		$filepath = '.';
	}
	my $seq = "";
	open (INPUT, "<$filepath/$file") or die print "Error opening seq file\n";
	while(<INPUT>) {
		my $line = $_;
		chomp($line);
		unless($line =~ /^\>.*/) {
			$seq = $seq . $line;
		}
	}
	close INPUT;
	$seq =~ s/\s*//g;
	unless ($silent) {
		printOut($seq, 2);
	}
	return $seq;
}
#################################
## choose the ortholog which reaches the highest score
sub getBestOrtholog {
	## max possible score is either one or two
	my $maxScore = 1;
	if ($fas_support){
		$maxScore += 1;
	}

	## get leavs to evaluate
	my @leaves = get_leaves($tree, $treeDelFlag);
	## sort by distance in taxonomy tree
	if (!$ignoreDistance){
		@leaves = sort_leaves(@leaves);
	}
	## don't sort by distance
	else{
		my @unsortedLeaves = @leaves;
		@leaves = qw();
		push @leaves, \@unsortedLeaves;
	}

	## create needed variables
	my $bestTaxon = '';
	my $rankScore = 0;
	my $header = '';
	my $seq = '';
	my $newNoRankDistNode;		## this will be the new Distance node, after a new candidate has been choosen
	my $newChildsToIgnoreNode;	## all leaves under this node will be ignored in future runs, after a new candidate has been choosen
	my $sufficientlyClose = 0;	## flag to break outer loop
	my $candidatesFile = $outputFa . ".extended";

	## iterate over each array with leaves of same distance
	foreach my $array (@leaves) {
		## break loop if a candidate was close to the max score and no more candidates remain with the same distance
		if ($sufficientlyClose){
			unless ($silent) {
				print "Best Taxon is sufficiently close to max score and no more candidates with same distance remain.\nStopping evaluation.\n";
			}
			last;
		}
		## iterate over each leaf with the same distance
		foreach my $key (@$array){
			my $keyName = @{$key->name('supplied')}[0];
			my $nodeId = $wholeTree->find_node(-ncbi_taxid => $refTaxa{$keyName})->id;
			unless ($silent) {
				print "fdog species: " . $key->scientific_name . " - " . @{$key->name('supplied')}[0] . "\n";
			}
			my $coreTaxon = @{$key->name('supplied')}[0];
			my $coreTaxonName = getTaxonName($coreTaxon);
			if (defined($coreTaxonName)) {
				unless ($silent) {
					print $coreTaxon, "\t", $coreTaxonName, "\n";
				} else {
					print $coreTaxonName, "\n";
				}
			}
			runHamstr($coreTaxon, $seqName, $outputFa, $refSpec, $core_hitlimit, $core_rep, $corestrict, $coremode, $eval_blast, $eval_hmmer, $aln);
			## check weather a candidate was found in the searched taxon
			if(-e $candidatesFile) {

				## get found candidates for one taxon in hash to iterate over
				my %candicontent = getCandicontent();

				## get scores in hashes because there might be more than one candidate sequence per taxon
				my %alnScores = getAlnScores();
				my %fas_box;
				my $gotFasScore = 0;
				## iterate over found candidates
				foreach my $candiKey (keys %candicontent){
					## candidates alnScore is high enought, that it would be better with a fasScore of one
					## -> evaluate
					if ($alnScores{$candiKey} > $rankScore * (1 + $distDeviation) - 1){
						if (!$gotFasScore and $fas_support){
							%fas_box = getFasScore();
							$gotFasScore = 1;
						}
						## get rankscore
						my $newRankScore = getFilteredRankScore($alnScores{$candiKey}, $fas_box{$candiKey});
						## candidate is significantly better, than the last one
						if ($newRankScore > $rankScore * (1 + $distDeviation)){ #uninit
						$bestTaxon = ">" . $candiKey;
						$rankScore = $newRankScore;
						($header, $seq) = getHeaderSeq($bestTaxon);
						$newNoRankDistNode = $currentNoRankDistNode;
						$newChildsToIgnoreNode = $currentChildsToIgnoreNode;
						my $newNodeId = $key->id;
						## set new distance nodes, which will replace the old ones given, that this candidate will remain the best
						while (!defined $hashTree{$newNoRankDistNode}{$newNodeId}){
							$newNoRankDistNode = $newNoRankDistNode->ancestor;
							$newChildsToIgnoreNode = $newChildsToIgnoreNode->ancestor;
						}
						unless ($silent) {
							print "New Best Taxon: $bestTaxon\n";
						}
					}
				}
				## candidate has the same distance, as the last one and could be better, with a fasScore of one
				elsif (defined $hashTree{$newNoRankDistNode}{$key->id} and $alnScores{$candiKey} > $rankScore - 1){
					if (!$gotFasScore and $fas_support){
						%fas_box = getFasScore();
						$gotFasScore = 1;
					}
					## get rankscore
					my $newRankScore = getFilteredRankScore($alnScores{$candiKey}, $fas_box{$candiKey});
					## candidate is better, than the last one
					if ($newRankScore > $rankScore){
						$bestTaxon = ">" . $candiKey;
						$rankScore = $newRankScore;
						($header, $seq) = getHeaderSeq($bestTaxon);
						printDebug("New Taxon has the same distance, choosing the one with higher score");
						unless ($silent) {
							print "New Best Taxon: $bestTaxon\n";
						}
					}
				}
			}
			## candidate reached the maximum score, no need to evaluate further
			if ($rankScore >= $maxScore){
				$sufficientlyClose = 1;
				printDebug("Rankscore is at maximum. Breaking loop...");
				last;
			}
			## rankscore got sufficiently close to the maximum, only evaluate candidates with the same distance now
			elsif ($rankScore >= $maxScore * (1 - $distDeviation) and !$ignoreDistance){
				printDebug("Sufficiently close to max score. Only evaluating leafs with same distance now.");
				unless ($silent) {
					print "MaxScore: $maxScore\n";
					print "RankScore: $rankScore\n";
				}
				$sufficientlyClose = 1;
			}
			clearTmpFiles();
		}
		## no candidate file was created -> so no candidate was found
		else{
			unless ($silent) {
				print "No Candidate was found for $keyName\n";
			}
		}
	}
}

my @best = (split '\|', $bestTaxon);
$currentNoRankDistNode = $newNoRankDistNode;
$currentChildsToIgnoreNode = $newChildsToIgnoreNode;
clearTmpFiles();

if ($bestTaxon ne ''){
	open (COREORTHOLOGS, ">>$outputFa") or die "Error: Could not open file: " . $outputFa . "\n";
	print COREORTHOLOGS "\n" . $header . "\n" . $seq;
	close COREORTHOLOGS;
	return $best[1];
}else{
	return '';
}
}

######################
## param: %candicontent - hashed information about candidates (id-> sequence)
## param: $scorefile - filename with alignment tool output
## cumulative alignment scores
## candidates vs sofar core ortholog set
## return: hash of scores (id->score)
sub cumulativeAlnScore{
	my $file = $_[0];
	my %content = %{$_[1]};

	my %cumscores;
	foreach my $key(keys%content) {
		my $gotScore = 0;
		open (RESULT, $file) or die "Error: Could not open file with candidate taxa\n";
		while(<RESULT>) {
			my $line = $_;
			$line =~ s/[\(\)]//g;
			my @line = split('\s+',$line);

			if($line[0] && ($line[0] eq $key)){
				if(exists $cumscores{$key}) {
					$gotScore = 1;
					$cumscores{$key} = $cumscores{$key} + $line[2];
				}else{
					$gotScore = 1;
					$cumscores{$key} = $line[2];
				}
			}
		}
		close RESULT;
		if ($gotScore == 0){
			$cumscores{$key} = 0;
		}
	}
	return %cumscores;
}

######################
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

#################################
## sorts given leaves by distance
## and delets all leaves to close to the current core
sub sort_leaves {
	my @leaves = @_;
	my $distNode = $currentChildsToIgnoreNode;
	my @candiLeaves;
	my @finalLeaves;

	for (@leaves) {
		if (!defined $hashTree{$distNode}{$_->id}){
			push @candiLeaves, $_ if defined($_->name('supplied'));
		}
	}
	while ($distNode->id != $tree->get_root_node->id and scalar @candiLeaves != 0){
		$distNode = $distNode->ancestor;
		my @nextCandiLeaves;
		my @sameDistLeaves;
		for (@candiLeaves){
			if (defined $hashTree{$distNode}{$_->id}){
				push @sameDistLeaves, $_ if defined($_->name('supplied'));
			}
			else{
				push @nextCandiLeaves, $_ if defined($_->name('supplied'));
			}
		}
		@sameDistLeaves = shuffle @sameDistLeaves;
		if (scalar @sameDistLeaves != 0){
			push @finalLeaves, \@sameDistLeaves;
		}
		@candiLeaves = @nextCandiLeaves;
	}
	return @finalLeaves;
}
####### get all taxa from the database (or the $genome_dir) where a genome is available
sub getTaxa {
	if ($dbmode) {
		my ($sql) = "select l.taxon_id, l.taxon_db, l.max_source_id, t.ncbi_id from cproteome_list.list l, taxon t where t.taxon_id = l.taxon_id and t.ncbi_id != 0";
		my ($query) = $dbHandle->prepare($sql);
		$query->execute();
		while(my @result = $query->fetchrow_array) {
			## modified by ingo: make sure to capture the max_source_id
			my $tax_src = $result[1] . '@' . $result[3] . '@' . $result[2];
			push @taxonlist, $tax_src;
			$taxa{$tax_src} = $result[3];
			printDebug("ncbiid of $tax_src is $taxa{$tax_src}");
			if ($getThemAll){
				getProteome($tax_src);
			}
		}
	}
	else {
		## removal of misplaced files in genome_dir
		if (-e "$genome_dir/query.sql"){
			unlink("$genome_dir/query.sql");
		}
		if (-e "$genome_dir/@@.fa"){
			unlink("$genome_dir/@@.fa");
		}
		@taxonlist = `ls $genome_dir`;
		chomp @taxonlist;
		for (my $i = 0; $i < @taxonlist; $i++) {
			my ($taxon_name, $ncbi_id, $src_id) = split /@/, $taxonlist[$i];
			if (!$src_id) {
				$src_id = '';
			}
			$taxon_name = $taxonlist[$i];
			$taxa{$taxon_name} = $ncbi_id;
		}
	}
	### if the blast option is chosen, we will need blast databases for all taxa
	### Baustelle: have one database including all taxa to run just a single instead of n blast searches
	if ($blast or $updateBlast_dir){
		for (my $i = 0; $i < @taxonlist; $i++){
			checkBlastDb($taxonlist[$i], $taxonlist[$i]);
		}
		if ($updateBlast_dir){
			print "\nMissing blast databases updated. Exiting.\n";
			exit;
		}
	}
	my $hashcount = keys(%taxa);
	printDebug("Returning $hashcount taxa from subroutine getTaxa");
	return(%taxa);
}
####### get all available reference taxa
sub getRefTaxa {
	@refTaxonlist = `ls $blastPath`;
	chomp @refTaxonlist;
	for (my $i = 0; $i < @refTaxonlist; $i++) {
		my ($taxon_name, $ncbi_id, $src_id) = split /@/, $refTaxonlist[$i];
		if (!$src_id) {
			$src_id = '';
		}
		$taxon_name = $refTaxonlist[$i];
		$refTaxa{$taxon_name} = $ncbi_id;
	}
	return(%refTaxa);
}
####################
sub getTree {
	# the full lineages of the species are merged into a single tree
	my $tree;
	foreach my $key (sort {lc $a cmp lc $b} keys %taxa) {
		my $node = $db->get_taxon(-taxonid => $taxa{$key});
		printDebug("\$key in sub getTree is $key and taxid is $taxa{$key}\n");
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
	if ($debug){
		print "\nTaxonomic Tree as text:\n";
		my $tree_as_string = $tree->as_text("tabtree");
		print $tree_as_string;
		print "\n";
	}
	return $tree;
}

sub getRefTree {
	# the full lineages of the species are merged into a single tree
	my $tree;
	foreach my $key (sort {lc $a cmp lc $b} keys %refTaxa) {
		my $node = $db->get_taxon(-taxonid => $refTaxa{$key});
		printDebug("\$key in sub getRefTree is $key and taxid is $refTaxa{$key}\n");
		if (!defined $node){
			print "ISSUE in sub getRefTree. No correspodence found in taxonomy file for $key and taxid $refTaxa{$key}. Skipping...\n";
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
	if ($debug){
		print "\nTaxonomic Tree as text:\n";
		my $tree_as_string = $tree->as_text("tabtree");
		print $tree_as_string;
		print "\n";
	}
	return $tree;
}

sub getTaxonName {
	my $taxAbbr = $_[0];
	my @tmp = split(/@/,$taxAbbr);
	my $taxon = $db_bkp->get_taxon($tmp[1]);
	if (defined($taxon)) {
		return($taxon->scientific_name);
	} else {
		return("Unk NCBI taxon for $taxAbbr");
	}
}

##################### perform the search for orthologs
# using the core-orthologs found in the previous steps
sub runHamstr {
	my ($taxon, $seqName, $outputFa, $refSpec, $hitlimit, $rep, $sub_strict, $subcoremode, $ev_blst, $ev_hmm, $aln) = (@_);
	my $taxaDir = $taxaPath . $taxon;
	printDebug("Running fdog: $taxon\t$seqName\t$outputFa\t$refSpec\t$taxaDir");
	if (! -e $taxaDir) {
		## backward compatibility. I used to name the dirs with the ending .dir
		if (-e "$taxaDir.dir"){
			$taxaDir = $taxaDir . '.dir';
		}
	}
	$taxaDir =~ s/\s*//g;
	if(! -e $taxaDir and $dbmode) {
		getProteome($taxon);
	}
	if (-e $taxaDir) {
		unless ($silent) {
			print "fdog for taxon: " . $taxon . "\n";
		}
		chdir($taxaDir) or die "Error: Directory for " . $taxon  . " does not exist!\n";
		my $seqfile = $taxon . ".fa";

		if(! -e $seqfile) {
			printOut("Could not find $seqfile. Check naming conventions of the files. Exiting...");
			exit;
		}

		if($seqFile ne "") {
			my $taxon_id = substr($taxon, 6, length($taxon));
			my @hamstr = ($hamstrPath, "-sequence_file=".$seqfile, "-fasta_file=".$outputFa, "-hmmpath=".$coreOrthologsPath , "-outpath=".$outputPath,
			"-blastpath=".$blastPath , "-protein", "-hmmset=".$seqName, "-taxon=".$taxon, "-force",
			"-eval_blast=".$ev_blst, "-eval_hmmer=".$ev_hmm, "-central", "-aligner=".$aln);

			my $resultFile;
			if (defined $autoLimit) {
				push(@hamstr, "-autoLimit");
			}
			elsif (defined $scoreThreshold) {
				push(@hamstr, "-scoreThreshold");
				push(@hamstr, "-scoreCutoff=$scoreCutoff");
			}
			elsif (defined $hitlimit) {
				push(@hamstr, "-hit_limit=$hitlimit");
			}
			if($sub_strict) {
				push(@hamstr, "-strict");
				$resultFile = $outputPath . "/fa_dir_" . $taxon . '_' . $seqName . "_strict/" . $seqName . ".fa";
			}
			else {
				push(@hamstr,  "-refspec=".$refSpec);
				$resultFile = $outputPath . "/fa_dir_" . $taxon . '_' . $seqName . "_" . $refSpec . "/" . $seqName . ".fa";
			}
			if($rep) {
				push(@hamstr, "-representative");
			}
			if ($checkcoorthologsref and $subcoremode==0){
				push @hamstr, '-checkCoorthologsRef';
			}
			if ($cccr and $subcoremode==1){
				push @hamstr, '-checkCoorthologsRef';
			}
			if ($rbh) {
				push @hamstr, "-rbh";
			}
			## added 2019-11-19
			if ($append) {
				push @hamstr, "-append";
			}
			##
			if ($silent) {
				push @hamstr, "-silent";
			}
			if ($debug) {
				push @hamstr, "-debug";
			}
			printDebug(@hamstr);
			system(@hamstr) == 0 or die "Error: fdog failed for " . $taxon . "\n";

			if ($outputFa !~ /extended/){
				$outputFa .= '.extended';
			}
			if(-e $resultFile) {
				unless (-e $outputFa) {
					open(EXTENDEDFA, ">$outputFa") or die "Cannot create $outputFa\n";
				} else {
					open(EXTENDEDFA, ">>$outputFa") or die "Cannot create $outputFa\n";
				}
				my $resultFa = Bio::SeqIO->new(-file => $resultFile, '-format' => 'Fasta');
				while(my $resultSeq = $resultFa->next_seq) {
					if ($resultSeq->id =~ /$taxon\|(.)+\|[01]$/) {
						my @tmpId = split("\\|", $resultSeq->id);
						print EXTENDEDFA ">$tmpId[0]\|$tmpId[-3]\|$tmpId[-2]\|$tmpId[-1]\n",$resultSeq->seq,"\n";
					}
				}
			} else {
				# add seed sequence to output extended.fa if no ortholog was found in refSpec
				if ($taxon eq $refSpec) {
					addSeedSeq($seqId, $seqName, $coreOrthologsPath, $refSpec, $outputFa);
				}
				printDebug("$resultFile not found");
			}
		}
		#remove the created folders and files
		#delete fa_dir
		my $delCommandFa;
		my $delCommandHmm;
		my $delCommandHam;
		# my $outputPathTmp = $outputPath; $outputPathTmp =~ s/\|/\\\|/g;
		# my $taxonTmp = $taxon; $taxonTmp =~ s/\|/\\\|/g;
		# my $seqNameTmp = $seqName; $seqNameTmp =~ s/\|/\\\|/g;
		if (!$strict) {
			$delCommandFa = "rm -rf  \"" . $outputPath . "/fa_dir_" . $taxon . "_" . $seqName . "_" . $refSpec . "\"";
			$delCommandHmm = "rm -rf \"" .  $outputPath . "/hmm_search_" . $taxon . "_" . $seqName . "\"";
			$delCommandHam = "rm -f \"" . $outputPath . "/hamstrsearch_" . $taxon . "_" . $seqName . ".out" . "\"";
		} else {
			$delCommandFa = "rm -rf \"" . $outputPath . "/fa_dir_" . $taxon . "_" . $seqName . "_strict" . "\"";
			$delCommandHmm = "rm -rf \"" .  $outputPath . "/hmm_search_" . $taxon . "_" . $seqName . "\"";
			$delCommandHam = "rm -f \"" . $outputPath . "/hamstrsearch_" . $taxon . "_" . $seqName . ".strict.out" . "\"";
		}
		printDebug("executing $delCommandFa", "executing $delCommandHmm", "executing $delCommandHam");
		if ($autoclean) {
			system ($delCommandFa) == 0 or die "Error deleting result files\n";
			system ($delCommandHmm) == 0 or die "Error deleting result files\n";
			system ($delCommandHam) == 0 or die "Error deleting result files\n";
		}
	}
	else {
		print "No protein set available for $taxon. Failed to fetch it from database and nothing at $taxaDir. Skipping!\n";
	}
}

# add seed sequence to output file if not exists
sub addSeedSeq {
	my ($seqId, $seqName, $coreOrthologsPath, $refSpec, $outputFa) = @_;
	unless (-e $outputFa) {
		system("touch $outputFa");
	}
	# get seed sequence and add it to the beginning of the fasta output
	open(TEMP, ">$outputFa.temp") or die "Cannot create $outputFa.temp!\n";
	my $seqio = Bio::SeqIO->new(-file => "$coreOrthologsPath/$seqName/$seqName.fa", '-format' => 'Fasta');
	while(my $seq = $seqio->next_seq) {
		my $id = $seq->id;
		if ($id =~ /$refSpec/) {
			print TEMP ">$id|1\n", $seq->seq, "\n";
			last;
		}
	}
	# then write other sequences
	my $seqio2 = Bio::SeqIO->new(-file => "$outputFa", '-format' => 'Fasta');
	while(my $seq = $seqio2->next_seq) {
		my $id = $seq->id;
		unless ($id =~ /$refSpec\|$seqId/) { # /$refSpec/) {
			print TEMP ">$id\n", $seq->seq, "\n";
		}
	}
	close(TEMP);
	system("mv $outputFa.temp $outputFa")
}

##########################
sub hmmbuild {
	# my @hmmbuild = ("hmmbuild", $_[0], $_[1]);
	# system(@hmmbuild) == 0 or die "hmmbuild failed";
	my $hmmbuild = `hmmbuild $_[0] $_[1] > /dev/null 2>&1`;
}

sub parseInput {
	my $node =  $_[0];
	my $level = $_[1];
	my $rank = $node->rank;
	printDebug("\nLEVEL:".$level."\n");
	printDebug("\nRANK:".$rank."\n");
	while($level > 0) {
		$node = $node->ancestor;
		$rank = $node->rank;
		--$level;
	}
	print "\nRETURN RANK: ".$rank."\n";
	return $rank;
}
##########################
sub parseTaxaFile {
	open (INPUT, "<$coreTaxa") or die print "Error opening file with taxa for core orthologs search\n";
	my @userTaxa;
	while(<INPUT>) {
		my $line = $_;
		chomp($line);
		if(!$taxa{$line}) {
			print "You specified " . $line . " in your core orthologs file but the taxon is not in the database!\n";
			exit;
		}
		else {
			push(@userTaxa, $line);
		}
	}
	close INPUT;
	return @userTaxa;
}
##########################
# sub printTaxa {
# 	my @result = qw();
# 	if ($dbmode) {
# 		print "taxon_schema\tsource_id\ttaxon name\n";
# 		print "------------\t---------\t----------\n";
# 		my ($sql) = "select t.name, c.taxon_db, c.max_source_id from taxon t, cproteome_list.list c where t.taxon_id=c.taxon_id";
# 		my ($query) = $dbHandle->prepare($sql);
# 		$query->execute();
# 		@result = $query->fetchrow_array;
# 		while(my @result = $query->fetchrow_array) {
# 			print $result[1] . " \t" . $result[2] . "\t" . $result[0] . "\n";
# 		}
# 	}
# 	else {
# 		print "Taxon_Name\tNCBI_ID\n";
# 		print "-------------\t------------\n";
# 		my $taxacall= "ls $genome_dir |$sedprog -e 's/@/\t/'";
# 		@result = `$taxacall`;
# 		chomp @result;
# 		print join "\n", @result;
# 		print "\n";
# 	}
# }
###########################
sub printTaxonomy {
	my $node = $_[0];
	my $i = 0;
	if($node->rank eq "species") {
		print "[" . $i . "]: " . $node->rank . " (" . $node->scientific_name . ")\n";
		while($node->ancestor) {
			$node = $node->ancestor;
			++$i;
			print "[" . $i . "]: " . $node->rank . " (" . $node->scientific_name . ")\n";
		}
	}
}
############################
sub remove_branch {
	my $node = $_[0];
	my $delFlag = 0;
	printDebug("Subroutine remove_branch\nNode is $node\nRank of node: ".$node->rank."\nNumber of leaves before removing branch ".get_leaves($tree)."\n\n");

	# undef the tree if there is only one leave left which must be removed
	if (get_leaves($tree) == 1){
		$delFlag = 1;
	}else{
		while (defined $node->ancestor) {
			last if $node->ancestor->each_Descendent > 1;
			$node = $node->ancestor;
		}
		$node->remove_all_Descendents;
		if(defined $node->ancestor) {
			$node->ancestor->remove_Descendent($node);
		}
	}
	printDebug("Subroutine remove_branch\nNode is $node\nRank of node: ".$node->rank."\nNumber of leaves after removing branch ".get_leaves($tree, $delFlag)."\n\n");
	return $delFlag;
}
############################
sub removeMaxDist {
	my $node = $tree->find_node(-ncbi_taxid => $refTaxa{$refSpec});
	my $root = $tree->get_root_node();

	if ($maxDist eq "no rank"){
		$tree->set_root_node($root);
	}else{
		while($node->rank ne $maxDist && $node != $root) {
			$node = $node->ancestor;
		}
		$tree->set_root_node($node);
	}
}
############################
# node determines the node in the tree in accordance to the given ncbi taxon id
sub removeMinDist {
	my $ncbiId = $_[0];
	my $node = $tree->find_node(-ncbi_taxid => $ncbiId);
	my $root = $tree->get_root_node();
	my $delFlag;

	printDebug("Subroutine removeMinDist\nncbiID is $ncbiId\nNode is  $node\nRank of node is ".$node->rank."\nroot is $root\nMinimal distance is $minDist\n");

	# increasing the rank of the node
	while($node->rank ne $minDist && $node != $root && defined($node->ancestor)) {
		if ($debug){
			print "Increasig the rank\nRank: ".$node->rank."\nNode: ".$node."\n\n";
		}

		$node = $node->ancestor;
	}

	#if the species has the same ranks as the references species
	if($node == $root) {
		my @toCompare = ();
		my $i = @defaultRanks - 1;
		while($i >= 0 && $defaultRanks[$i] ne $minDist) {
			push(@toCompare, $defaultRanks[$i]);
			--$i;
		}
		$node = $tree->find_node(-ncbi_taxid => $ncbiId);
		my $lastToCompare = $toCompare[$#toCompare];
		foreach(@toCompare) {
			while($node->rank eq "no rank") {
				$node = $node->ancestor;
			}
			if($node->rank ne $lastToCompare && $node->rank eq $_) {
				$node = $node->ancestor;
			}
		}
	}
	$delFlag = remove_branch($node);
	return $delFlag;
}

############################
## builds a 2 dimensional hash in which you can check for a node,
## wheather there is a path down the tree to a given species
sub buildHashTree {
	unless ($silent) {
		print "Building hash tree\n";
	}

	printDebug("Creating variables...");
	my %hashTree;
	my %nextNodes;
	my %processed;
	my @ancestors;
	my $rootNode = $wholeTree->get_root_node();

	unless ($silent ){
		print "Processing leafs...\n";
	}
	## create entry for leafes
	foreach my $leaf (get_leaves($wholeTree)){
		my $key = $leaf->id;
		my %leafHash;
		$leafHash{$key} = "exists";
		$hashTree{$leaf}{$key} = "exists";
		my $nextNode = $leaf->ancestor;
		my $nextNodeKey = $nextNode->id;
		my $test = $hashTree{$leaf}{$key};
		my $nodeTest = $nextNodeKey;
		printDebug("Leaf $key set to $test");
		## queue ancestor node for processing, if it hasn't been queued already
		if (!$nextNodes{$nextNodeKey}){
			$nextNodes{$nextNodeKey} = $nextNode;
			push @ancestors, $nextNode;
			printDebug("Queuing ancestor $nextNodeKey for processing...\n");
		}
		$processed{$leaf} = 1;
	}
	unless ($silent) {
		print "Finished leafs\n";
	}

	## create entries for all other nodes
	unless ($silent) {
		print "Processing ancestor nodes\n";
	}
	foreach my $node (@ancestors){
		my $test = $node->id;
		printDebug("Processing node: $test\n");
		my $bool = 1;
		## check, weather all childs have already been processed
		foreach  my $child ($node->each_Descendent()){
			if (!defined $processed{$child}){
				$bool = 0;
			}
		}
		## if all childs have been processed, process this node
		if ($bool == 1){
			printDebug("All children processed for node: $test");
			## node is not root
			if ($node != $rootNode){
				printDebug("Node $test is not root");
				foreach my $child ($node->each_Descendent()){
					while (my ($key, $value) = each %{$hashTree{$child}}){
						$hashTree{$node}{$key} = $value;
						printDebug("Node $key $value in node $test");
					}
				}
				my $nextNode = $node->ancestor;
				my $nextNodeKey = $nextNode->id;
				## queue ancestor node for processing, if it hasn't been queued already
				if (!$nextNodes{$nextNodeKey}){
					$nextNodes{$nextNodeKey} = $nextNode;
					push @ancestors, $nextNode;
					printDebug("Queuing ancestor $nextNodeKey for processing...");
				}
			}
			## node is root
			else{
				printDebug("Node $test is root");
				foreach my $child ($node->each_Descendent()){
					while (my ($key, $value) = each %{$hashTree{$child}}){
						$hashTree{$node}{$key} = $value;
						printDebug("Node $key $value in node $test");
					}
				}
			}
			## mark node as processed
			$processed{$node} = 1;
			printDebug("Node $test has been processed\n\n");
		}
		## not all childs have been processed
		## queue node again
		else{
			push @ancestors, $node;
			printDebug("Not all children processed for node: $test");
			printDebug("Queuing $test again...\n\n");
		}
	}
	unless ($silent) {
		print "Finished processing ancestor nodes\n";
		print "Finished building hash tree\n";
		print "Returning hash tree...\n";
	}
	return %hashTree;
}

##########################
sub getProteome {
	my $taxstring = shift;
	my $outdir = $taxstring;
	my $outfile = $taxstring;
	my @outfile;
	$taxstring =~ /(.*)@([0-9]+)@([0-9]+)/;
	my ($schema, $ncbi_id, $src_id) = ($1, $2, $3);
	print "\n\nAttempting to fetch information for $schema using source id $src_id\n\n";

	## create the relevant directory
	if (!-e "$taxaPath/$outdir"){
		print "creating directory $taxaPath/$outdir\n";
		mkpath($taxaPath."/".$outdir);
		if (-e "$taxaPath/$outdir") {
			print "succeeded\n";
		}
		else {
			print "create directory failed\n";
		}
	}

	### This is the sql statement required for fetching the sequence information from the database
	## Using Here Documents #######
	my $sql = <<"________END_OF_STATEMENT";
	use $schema;
	select concat('>',i.id, '\\n', p.seq) from ids i, protein p
	where
	i.protein_id = p.id and
	i.representative = 1 and
	i.src_id = $src_id and
	length(p.seq) > 30;
________END_OF_STATEMENT
	## the previous line must be exactly like it is to match the end note (Here Doc)

	printDebug("$sql\n");
	open (OUTQUERY, ">$taxaPath/$outdir/query.sql") or die "Could neither find nor create query.sql in $taxaPath/$outdir";
	print OUTQUERY $sql;
	close OUTQUERY;
	print "attempting to enter $taxaPath/$outdir\n";
	chdir("$taxaPath/$outdir") or die "could not enter $taxaPath/$outdir";
	`$homeDir/bin/run-query.sh $schema $ncbi_id $src_id`;
}
############
#Baustelle: run generation of BlastDb in a sub routine
## now create the relevant blast directories if necessary
sub checkBlastDb {
	my ($taxstring, $filename) = @_;
	## $taxstring identifies the species directory, $filename identifies the name of the file containing the protein set
	while (! -e "$taxaPath/$taxstring/$filename.fa"){
		my $count = 0; ## avoid endless loop
		printDebug("could not find $taxaPath/$taxstring/$filename.fa\n");
		getProteome($taxstring);
		if ($count == 5){
			die "could not find $taxaPath/$taxstring/$filename.fa and could not retrieve this information from the database.\nTerminating...\n\n";
		}
	}
	if (! -e "$blastPath/$taxstring" or $updateBlast_dir){
		`mkdir $blastPath/$taxstring`;
	}
	if (! -e "$blastPath/$taxstring/$filename.fa" or $updateBlast_dir){
		`ln -s $taxaPath/$taxstring/$filename.fa $blastPath/$taxstring/$filename.fa`;
	}
	if (! -e "$blastPath/$taxstring/$filename.pin" or $updateBlast_dir){
		chdir("$blastPath/$taxstring") or die "failed to change to dir\n";
		if ($blast_prog eq 'blastall'){
			`formatdb -i $filename.fa -t $filename -n $filename`;
		}
		elsif ($blast_prog eq 'blastp'){
			printOut("attempting to run makeblastdb", 2);
			`makeblastdb -in $filename.fa -dbtype prot -title $filename -out $filename`;
		}
	}
}
#################
sub printDebug{
	my @message = @_;
	if ($debug){
		print join "\n", @message;
		print "\n";
	}
}
sub printVariableDebug{
	my @values = @_;
	print "\n\nDEBUG\n";
	foreach (@values){
		print $_."\n";
	}
	print "\nEND OF DEBUG\n\n";
}
#################
sub getInput {
	my ($message, $dieopt) = @_;
	if ($dieopt){
		$message .= ', or type \'q\' to quit';
	}
	print ("\n" . $message . ": ");
	my $input = <STDIN>;
	chomp $input;
	if ($input =~ /^q$/i and $dieopt) {
		die "Quitting!\n";
	}
	else {
		return ($input);
	}
}
#################
sub runBlast {
	my ($query, $inpath, $outname, $outpath, $blastdb) = @_;
	printDebug("running $blast_prog on database $blastdb using input $inpath/$query and writing to $outpath/$outname.blast");

	if ($blast_prog =~ /blast[px]/) {
		!`$blast_prog -db $blastdb -seg $filter -max_target_seqs 10 -evalue $eval_blast_query -outfmt 5 -query $inpath/$query -out $outpath/$outname.blast` or die "Problem running $blast_prog\n";
	}
	elsif ($blast_prog =~ /blastall/) {
		!`$blast_prog -p $algorithm -d $blastdb -F $filter -e $eval_blast_query -m7 -i $inpath/$query -o $outpath/$outname.blast` or die "Problem running $blast_prog\n"
	}
	else {
		`$blast_prog -ublast $inpath/$query -db $blastdb -accel $accel -evalue $eval_blast_query -blast6out $outpath/$outname.blast` or die "Problem running $blast_prog\n";

		## sort the output as ublast does not do it (at least not for ESTs)
		`sort -n -r -k 12 $outpath/$outname.blast >$outpath/blastsort.tmp`;
		`mv $outpath/blastsort.tmp $outpath/$outname.blast`;
	}
	printDebug("returning $outname.blast for subroutine runBlast\n");
	return("$outname.blast");
}
#############
sub getBestBlasthit {
	my $hits;
	my $frame;
	my $count = 0;
	my ($inpath, $resultfile) = @_;
	printDebug("Sub getBestBlasthit running on $inpath/$resultfile");
	my $searchio = Bio::SearchIO->new(
	-file        => "$inpath/$resultfile",
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
			else {
				$frame = 'na';
			}
			## now I enter all top hits having the same score into the result
			$sig = $hit->score;
			if (!defined $sig_old) {
				$sig_old = $sig;
			}
			if ($sig == $sig_old) {
				$hits->[$count]->{name} = $hit->name;
				$hits->[$count]->{score} = $sig;
				$hits->[$count]->{evalue} = $hit->significance;
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
##########################
sub printOut {
	my ($message, $mlevel) = @_;
	if ($mlevel <= $vlevel){
		print "$message\n";
	}
	######################
	sub checkInt {
		my $number = shift;
		if ($number =~ /[^0-9]/){
			return();
		}
		else{
			return($number);
		}
	}
}

###########################
sub initialCheck {
	my ($seed, $ogName, $blastDir, $genomeDir, $weightDir, $fasoff) = @_;
	# check tools exist
	my @tools = ("hmmsearch", "muscle", "mafft", $globalaligner, $localaligner, $glocalaligner);
	if ($^O eq "darwin") {
		push(@tools, "clustalw2")
	} else {
		push(@tools, "clustalw")
	}
	my $flag = 1;
	foreach my $tool (@tools) {
		my $check = `which $tool`;
		if (length($check) < 1) {
			print "$tool not found\n";
			$flag = 0;
		}
	}
	if ($flag < 1) {
		die "ERROR: Some required tools not found! Please install fdog again!\n";
	}

	# check executable FAS
	my $fasCheckMsg = `setupFAS -t ./ -c 2>&1`;
	if ($fasoff != 1 && $fasCheckMsg =~ /ERROR/) {
		die "ERROR: greedyFAS not ready to use! Please check https://github.com/BIONF/FAS/wiki/prepareFAS\n";
	}

	# check seed fasta file
	unless (-e $seed) {
		$seed = "$dataDir/$seed";
	}
	my $seqio = Bio::SeqIO->new(-file => $seed, '-format' => 'Fasta');
	while(my $seq = $seqio->next_seq) {
		my $string = $seq->seq;
		if ($string =~ /[^a-zA-Z]/) {
			die "ERROR: $seed contains special characters!\n";
		}
	}

	# check ortholog group name
	if (!defined $ogName) {
		die "ERROR: Ortholog group name (-seqName) invalid!\n";
	} else {
		if ($ogName =~ /[\|\s+\"\'\`\´\!\^]/) {
			die "ERROR: Ortholog group name (-seqName) cannot contain PIPE|space or \" \' \` \´ \! \^\n";
		}
	}

	# check genome_dir
	my @genomeDir = checkValidFolderName($genomeDir);
	foreach my $genomeFd (@genomeDir) {
		unless ($genomeFd =~ /^\./) {
			my $genome = getGenomeFile("$genomeDir/$genomeFd", $genomeFd);
			unless (-e "$genome.checked") {
				die "ERROR: $genome.checked not found!\nPlease run fdog.checkData before running fdog!\n";
			}
		}
	}
	# check blast_dir
	my @blastDir = checkValidFolderName($blastDir);
	foreach my $blastFd (@blastDir) {
		unless ($blastFd =~ /^\./) {
			my $genome = getGenomeFile("$blastDir/$blastFd", $blastFd);
			unless (-e "$genome.checked") {
				die "ERROR: $genome.checked not found!\nPlease run fdog.checkData before running fdog!";
			}
		}
	}
	# check weight_dir
	if ($fasoff != 1) {
		my %seen;
		my @allTaxa = grep( !$seen{$_}++, @genomeDir, @blastDir);
		chomp(my $allAnno = `ls $weightDir | $sedprog \'s/\\.json//\'`);
		my @allAnno = split(/\n/, $allAnno);
		my @missingAnno = array_minus(@allTaxa, @allAnno);
		if (scalar @missingAnno > 0) {
			my $missingAnno = join("\n", @missingAnno);
			die "ERROR: Some taxa do not have annotation! Please turn off FAS calculation (with -fasoff), or annotate their genomes before continue.\n$missingAnno\n";
		}
	}
}

sub getGenomeFile {
	my ($folder, $filename) = @_;
	chomp(my $faFile = `ls $folder/$filename.fa* | $grepprog -v \"\\.checked\\|\\.mod\\|\\.tmp\"`);
	my $out = $faFile;
	chomp(my $link = `$readlinkprog -f $faFile`);
	if ($link ne "") {
		$out = $link;
	}
	return($out);
}

sub checkValidFolderName {
	my $folder = $_[0];
	# check if folder and its subfolders contain illegal character (e.g. pipe)
	opendir(my $dh, $folder) || die "Can't open $folder: $!";
	if ($folder =~ /[\|\s+]/) {
		die "ERROR: $folder contains illegal character (e.g. PIPE or space)!\n";
	}
	my @folders = readdir($dh);
	foreach my $fd (@folders) {
		next if ($fd eq "." or $fd eq "..");
		if ($fd =~ /[\|\s+]/) {
			die "ERROR: $folder/$fd contains illegal character (e.g. PIPE or space)!\n";
		}
	}
	closedir $dh;
	my @notFd = (".", "..");
	return(array_minus(@folders, @notFd));
}

sub gettime { sprintf"%d.%03d",Time::HiRes::gettimeofday }
sub roundtime { sprintf("%.2f", $_[0]); }

###########################
sub helpMessage {
	my $helpmessage = "
${bold}YOU ARE RUNNING $version on $hostname$norm

This program is freely distributed under a GPL.
Copyright (c) GRL limited: portions of the code are from separate copyrights

\n${bold}USAGE:${norm} oneSeq.pl -seqFile=<> -seqId=<>  -seqName=<> -refSpec=<> -minDist=<> -maxDist=<> [OPTIONS]

${bold}OPTIONS:$norm

${bold}GENERAL$norm

-h
	Invoke this help method
-version
	Print the program version

${bold}REQUIRED$norm

-seqFile=<>
	Specifies the file containing the seed sequence (protein only) in fasta format.
	If not provided the program will ask for it.
-seqId=<>
	Specifies the sequence identifier of the seed sequence in the reference protein set.
	If not provided, the program will attempt to determin it automatically.
-refSpec=<>
	Determines the reference species for the ortholog search. It should be the species the seed sequence was derived from.
	If not provided, the program will ask for it.
-minDist=<>
	specify the minimum systematic distance of primer taxa for the core set compilation.
	If not provided, the program will ask for it.
-maxDist=<>
	specify the maximum systematic distance of primer taxa to be considered for core set compilation.
	If not provided, the program will ask for it.
-coreOrth=<>
	Specify the number of orthologs added to the core set.

${bold}USING NON-DEFAULT PATHS$norm

-outpath=<>
	Specifies the path for the output directory. Default is $outputPath;
-hmmpath=<>
	Specifies the path for the core ortholog directory. Default is $coreOrthologsPath;
-blastpath=<>
	Specifies the path for the blastDB directory. Default is $blastPath;
-searchpath=<>
	Specifies the path for the search taxa directory. Default is $genome_dir;
-weightpath=<>
	Specifies the path for the pre-calculated feature annotion directory. Default is $weightPath;

${bold}ADDITIONAL OPTIONS$norm

-append
	Set this flag to append the output to existing output files
-seqName=<>
	Specifies a name for the search. If not set a random name will be set.
-db
	Run in database mode. Requires a mySql database. Only for internatl use.
-filter=[T|F]
	Switch on or off the low complexity filter for the blast search. Default: T
-silent
	Surpress output to the command line
-coreTaxa=<>
	You can provide a list of primer taxa that should exclusively be used for the compilation
	of the core ortholog set
-strict
	Run the final ortholog search in 'strict mode'. An ortholog is only then accepted when the reciprocity is fulfilled
	for each sequence in the core set.
-force
	Force the final ortholog search to create output file. Existing files will be overwritten.
-coreStrict
	Run the compilation of the core set in strict mode.
-checkCoorthologsRef
	During the final ortholog search, accept an ortholog also when its best hit in the reverse search is not the
	core ortholog itself, but a co-ortholog of it.
-CorecheckCoorthologsRef
	Invokes the 'checkCoorthologsRef' behavior in the course of the core set compilation.
-rbh
	Requires a reciprocal best hit during the ortholog search to accept a new ortholog.
-evalBlast=<>
	This option allows to set the e-value cut-off for the Blast search. Default: 1E-5
-evalHmmer=<>
	This options allows to set the e-value cut-off for the HMM search. Default: 1E-5
-evalRelaxfac=<>
	This options allows to set the factor to relax the e-value cut-off (Blast search and HMM search) for the final ortholog run. Default: 10
-hitLimit=<>
	Provide an integer specifying the number of hits of the initial pHMM based search that should be evaluated
	via a reverse search. Default: 10
-coreHitLimit=<>
	Provide an integer specifying the number of hits of the initial pHMM based search that should be evaluated
	via a reverse search. Default: 3
-autoLimit
	Setting this flag will invoke a lagPhase analysis on the score distribution from the hmmer search. This will determine automatically
	a hit limit for each query. Note, when setting this flag, it will be effective for both the core ortholog compilation
	and the final ortholog search.
-scoreThreshold
	Instead of setting an automatic hit limit, you can specify with this flag that only candidates with an hmm score no less
	than x percent of the hmm score of the best hit are further evaluated. Default is x = 10.
	You can change this cutoff with the option -scoreCutoff. Note, when setting this flag, it will be effective for
	both the core ortholog compilation and the final ortholog search.
-scoreCutoff=<>
	In combination with -scoreThreshold you can define the percent range of the hmms core of the best hit up to which a
	candidate of the hmmsearch will be subjected for further evaluation. Default: 10%.
-coreOnly
	Set this flag to compile only the core orthologs. These sets can later be used for a stand alone ortholog search.
-reuseCore
	Set this flag if the core set for your sequence is already existing. No check currently implemented.
-ignoreDistance
	Set this flag to ignore the distance between Taxa and to choose orthologs only based on score
-distDeviation=<>
	Specify the deviation in score in percent (1=100%, 0=0%) allowed for two taxa to be considered similar
-blast
	Set this flag to determine sequence id and refspec automatically. Note, the chosen sequence id and reference species
	does not necessarily reflect the species the sequence was derived from.
-rep
	Set this flag to obtain only the sequence being most similar to the corresponding sequence in the core set rather
	than all putative co-orthologs.
-coreRep
	Set this flag to invoke the '-rep' behaviour for the core ortholog compilation.
-cpu
	Determine the number of threads to be run in parallel
-hyperthread
	Set this flag to use hyper threading
-batch=<>
	Currently has NO functionality.
-group=<>
	Allows to limit the search to a certain systematic group
-cleanup
	Temporary output will be deleted.
-aligner
	Choose between mafft-linsi or muscle for the multiple sequence alignment. DEFAULT: muscle
-local
	Specify the alignment strategy during core ortholog compilation. Default is local.
-glocal
	Set the alignment strategy during core ortholog compilation to glocal.
-searchTaxa
	Input file containing list of search taxa.
${bold}SPECIFYING FAS SUPPORT OPTIONS$norm

-fasoff
	Turn OFF FAS support. Default is ON.
-coreFilter=[relaxed|strict]
	Specifiy mode for filtering core orthologs by FAS score. In 'relaxed' mode candidates with insufficient FAS score will be disadvantaged.
	In 'strict' mode candidates with insufficient FAS score will be deleted from the candidates list. Default is None.
	The option '-minScore=<>' specifies the cut-off of the FAS score.
-minScore=<>
	Specify the threshold for coreFilter. Default is 0.75.
-countercheck
	Set this flag to counter-check your final profile. The FAS score will be computed in two ways (seed vs. hit and hit vs. seed).

${bold}SPECIFYING EXTENT OF OUTPUT TO SCREEN$norm

-debug
	Set this flag to obtain more detailed information about the programs actions
-silent
	Surpress output to screen as much as possbile
\n\n";
	return($helpmessage);
}
