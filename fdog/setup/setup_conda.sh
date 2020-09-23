#!/bin/bash

sys="$(uname)" # Linux for Linux or Darwin for MacOS
echo "Current OS system: $sys"

CURRENT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
CURRENT="${CURRENT/\/setup/}"
BIN="$CURRENT/bin"

flag=0
outDir=$CURRENT

while getopts o: opt; do
  case ${opt} in
    o )
    echo "Data output path: $OPTARG"
    outDir=$OPTARG
    ;;
    \? )
    echo "Usage: setup.sh [-l] [-f]"
    exit 1
    ;;
  esac
done
if [ ! -d "$outDir" ]; then mkdir "$outDir"; fi

### check grep, sed and wget availability
echo "-------------------------------------"
echo "Checking .bash_profile/.bashrc, grep, sed/gsed and wget availability..."
grepprog='grep'
sedprog='sed'
readlinkprog='readlink'
wgetprog='wget'
bashFile='.bashrc'
rprofile='.Rprofile'

if [ "$sys" == "Darwin" ]; then
  sedprog='gsed'
  grepprog='ggrep'
  readlinkprog='greadlink'
  shell=$(echo $SHELL)
  if [ $shell == "/bin/zsh" ]; then
    bashFile='.zshrc'
  else
    bashFile='.bash_profile'
  fi
fi

# NOTE: install only available for Linux!
if [ -z "$(which $sedprog)" ]; then
  if [ "$sys" == "Darwin" ]; then
    echo -e "\e[31m$sedprog not found. Please install it first (e.g. using brew)!\e[0m"
    flag=1
  fi
  conda install -c conda-forge sed
fi

if [ -z "$(which $grepprog)" ]; then
  if [ "$sys" == "Darwin" ]; then
    echo -e "\e[31m$grepprog not found. Please install it first (e.g. using brew)!\e[0m"
    flag=1
  fi
  conda install -c bioconda grep
fi

if [ -z "$(which $wgetprog)" ]; then
  if [ "$sys" == "Darwin" ]; then
    echo -e "\e[31m$wgetprog not found. Please install it first (e.g. using brew)!\e[0m"
    flag=1
  fi
  conda install -c anaconda wget
fi

if [ -z "$(which $readlinkprog)" ]; then
  if [ "$sys" == "Darwin" ]; then
    echo -e "\e[31m$readlinkprog not found. Please install it first (e.g. using brew)!\e[0m"
    flag=1
  fi
fi

if ! [ -f ~/$bashFile ]; then
  touch ~/$bashFile
fi
if ! [ -f ~/$rprofile ]; then
  touch ~/$rprofile
fi
if [ "$flag" == 1 ]; then exit 1; fi
echo "done!"

### check dependencies
echo "-------------------------------------"
echo "Installing dependencies..."

if [ -z "$(which R)" ]; then
  echo "R"
  conda install -y r
fi

if [[ -z $(conda list | $grepprog "pkg-config ") ]]; then
  echo "pkg-config"
  conda install -y pkg-config
fi

if [[ -z $(conda list | $grepprog "perl-bioperl ") ]]; then
  echo "perl-bioperl"
  conda install -y -c bioconda perl-bioperl
  conda install -y -c bioconda perl-bioperl-core
  conda install -y -c bioconda perl-bioperl-run
fi

dependencies=(
  blastp # blast
  genewise # wise2
  hmmsearch # hmmer (for both hmmsearch and hmmbuild)
  clustalw
  mafft # for linsi
  muscle
  fasta36
)

for i in "${dependencies[@]}"; do
  if [ -z "$(which $i)" ]; then
    echo $i
    tool=$i
    if [ "$tool" = "blastp" ]; then
      conda install -y -c bioconda blast
    elif [ "$tool" = "hmmsearch" ]; then
      conda install -y -c bioconda hmmer
    elif [ "$tool" = "genewise" ]; then
      conda install -y -c bioconda wise2
      wisePath=$(which "genewise")
      if [ -z "$($grepprog WISECONFIGDIR=$wisePath ~/$bashFile)" ]; then
        echo "export WISECONFIGDIR=${wisePath}" >> ~/$bashFile
      fi
    elif [ "$tool" = "fasta36" ]; then
      conda install -y -c bioconda fasta3
    else
      conda install -y -c bioconda $i
    fi
  fi
done

for i in "${dependencies[@]}"; do
  if [ -z "$(which $i)" ]; then
    echo -e "\e[31m$i not found. Please install it to use fdog!\e[0m"
    flag=1
  fi
done
if [ "$flag" == 1 ]; then exit 1; fi

perlModules=(
  Array::Utils
  Capture::Tiny
  DBI
  DB_File
  File::Copy
  File::Path
  File::Basename
  File::Which
  List::Util
  Parallel::ForkManager
  POSIX
  XML::SAX
  XML::NamespaceSupport
  XML::Parser
  Getopt::Long
  IO::Handle
  IPC::Run
  Statistics::R
  Term::Cap
  Time::HiRes
  Bio::AlignIO
  Bio::Align::ProteinStatistics
  Bio::DB::Taxonomy
  Bio::SearchIO
  Bio::SearchIO::blastxml
  Bio::Search::Hit::BlastHit
  Bio::Seq
  Bio::SeqIO
  Bio::SeqUtils
  Bio::Tree::Tree
  Bio::Tools::Run::StandAloneBlast
)

for i in "${perlModules[@]}"; do
  msg=$((perldoc -l $i) 2>&1)
  if [[ "$(echo $msg)" == *"No documentation"* ]]; then
    cpanm ${i} --quiet --force
  fi
done

echo "done!"

### prepare folders
echo "-------------------------------------"
echo "Preparing folders..."

# create required folders
if [ ! -d "$CURRENT/taxonomy" ]; then mkdir "$CURRENT/taxonomy"; fi
if [ ! -d "$CURRENT/bin" ]; then mkdir "$CURRENT/bin"; fi
if [ ! -d "$CURRENT/bin/aligner" ]; then mkdir "$CURRENT/bin/aligner"; fi
echo "done!"

### download tools
echo "-------------------------------------"
echo "Downloading and installing annotation tools/databases:"

cd "$CURRENT/taxonomy"
if ! [ -f "nodes" ]; then
  wget "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
  tar xfv taxdump.tar.gz
  rm taxdump.tar.gz
  echo "Taxonomy database indexing. It can take a while, please wait..."
  perl $CURRENT/setup/indexTaxonomy.pl $CURRENT/taxonomy
  rm citations.dmp
  rm delnodes.dmp
  rm division.dmp
  rm gencode.dmp
  rm merged.dmp
  rm gc.prt
  rm readme.txt
fi
cd $CURRENT
if ! [ -f "$CURRENT/taxonomy/nodes" ]; then
  echo -e "\e[31mError while indexing NCBI taxonomy database! Please check $CURRENT/taxonomy/ folder and run this setup again!\e[0m"
  exit
fi

cd "$CURRENT/bin"
fasPrepare=0
if [ -z "$(which annoFAS)" ]; then
  echo "FAS"
  conda install -y -c BIONF fas
  if [ -z "$(which annoFAS)" ]; then
    echo -e "\e[31mInstallation of FAS failed! Please try again!\e[0m"
    exit
  fi
  fasPrepare=1
else
  if ! [ -z "$(prepareFAS -t ./ --check 2>&1 | grep ERROR)" ]; then
    fasPrepare=1
  fi
fi

if [ -z "$(which annoFAS)" ]; then
  echo -e "Installation of FAS failed! Please try again or install FAS by yourself using \e[91mconda install -c BIONF fas\e[0m or \e[91mpip install greedyFAS\e[0m"
  echo -e "For more info, please check FAS website at \e[91mhttps://github.com/BIONF/FAS\e[0m"
  exit
else
  if ! [ -z "$(prepareFAS -t ./ --check 2>&1 | grep ERROR)" ]; then
    fasPrepare=1
  fi
fi
cd $CURRENT
echo "done!"

### download data
data_fdog_file="data_HaMStR-2019c.tar.gz"
checkSumData="1748371655 621731824 $data_fdog_file"
cd $outDir
if [ ! -d "$outDir/core_orthologs" ]; then mkdir "$outDir/core_orthologs"; fi

if ! [ "$(ls -A $outDir/genome_dir)" ]; then
  echo "-------------------------------------"
  echo "Getting pre-calculated data"

  echo "Processing $outDir ..."
  if [ ! -f $outDir/$data_fdog_file ]; then
    echo "Downloading data from https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/$data_fdog_file"
    wget --no-check-certificate https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/$data_fdog_file
  else
    CHECKSUM=$(cksum $data_fdog_file)
    echo "Checksum: $CHECKSUM"
    if ! [ "$CHECKSUM" == "$checkSumData" ]; then
      rm $outDir/$data_fdog_file
      echo "Downloading data from https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/$data_fdog_file"
      wget --no-check-certificate https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/$data_fdog_file
    fi
  fi

  if [ ! -f $outDir/$data_fdog_file ]; then
    echo "File $data_fdog_file not found! Please try to download again from"
    echo "https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/$data_fdog_file"
    exit
  fi

  CHECKSUM=$(cksum $data_fdog_file)
  if [ "$CHECKSUM" == "$checkSumData" ]; then
    echo "Extracting archive $data_fdog_file..."
    tar xf $outDir/$data_fdog_file
    rm $outDir/$data_fdog_file
    for i in $(ls "$outDir/genome_dir"); do rm -f "$outDir/genome_dir/$i/$i.fa.mod"; done

    if [ "$(ls -A $outDir/blast_dir)" ]; then
      echo "Data should be in place to run fdog.\n"
    else
      echo -e "\e[31mSomething went wrong with the download. Data folders are empty.\e[0m"
      echo "Please try to download again from"
      echo "https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/$data_fdog_file"
      echo "Or contact us if you think this is our issue!"
      exit
    fi
  else
    echo -e "\e[31mSomething went wrong with the download. Checksum does not match.\e[0m"
    echo "Please try to download again from"
    echo "https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo/$data_fdog_file"
    echo "Please put it into $outDir folder and run this setup again!"
    exit
  fi
fi
# write data path to pathConfig file
if [ -f $BIN/pathconfig.txt ]; then
  rm $BIN/pathconfig.txt
fi
touch $BIN/pathconfig.txt
echo $outDir >> $BIN/pathconfig.txt

### add paths to bash profile file
echo "-------------------------------------"
echo "Adding WISECONFIGDIR to ~/$bashFile"

wisePath=$(which "genewise")
if [ -z "$($grepprog WISECONFIGDIR=$wisePath ~/$bashFile)" ]; then
  echo "export WISECONFIGDIR=${wisePath}" >> ~/$bashFile
fi

# echo "Adding paths to ~/$rprofile"
# if [ -z "$($grepprog $CURRENT/bin ~/$rprofile)" ]; then
#   echo "Sys.setenv(PATH = paste(\"$CURRENT/bin\", Sys.getenv(\"PATH\"), sep=\":\"))" >> ~/$rprofile
# fi

echo "done!"

### adapt paths in fdog scripts
echo "-------------------------------------"
echo "Adapting paths in fdog scripts"
# update the sed and grep commands
$sedprog -i -e "s/\(my \$sedprog = '\).*/\1$sedprog';/" $CURRENT/bin/hamstr.pl
$sedprog -i -e "s/\(my \$grepprog = '\).*/\1$grepprog';/" $CURRENT/bin/hamstr.pl
$sedprog -i -e "s/\(my \$readlinkprog = '\).*/\1$readlinkprog';/" $CURRENT/bin/hamstr.pl
$sedprog -i -e "s/\(my \$sedprog = '\).*/\1$sedprog';/" $CURRENT/bin/oneSeq.pl
$sedprog -i -e "s/\(my \$grepprog = '\).*/\1$grepprog';/" $CURRENT/bin/oneSeq.pl
$sedprog -i -e "s/\(my \$readlinkprog = '\).*/\1$readlinkprog';/" $CURRENT/bin/oneSeq.pl

# localize the perl installation
path2perl=`which perl`
echo "path to perl: $path2perl"
$sedprog -i -e "s|\#\!.*|\#\!$path2perl|g" $CURRENT/bin/hamstr.pl
$sedprog -i -e "s|\#\!.*|\#\!$path2perl|g" $CURRENT/bin/translate.pl
$sedprog -i -e "s|\#\!.*|\#\!$path2perl|g" $CURRENT/bin/oneSeq.pl

### final check
echo "-------------------------------------"
echo "Final check..."
flag=0

echo "Conda packages"
condaPkgs=(
perl-bioperl
perl-bioperl-core
blast
hmmer
wise2
clustalw
mafft
muscle
fasta3
)
for i in "${condaPkgs[@]}"; do
  if [[ -z $(conda list | $grepprog "$i ") ]]; then
    progname=$i
    if [ "$i" == "blast" ]; then
      progname="blastp"
    elif [ "$i" == "wise2" ]; then
      progname="genewise"
    elif [ "$i" == "hmmer" ]; then
      progname="hmmsearch"
    elif [ "$i" == "fasta3" ]; then
      progname="fasta36"
    fi
    if [ -z "$(which $progname)" ]; then
      echo -e "\t\e[31m$i could not be installed\e[0m"
      flag=1
    fi
  fi
done
echo "done!"

echo "Perl modules"
for i in "${perlModules[@]}"; do
  msg=$((perl -e "use $i") 2>&1)
  if ! [[ -z ${msg} ]]; then
    echo -e "\t\e[31m$i could not be installed\e[0m"
    flag=1
  fi
done
echo "done!"

echo "Environment paths"
envPaths=(
WISECONFIGDIR
)
for i in "${envPaths[@]}"; do
  if [ -z "$($grepprog $i ~/$bashFile)" ]; then
    echo -e "\t\e[31m$i was not added into ~/$bashFile\e[0m"
    flag=1
  fi
done
echo "done!"

if [ "$flag" == 1 ]; then
  echo "Some tools/libraries counld not installed correctly or paths were not added into ~/$bashFile."
  echo "Please manually install the missing dependencies using using fdog.setup with --lib option (ask your admin if you don't have root privileges)."
  echo "Then run this setup again to try one more time!"
  exit
else
  echo "Generating symbolic links"
  ln -s -f $CURRENT/bin/hamstr.pl $CURRENT/bin/hamstr
  ln -s -f $CURRENT/bin/oneSeq.pl $CURRENT/bin/oneSeq
  echo "Sourcing bash profile file"
  source ~/$bashFile
  echo "-------------------------------------"
  $sedprog -i -e 's/my $configure = .*/my $configure = 1;/' $CURRENT/bin/hamstr.pl
  $sedprog -i -e 's/my $configure = .*/my $configure = 1;/' $CURRENT/bin/oneSeq.pl
  if [ "$fasPrepare" == 1 ]; then
    echo "All tests succeeded."
    echo -e "\e[91mPLEASE RUN\e[0m \e[96msetupFAS\e[0m \e[91mTO CONFIGURE FAS BEFORE USING fdog!\e[0m"
    echo "Then you can test fdog with:"
  else
    echo "All tests succeeded, fdog should be ready to run. You can test it with:"
  fi
  echo -e "\e[96mfdog.run --seqFile infile.fa --seqName test --refspec HUMAN@9606@3\e[0m"
  echo "Output files with prefix \"test\" will be found at your current working directory!"
  echo -e "For more details, use \e[96mfdog.run -h\e[0m or visit https://github.com/BIONF/fDOG/wiki"
  echo "Happy using fdog! ;-)"
fi
exit 1
