#!/bin/bash

sys="$(uname)" # Linux for Linux or Darwin for MacOS
echo "Current OS system: $sys"

CURRENT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
CURRENT="${CURRENT/\/setup/}"
BIN="$CURRENT/bin"

flag=0
fas=1
installLib=0
homedir="$(echo $HOME)"
outDir=$CURRENT

while getopts lfo: opt; do
  case ${opt} in
    o )
    echo "Data output path: $OPTARG"
    outDir=$OPTARG
    ;;
    l )
    echo "INSTALL LIB"
    installLib=1
    ;;
    f )
    echo "NO FAS!"
    fas=0
    ;;
    \? )
    echo "Usage: setup.sh [-l] [-f]"
    exit 1
    ;;
  esac
done

### install dependencies
if [ $installLib == 1 ]; then
    if [ "$sys" == "Darwin" ]; then
      $CURRENT/setup/install_lib.sh
    else
      echo "Enter sudo password to install required libraries..."
      sudo $CURRENT/setup/install_lib.sh
    fi
    exit
fi

### check grep, sed, readlink and wget availability
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

if [ -z "$(which $sedprog)" ]; then
  echo -e "\e[31m$sedprog not found!\e[0m"
  echo "Please run fdog.setup with --lib first!"
  exit
fi

if [ -z "$(which $grepprog)" ]; then
  echo -e "\e[31m$grepprog not found!\e[0m"
  echo "Please run fdog.setup with --lib first!"
  exit
fi

if [ -z "$(which $wgetprog)" ]; then
  echo -e "\e[31m$wgetprog not found!\e[0m"
  echo "Please run fdog.setup with --lib first!"
  exit
fi

if [ -z "$(which $readlinkprog)" ]; then
  echo -e "\e[31m$readlinkprog not found!\e[0m"
  echo "Please run fdog.setup with --lib first!"
  exit
fi

if ! [ -f ~/$bashFile ]; then
  touch ~/$bashFile
fi
if ! [ -f ~/$rprofile ]; then
  touch ~/$rprofile
fi
echo "done!"

### prepare folders
echo "-------------------------------------"
echo "Preparing folders..."
if [ ! -d "$CURRENT/taxonomy" ]; then mkdir "$CURRENT/taxonomy"; fi
if [ ! -d "$CURRENT/bin" ]; then mkdir "$CURRENT/bin"; fi
if [ ! -d "$CURRENT/bin/aligner" ]; then mkdir "$CURRENT/bin/aligner"; fi
echo "done!"

### download tools
echo "-------------------------------------"
echo "Downloading and installing annotation tools/databases:"

fasta36="yes"
if [ -z "$(which fasta36)" ]; then
  fasta36="no"
  fasta36v="fasta-36.3.8h"
  if ! [ -f "bin/aligner/bin/fasta36" ]; then
    echo "fasta-36"
    wget "http://faculty.virginia.edu/wrpearson/fasta/fasta36/${fasta36v}.tar.gz"
    tar xf $fasta36v.tar.gz
    rm "${fasta36v}.tar.gz"
    mv $fasta36v/* $CURRENT/bin/aligner/
    rm -rf $fasta36v
    cd "$CURRENT/bin/aligner/src"
    if [ $sys=="Linux" ]; then
      make -f ../make/Makefile.linux64_sse2 all
    elif [ $sys=="Darwin" ]; then
      make -f ../make/Makefile.os_x86_64 all
    fi
  fi
  if [ -z "$($grepprog PATH=$CURRENT/bin/aligner/bin ~/$bashFile)" ]; then
    echo "export PATH=$CURRENT/bin/aligner/bin:\$PATH" >> ~/$bashFile
  fi
fi
cd $CURRENT
if [ -z "$(which fasta36)" ]; then
  if ! [ -f "$CURRENT/bin/aligner/bin/fasta36" ]; then
    echo -e "\e[31mfasta36 tool could not be found in $CURRENT/bin/aligner/. Please check again!\e[0m"
    exit
  fi
fi

cd "$CURRENT/taxonomy"
if ! [ -f "nodes" ]; then
  wget "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
  tar xf taxdump.tar.gz
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

fasPrepare=0
if [ $fas == 1 ]; then
  cd "$CURRENT/bin"
  if [ -z "$(which annoFAS)" ]; then
    echo "FAS"
    pip install --user greedyFAS
    if [ -z "$($grepprog \$HOME/.local/bin:\$PATH ~/$bashFile)" ]; then
      echo "export PATH=\$HOME/.local/bin:\$PATH" >> ~/$bashFile
    fi
    if [ -z "$($grepprog $homedir/.local/bin ~/$rprofile)" ]; then
      echo "Sys.setenv(PATH = paste(\"$homedir/.local/bin\", Sys.getenv(\"PATH\"), sep=\":\"))" >> ~/$rprofile
    fi
    fasPrepare=1
  else
    if ! [ -z "$(prepareFAS -t ./ --check 2>&1 | grep ERROR)" ]; then
      fasPrepare=1
    fi
  fi

  cd $CURRENT
  source ~/$bashFile
  if [ -z "$(which annoFAS)" ]; then
    echo -e "Installation of FAS failed! Please try again or install FAS by yourself using \e[91mpip install greedyFAS\e[0m!"
    echo -e "For more info, please check FAS website at \e[91mhttps://github.com/BIONF/FAS\e[0m"
    exit
  else
    if ! [ -z "$(prepareFAS -t ./ --check 2>&1 | grep ERROR)" ]; then
      fasPrepare=1
    fi
  fi
  echo "done!"
fi

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
    if [ -d "$outDir/genome_dir" ]; then
      for i in $(ls "$outDir/genome_dir"); do rm -f "$outDir/genome_dir/$i/$i.fa.mod"; done
    fi

    if [ "$(ls -A $outDir/blast_dir)" ]; then
      echo "Data should be in place to run fdog."
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

echo "done!"

### final check
echo "-------------------------------------"
echo "Final check..."
flag=0

echo "Tools"
dependencies=(
genewise
hmmsearch
hmmscan
hmmbuild
mafft
muscle
clustalw
blastp
)

for i in "${dependencies[@]}"; do
  tool=$i
  if [ $tool == "clustalw" ]; then
    if [ "$sys" == "Darwin" ]; then
      tool="clustalw2"
    fi
  fi
  if [ -z "$(which $tool)" ]; then
    echo -e "\t\e[31mWARNING $tool not found!\e[0m"
    flag=1
  fi
done

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

echo "Perl modules"
for i in "${perlModules[@]}"; do
  msg=$((perl -e "use $i") 2>&1)
  if ! [[ -z ${msg} ]]; then
    echo -e "\t\e[31mWARNING $i could not be installed\e[0m"
    flag=1
  fi
done

echo "Environment paths"
envPaths=(
WISECONFIGDIR
)
for i in "${envPaths[@]}"; do
  if [ -z "$($grepprog $i ~/$bashFile)" ]; then
    echo -e "\t\e[31mWARNING $i was not added into ~/$bashFile\e[0m"
    flag=1
  fi
done
if [ "$fasta36" == "no" ]; then
  if [ -z "$($grepprog PATH=$CURRENT/bin/aligner/bin ~/$bashFile)" ]; then
    echo -e "\t\e[31mWARNING $CURRENT/bin/aligner/bin was not added into ~/$bashFile\e[0m"
    flag=1
  fi
fi
echo "done!"

if [ "$flag" == 1 ]; then
  echo "Some tools/libraries counld not installed correctly or paths were not added into ~/$bashFile."
  echo "Please manually install the missing dependencies using fdog.setup with --lib option (ask your admin if you don't have root privileges)."
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
    echo -e "\e[91mPLEASE RUN\e[0m \e[96mprepareFAS\e[0m \e[91mTO CONFIGURE FAS BEFORE USING fdog!\e[0m"
    echo "Then you can test fdog with:"
  else
    echo "All tests succeeded, fdog should be ready to run. You can test it with:"
  fi
  echo -e "\e[96mfdog.run --seqFile infile.fa --seqName test --refspec HUMAN@9606@3\e[0m"
  echo "Output files with prefix \"test\" will be found at your current working directory!"
  echo -e "For more details, use \e[96mfdog.run -h\e[0m or visit https://github.com/BIONF/fDOG/wiki"
  echo "Happy running fdog! ;-)"
fi
exit 1
