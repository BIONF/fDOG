#!/bin/bash

sys="$(uname)" # Linux for Linux or Darwin for MacOS

flag=0

### update GPG key (Google signature key for signing and authenticating packages)
if ! [ "$sys" == "Darwin" ]; then
  wget -q -O - https://dl.google.com/linux/linux_signing_key.pub | sudo apt-key add -
fi

### check grep, sed and wget availability
grepprog='grep'
sedprog='sed'
readlinkprog='readlink'
wgetprog='wget'
bashFile='.bashrc'
if [ "$sys" == "Darwin" ]; then
  if [ -z "$(which brew)" ]; then
    echo "Please install homebrew to install dependencies tools and libraries!"
    echo "Check https://brew.sh"
    exit
  fi
  sedprog='gsed'
  grepprog='ggrep'
  readlinkprog='greadlink'
  shell=$(echo $SHELL)
  if [ $shell == "/bin/zsh" ]; then
    bashFile='.zshrc'
  else
    bashFile='.bash_profile'
  fi
else
  if [ "$EUID" -ne 0 ]; then
    echo "You must run this setup as a root user!"
    exit
  fi
fi

if [ -z "$(which $readlinkprog)" ]; then
  if [ "$sys" == "Darwin" ]; then
    brew install coreutils
  fi
fi

if [ -z "$(which $sedprog)" ]; then
  if [ "$sys" == "Darwin" ]; then
    brew install gnu-sed
  fi
fi

if [ -z "$(which $grepprog)" ]; then
  if [ "$sys" == "Darwin" ]; then
    brew install grep
  fi
fi

if [ -z "$(which $wgetprog)" ]; then
  if [ "$sys" == "Darwin" ]; then
    brew install wget
  fi
fi

if ! [ -f ~/$bashFile ]; then
  touch ~/$bashFile
fi

### check dependencies
echo "-------------------------------------"
echo "Installing dependencies..."

dependenciesUbuntu=(
  build-essential # for make
  curl
  r-base # for Statistics::R
  wise
  hmmer # hmmer (for both hmmsearch and hmmbuild)
  clustalw
  mafft
  muscle
  blast2 # blast
  ncbi-blast+
  libdbi-perl
  libipc-run-perl
  perl-doc
  locales
  lib32z1
)

dependenciesMac=(
  brewsci/bio/genewise
  hmmer # hmmer (for both hmmsearch and hmmbuild)
  brewsci/bio/clustal-w
  mafft
  brewsci/bio/muscle
  blast
)

if [ "$sys" == "Darwin" ]; then
  for i in "${dependenciesMac[@]}"; do
    echo $i
    brew install $i
  done
  if [ -z "$(grep clustalw ~/$bashFile)" ]; then
    echo "alias clustalw='clustalw2'" >> ~/$bashFile
  fi
else
  sudo apt-get update -y
  for i in "${dependenciesUbuntu[@]}"; do
    echo $i
    sudo apt-get install -y -qq $i > /dev/null
  done
fi

dependencies=(
  genewise
  hmmsearch
  hmmbuild
  mafft
  muscle
  blastn
)

for i in "${dependencies[@]}"; do
  if [ -z "$(which $i)" ]; then
    echo "$i not found / cannot be automatically installed. Please install it manually and run this setup again!"
    flag=1
  fi
done
if [ "$flag" == 1 ]; then exit 1; fi

wisePath=$(which "genewise")
if [ -z "$(grep WISECONFIGDIR=$wisePath ~/$bashFile)" ]; then
  echo "export WISECONFIGDIR=${wisePath}" >> ~/$bashFile
fi

echo "Installing Perl modules..."
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

if [ -z "$(which cpanm)" ]; then
  curl -L http://cpanmin.us | perl - --sudo App::cpanminus
fi
for i in "${perlModules[@]}"; do
  msg=$((perldoc -l $i) 2>&1)
  if [[ "$(echo $msg)" == *"No documentation"* ]]; then
    sudo cpanm ${i} --quiet --force
  fi
done
echo "done!"

echo "-------------------------------------"
CURRENT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "Please run fdog.setup without --lib option to continue setup fdog!"
