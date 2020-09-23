package run_genewise_hamstr;
use strict;
#$ENV{'WISECONFIGDIR'} =  "/usr/local/src/wise2.2.0/wisecf/";
# this module runs genewise on a DNA sequence and a protein sequence
# and then allows to parse this result.
# the constructor creates an object containing a reference to an array
# containing the file content

# Modified 11.01.2010 renamed the file names for the genewise run to avoid overwriting of files when multipe runs are performed in parallel on the same sequence file
# LAST Modified: 31.07.2015. Added the option to keep, mask or remove partial codons and introns from
# the transcript. 

1;
sub new {
    my $self_tmp = [];
    my $self;
    my ($class, $dna, $prot, $path, $keepintron) = @_;
    if (!defined $path) {
	$path = '/tmp';
    }
    if (!defined $keepintron) {
	$keepintron = 2;
    }
    my $pid=$$;
    # the file names
    my $protname = $pid.'_protein';
    my $dnaname = $pid . '_dna';
    ## print the two sequences to default path /tmp/
    open (DNA, ">$path/$dnaname") or die "could not open $path/$dnaname for writing\n";
    print DNA ">$dnaname\n$dna";
    close DNA;
    open (PROTEIN, ">$path/$protname") or die "could not open $path/$protname for writing\n";
    print PROTEIN ">$protname\n$prot";
    close PROTEIN;

    ## run genewise on the two sequences
  `echo \$WISECONFIGDIR`;
    
    $self_tmp = [`genewise -trans -cdna -pep -sum $path/$protname $path/$dnaname`];
    for (my $i = 0; $i < @$self_tmp; $i++) {
	$self_tmp->[$i] =~ s/\s{1,}$//;
    }
    $self->{gw} = $self_tmp;
    $self->{nt_seq} = $dna;
    $self->{prot_seq} = $prot;
    $self->{protname} = $protname;
    $self->{dnaname} = $dnaname;
    $self->{gw_count} = @$self_tmp;

    if ($keepintron =~ /^k/i ) {
    	$self->{get_indel} = 2; ## per default the indel-part is recovererd in lower case letters rather than masked or removed. See code for details
    }
    elsif ($keepintron =~ /^m/i) {
	$self->{get_indel} = 1; ## The indel-part is masked. See code for details;
   }
   else {
	$self->{get_indel} = 0; ## the indel-part is removed making the cDNA consistent with the translaton. See code for details;
   }
   print "intron is $self->{get_indel}\n";

    $self->{indels} = _GetIndels($self_tmp);
    bless ($self, $class);
    return $self;}
#################
## sub score extract the score for the alignment
sub score {
    my $self = shift;
    my $score;
    for (my $i = 0; $i < $self->{gw_count}; $i ++) {
	if ($self->{gw}->[$i] =~ /^(\d{1,}\.{0,1}\d{0,}).*/) {
	    $score = $1;
	    last;
	}
    }
    return ($score);
}
##################
sub protein {
    my $self = shift;
    my $gw = $self->{gw};
    my $prot = '';
    for (my $i = 0; $i < @$gw; $i++) {
      if ($gw->[$i] =~ />.*\.pep/) { #the protein seq starts
	my $count = 1;
	while ($gw->[$i+$count] ne '//') {
	  my $protpart = $gw->[$i+$count];
	  chomp $protpart;
	  $prot .= $protpart;
	  $count ++;
	}
      }
      elsif (length $prot > 0) {
	last;
      }
    }
    return($prot);
 }
##################
sub translation {
    my $self = shift;
    my $finish = 0;
    my $translated_seq = '';
    my @transtmp;

    ## step 1: extract the relevant info from the genewise output
    for (my $i = 0; $i < $self->{gw_count}; $i++) {
      if ($self->{gw}->[$i] =~ />.*.tr/) {# a translated bit starts
	while ($self->{gw}->[$i] !~ '//') {
	  push @transtmp, $self->{gw}->[$i];
	  $i++;
	}
	last; # end the for loop since nothing left to be done
      }
    }
    
    ## step two: get the sequences
    my $count = -1;
    my $trans;
    for (my $i = 0; $i < @transtmp; $i++) {
      if ($transtmp[$i] =~ />/) {
	$count++;
	$trans->[$count]->{seq} = ''; # initialize
	if ($transtmp[$i] =~ /.*\[([0-9]{1,}):([0-9]{1,})\].*/) {
	  $trans->[$count]->{start} = $1;
	  $trans->[$count]->{end} = $2;
	  }
      }
      else {
	$trans->[$count]->{seq} .= $transtmp[$i];
      }
    }

    ## step 3: connect the fragments
    if (@$trans == 1) {
      $translated_seq = $trans->[0]->{seq};
    }
    else {
      for (my $i = 0; $i < @$trans; $i++) {
	$translated_seq .= $trans->[$i]->{seq};
	if ($i < (@$trans - 1)) {
	  my $missing = $trans->[$i+1]->{start} - $trans->[$i]->{end} -1;
	  if ($self->{get_indel} > 0) {
	  	$translated_seq .= 'X';
	  }
	}
      }
    }
    return($translated_seq);
  }

##################
sub codons {
    my $self = shift;
    my $finish = 0;
    my $codon_seq = '';
    my @transtmp;

    ## step 1: extract the relevant info from the genewise output
    for (my $i = 0; $i < $self->{gw_count}; $i++) {
      if ($self->{gw}->[$i] =~ />.*sp$/) {# the codons set starts
	while ($self->{gw}->[$i] !~ '//') {
	  push @transtmp, $self->{gw}->[$i];
	  $i++;
	}
	last; # end the for loop since nothing left to be done
      }
    }
    
    ## step two: get the sequences
    my $count = -1;
    my $trans;
    for (my $i = 0; $i < @transtmp; $i++) {
      if ($transtmp[$i] =~ />/) {
	$count++;
	$trans->[$count]->{seq} = ''; # initialize
	if ($transtmp[$i] =~ /.*\[([0-9]{1,}):([0-9]{1,})\].*/) {
	  $trans->[$count]->{start} = $1;
	  $trans->[$count]->{end} = $2;
	  }
      }
      else {
	$transtmp[$i] =~ tr/a-z/A-Z/;
	$trans->[$count]->{seq} .= $transtmp[$i];
      }
    }

    ## step 3: connect the fragments
    if (@$trans == 1) {
      $codon_seq = $trans->[0]->{seq};
    }
    else {
      for (my $i = 0; $i < @$trans; $i++) {
	$codon_seq .= $trans->[$i]->{seq};
	if ($i < (@$trans - 1)) {
	  my $indel = '';
	  my $missing = $trans->[$i+1]->{start} - $trans->[$i]->{end} -1;
	  ## now decide whether the nts that did not got translated are masked by
	  ## 'N' or whether they will be represented as lower case letters
	  if ($self->{get_indel} == 2) {
	    $indel = substr($self->{nt_seq}, $trans->[$i]->{end}, $missing);
	    $indel =~ tr/A-Z/a-z/;
	  }
	  elsif ($self->{get_indel} == 1) {
	    $indel = 'N' x $missing;
	  }
	  else {
	    $indel = '';
	  }
	  ## now append gap characters until the frame is recovered. Note that the gap
	  ## characters are added to the end of the indel-part. Thus, the codons are
	  ## not considered.
	  while (length($indel)%3 != 0) {
	    $indel .= '-';
	  }

	  $codon_seq .= $indel;
	}
      }
    }
    return ($codon_seq);
  }
###########################
sub protein_borders {
  my $self = shift;
  my $gw = $self->{gw};
  for (my $i = 0; $i < @$gw; $i++) {
    if ($gw->[$i] =~ /Bits.*introns$/) {
      my ($start, $end) = $gw->[$i+1] =~ /.*$self->{protname}\s{1,}([0-9]{1,})\s{1,}([0-9]{1,}).*/;
      return($start, $end);
    }
    else {
      die "no protein-start and end could not be determnined. Check genewise command\n";
    }
  }
}
##########################
sub cdna_borders {
  my $self = shift;
  my $gw = $self->{gw};
  for (my $i = 0; $i < @$gw; $i++) {
    if ($gw->[$i] =~ /Bits.*introns$/) {
      my ($start, $end) = $gw->[$i+1] =~ /.*$self->{dnaname}\s{1,}([0-9]{1,})\s{1,}([0-9]{1,}).*/;
      return($start, $end);
    }
    else {
      die "no cdna-start and end could not be determnined. Check genewise command\n";
    }
  }
}
##########################
sub _GetIndels {
  my $gw = shift;
  my $indel;
  for (my $i = 0; $i < @$gw; $i++) {
    if ($gw->[$i] =~ /Bits/) {
      $indel = $gw->[$i+1] =~ /.*([0-9]{1,})/;
      return($indel);
    }
  }
}
