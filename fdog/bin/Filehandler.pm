package Filehandler;
use strict;
# PROGRAM NAME: Filehandler.pm

# AUTHOR: INGO EBERSBERGER, ingo.ebersberger@univie.ac.at 

# PROGRAM DESCRIPTION: A module that retrieves a filename, path and
# input separator. It opens a file and hands back an object where via
# the command 'next' the next line of the file is fetched.

# DATE: 19.08.2003


# DATE LAST MODIFIED:


##################### start subroutine #######################
## blessing the variable:
## constructor that returns a file handle
sub TIEHANDLE {
    my $class = shift;
    my $name = shift;
    my $path = shift;
    $/ = shift;
    $path =~ s/\/$//;
    my $self;
    open ($self, "$path/$name") or die "could not open $path/$name\n";
    bless($self, $class);
    return ($self);
}

sub READLINE {
    my ($self) = shift;
    return <$self>;
}
sub CLOSE {
    my $self = shift;
    close ($self) or die "could not close filehandle\n";
}

sub PRINT {
    my $self = shift;
    print $self @_;
}
1;
