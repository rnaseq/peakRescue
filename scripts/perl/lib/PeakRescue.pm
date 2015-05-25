package PeakRescue;
use strict;
use Const::Fast qw(const);

use base 'Exporter';
our $VERSION = '5.2.0';
our @EXPORT = qw($VERSION);
const my $LICENSE =>
"#################
# PeakRescue::GlobalTranscript version %s, Copyright (C) 2015 
# This scripts comes with ABSOLUTELY NO WARRANTY
# See LICENSE for full details.
#################";

sub license {
  return sprintf $LICENSE, $VERSION;
}

1;

