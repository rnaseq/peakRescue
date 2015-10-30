package PeakRescue::Base;
use PeakRescue;
our $VERSION = PeakRescue->VERSION;

use strict;
use Tabix;
use File::Path qw(mkpath remove_tree);
use File::Basename;
use File::Spec;
use FindBin qw($Bin);
use List::Util qw(max min);
use Capture::Tiny qw(:all);
use Config::IniFiles;
use Data::Dumper;
use Log::Log4perl;
use Const::Fast qw(const);
Log::Log4perl->init("$Bin/../config/log4perl.gt.conf");
my $log = Log::Log4perl->get_logger(__PACKAGE__);
const my $INI_FILE => "$Bin/../config/peakrescue.ini";

# recursively cleanups folder and underlying substructure
sub cleanup_dir {
    my ($self,$dir)=@_;
    $log->logcroak("Unable to find cleanup dir:$dir") if(! -d $dir);
    remove_tree($dir);
    $log->logcroak("Unable to remove cleanup dir:$dir") if( -d $dir);
    $log->debug("Dir: $dir cleaned successfully");
    return;
}


sub get_paths {
	my($self)=@_;
	tie my  %ini , 'Config::IniFiles', ( -file => $INI_FILE);
  \%{$ini{'PATHS'}}
}


sub cfg_path {
 shift->{'cfg_path'};
}


=head2 _create_fh
Create file handlers 
Inputs
=over 2
=item file_names -array of file names
=item overwrite - flag to overwrite existing file
=back
=cut

sub _create_fh {
	my ($self,$file_names,$overwrite)=@_;
	my $fh_array;
	foreach my $file(@$file_names) {
		if( -e $file && !$overwrite) {
			$log->debug("Output file $file already exists at:".$self->options->{'o'});
			next;
		}
		local *tmp_fh;
		open(tmp_fh,'>',$file) || $log->logcroak("Unable to create output file $!");
		push(@$fh_array,*tmp_fh);
	}
	
	return $fh_array;
}


=head2 _run_cmd
runs external command
Inputs
=over 2
=item cmd - command to run
=back
=cut

sub _run_cmd {
	my($self,$cmd)=@_;
	my ($out,$stderr,$exit)=capture{system($cmd)};
	if($exit) {
			$log->logcroak("Failed to run <<<<<<< \n $cmd  <<<<<< \n with status <<<<<< \n OUT:\n $out  :ERR:\n $stderr EXIT:\n $exit \n <<<<<<< \n");
	}
	#else {
	#	$log->debug("\ncommand <<<<<< \n $cmd \nrun successfully <<<<<<<<< ");
	#}
	return $out;
}

=head2 _close_fh
Close file handlers 
Inputs
=over 2
=item file_names -array of file handlers
=back
=cut

sub _close_fh {
	my ($file_names)=@_;;
	foreach my $fh(@$file_names) {
		close($fh);
	}
}













1;
