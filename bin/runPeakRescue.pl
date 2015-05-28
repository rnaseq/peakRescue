#!/usr/bin/perl
BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  unshift (@INC,dirname(abs_path($0)).'/../lib');
  $SIG{__WARN__} = sub {warn $_[0] unless(( $_[0] =~ m/^Subroutine Tabix.* redefined/) || ($_[0] =~ m/^Use of uninitialized value \$buf/))};
};

use strict;
use FindBin qw($Bin);
use English qw( -no_match_vars );
use Pod::Usage qw(pod2usage);
use Getopt::Long;
use Try::Tiny qw(try catch finally);

use PeakRescue::RunPeakRescue;
my $log = Log::Log4perl->get_logger(__PACKAGE__);

try {
  my ($options) = option_builder();
  my $peak_rescue=PeakRescue::RunPeakRescue->new($options);
     $peak_rescue->run_pipeline;
  
}
catch {
 $log->logcroak($_);
};


sub option_builder {
	my ($factory) = @_;
	my %opts;
	&GetOptions (
					'h|help'    => \$opts{'h'},
					'bam|sampleBam=s' => \$opts{'bam'},
					'gtf|gtfFile=s' => \$opts{'gtf'},
					'g|genomeFasta=s' => \$opts{'g'},
					'alg|algorithm=s' => \$opts{'alg'},
					'o|outdir=s'  => \$opts{'o'},
					'v|version'  => \$opts{'v'},
	);

  pod2usage(-message => PeakRescue::license, -verbose => 1) if(defined $opts{'h'});
	
	if(defined $opts{'v'}){
		my $version = PeakRescue->VERSION;
		print "$version\n";
		exit;
	}
	pod2usage(q{'-gtf' gtf file must be specified.}) unless(defined $opts{'gtf'}) ;
	pod2usage(q{'-bam' bam file must be specified.}) unless(defined $opts{'bam'}) ;
	pod2usage(q{'-g' genome file must be specified.}) unless(defined $opts{'g'}) ;
	pod2usage(q{'-alg' at least one algorithm must be specified.}) unless(defined $opts{'alg'});
	pod2usage(q{'-o' output location must be specified}) unless(defined $opts{'o'});

	return \%opts;
}

__END__

=head1 NAME

runPeakRescue.pl - run PeakRescue pipeline 

=head1 SYNOPSIS

runPeakRescue.pl  -bam -gtf -g -a -o [ -h -v ]

Required Options (bam and bed interval files must be defined):

  --sampleBam         (-bam) sample bam file 
  --gtfFile           (-gtf) genome gtf file 
  --genomeFasta       (-g) fasta reference genome file 
  --algorithm         (-alg) algorithm to be use for coverage calculation [ biodbsam, clipover, mpileup, gatk ]
  --outdir            (-o) outdir [ Path to output directory ]
  
Optional :
  --help             (-h)  This message
  --version          (-v) displays version number of this software

  Example:
      perl runPeakRescue.pl -bam test.bam -gtf test.gtf -g test.fa -o testdir
=cut

