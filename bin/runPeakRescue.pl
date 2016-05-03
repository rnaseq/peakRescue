#!/usr/bin/perl
BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  $ENV{POSIXLY_CORRECT}=1;
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
					'so|sortedByName=s' => \$opts{'so'},
					'gtf|gtfFile=s' => \$opts{'gtf'},
					'g|genomeFasta=s' => \$opts{'g'},
					'alg|algorithm=s' => \$opts{'alg'},
					'st|stranded=s' => \$opts{'st'},
					'o|outdir=s'  => \$opts{'o'},
					'v|version'  => \$opts{'v'},
	);

  pod2usage(-message => PeakRescue::license, -verbose => 1) if(defined $opts{'h'});
	
	if(defined $opts{'v'}){
		my $version = PeakRescue->VERSION;
		print "$version\n";
		exit(0);
	}
	if(!defined $opts{'so'}) {
		$opts{'so'}='no';
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

runPeakRescue.pl  -bam -gtf -g -o -alg [-st -h -v ]

Required Options (bam and gtf and genome files must be defined):

  --sampleBam         (-bam) sample bam file 
  --sortedByName      (-so) is input bam already sorted by name [yes,no :default is no] speedup analysis using name sorted bam file
  --gtfFile           (-gtf) genome gtf file 
  --genomeFasta       (-g) fasta reference genome file 
  --algorithm         (-alg) algorithm for coverage calculation [ biodbsam, clipover, mpileup, gatk ]
  --stranded          (-st) Stranded data [yes,no,reverse : default no ]
  --outdir            (-o) outdir [ Path to output directory ]
  
Optional:

  --help             (-h)  This message
  --version          (-v) displays version number of this software

  Example:
      perl runPeakRescue.pl -bam test.bam -gtf test.gtf -g test.fa -alg clipover -o testdir -so yes
=cut

