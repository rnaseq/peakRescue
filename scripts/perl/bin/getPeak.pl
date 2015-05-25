###!/software/perl-5.16.3/bin/perl
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

use PeakRescue::GetPeak;
my $log = Log::Log4perl->get_logger(__PACKAGE__);

try {
  my ($options) = option_builder();
  PeakRescue::GetPeak->new($options);
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
					'bed|intervalBed=s' => \$opts{'bed'},
					'gt|globalTranscript=s' => \$opts{'gt'},
					'alg|algorithm=s' => \$opts{'alg'},
					'g|genomeFasta=s' => \$opts{'g'},
					'o|outdir=s'  => \$opts{'o'},
					'v|version'  => \$opts{'v'},
	);

  pod2usage(-message => PeakRescue::license, -verbose => 1) if(defined $opts{'h'});
	
	if(defined $opts{'v'}){
		my $version = PeakRescue->VERSION;
		print "$version\n";
		exit;
	}
	pod2usage(q{'-bed' bed file must be specified.}) unless(defined $opts{'bed'}) ;
	pod2usage(q{'-bam' bam file must be specified.}) unless(defined $opts{'bam'}) ;
	pod2usage(q{'-gt' global transcript must be specified.}) unless(defined $opts{'gt'}) ;
	pod2usage(q{'-g' genome file must be specified.}) unless(defined $opts{'g'}) ;
	pod2usage(q{'-alg' at least one algorithm must be specified.}) unless(defined $opts{'alg'}) ;

	return \%opts;
}

__END__

=head1 NAME

getPeak.pl - process bam and gene boundaries bed file 

Get max coverage peak per gene --this file is used in later step of PeakRescue:

=head1 SYNOPSIS

getPeak.pl -bam -bed   [ -o -h -v ]

Required Options (bam and bed interval files must be defined):

  --intervalBed       (-bed) gene interval bed file 
  --globalTranscript  (-gt) global transcript interval bed file 
  --sampleBam         (-bam) sample bam file 
  --algorithm         (-alg) algorithm to be use for coverage calculation [ biodbsam, clipover, mpileup, gatk ]
  --genomeFasta       (-g) fasta reference genome file 
  
Optional :

  --outdir           (-o) outdir [ Path to output directory ]
  --help             (-h)  This message
  --version          (-v) displays version number of this software

  Example:
    - process bam file <test.bam> and bed file <test.bed>
      perl getPeak.pl -bam test.bam -bed test.bed -o testdir
=cut

