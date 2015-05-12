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
use Carp;
use Getopt::Long;
use Try::Tiny qw(try catch finally);
#use warnings FATAL => 'all';
#use autodie qw(:all);

use PeakRescue::GlobalTranscript;


try {
  my ($options) = option_builder();
  PeakRescue::GlobalTranscript->new($options);
}
catch {
 croak "\n\n".$_."\n\n" if($_);
};

sub option_builder {
	my ($factory) = @_;

	my %opts;

	&GetOptions (
					'h|help'    => \$opts{'h'},
					'gtf|gtf_file=s' => \$opts{'gtf'},
					'chr|include_chr=s' => \$opts{'chr'},
					'a|analyse_all_chr' => \$opts{'a'},
					'o|outdir=s'  => \$opts{'o'},
					'v|version'  => \$opts{'v'},
	);

  pod2usage(-message => PeakRescue::license, -verbose => 1) if(defined $opts{'h'});
	
	if(defined $opts{'v'}){
		my $version = PeakRescue->VERSION;
		print "$version\n";
		exit;
	}
	pod2usage(q{'-gtf' gtf must be specified.}) unless(defined $opts{'gtf'}) ;
	return \%opts;
}

__END__

=head1 NAME

processGTF.pl - process GTF file

Process GTF file and create following files to be used in later steps of PeakRescue:

1. unique_regions.bed -  unique and non overlapping intervals per gene
2. global_transcript.bed - A single transcript created after merging all known transcripts for a gene
3. geneboundaries.bed - Start and stop of global transcript
4. global_transcript_gene_length.bed - Gene length based on all exons in a global transcript 
5. unique_segment_gene_length.bed - Gene length based on unique segments of a gene


=head1 SYNOPSIS


processGTF.pl -gtf   [ -o -a -h -v ]

Required Options (GTF file must be defined):

  --gtf_file       (-gtf) gtf_file to use
  
Optional :

  --outdir           (-o) outdir [ Path to output directory ]
  --include_chr      (-chr) include chromosomes [File containing one chromosme per line - default to process standard chromosomes 1..22,X,Y,MT in Human genome]
  --analyse_all_chr  (-a) analyse all chromosomes
  --help             (-h)  This message
  --version          (-v) displays version number of this software

  Example:
    - process GTF file <test.gtf.gz> 
      perl processGTF.pl -gtf test.gtf.gz -o testdir
=cut

