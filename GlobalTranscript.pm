#Package- create Global Transcript intervals and Uniques Regions per genes
package PeakRescue::GlobalTranscript; 
use PeakRescue;
our $VERSION = PeakRescue->VERSION;

BEGIN {
  $SIG{__WARN__} = sub {warn $_[0] unless(( $_[0] =~ m/^Subroutine Tabix.* redefined/) || ($_[0] =~ m/^Use of uninitialized value \$buf/))};
};

#use lib qq(/software/CGP/canPipe/test/lib/perl5/);
use strict;
use Tabix;
use File::Path qw(mkpath);
use File::Basename;
use File::Spec;
use FindBin qw($Bin);
use List::Util qw(max min);
use Capture::Tiny qw(:all);
use warnings FATAL => 'all';
use Data::Dumper;
use Log::Log4perl;
Log::Log4perl->init("$Bin/../config/log4perl.test.conf");

my @std_chr=(1..22,'X','Y','MT');

my $log = Log::Log4perl->get_logger(__PACKAGE__);

use PeakRescue::GlobalTranscript;


=head2 new
Initialise the object and populate hash with options, create tabix object and do GTF processing
Inputs
=over 2
=item class
=item options -user input paramaters
=back
=cut

sub new {
	my ($class,$options)=@_;
	$log->trace('new');
	my $self={};
	bless $self, $class;
	$self->_init($options);
	$self->_create_sorted_gtf_tabix_object();
	$self->_check_tabix_overlap();
	return $self;
}

=head2 _init
populate object with  necessary file names
Inputs
=over 2
=item options -user input paramaters
=back
=cut

sub _init {
	my ($self,$options)=@_;
	my ($file_name,$dir_name,$suffix) = fileparse($options->{'gtf'},qr/\.[^.]*/);
	$options->{'s'}=$suffix;
	$options->{'d'}=$dir_name;
  if ($options->{'o'} && ! (-d $options->{'o'})) {
  	$log->debug("Creating dir:".$options->{'o'});
		mkpath($options->{'o'});
  }
	elsif(!$options->{'o'}) {
   $options->{'o'}=$dir_name;
  }
  elsif (!(-d $options->{'o'})) {
  	$log->logcroak("Unable to create directory");
  }
	$options->{'f'}=$options->{'o'}.'/'.$file_name;	
	$self->{'options'} = $options;
	
	return;
}

=head2 _create_sorted_gtf_tabix_object
Sort gtf file and create bgzip and tabix index for gtf file, using indexed gtf file create TABIX object
Inputs
=over 2
=item options -user input paramaters
=back
=cut

sub _create_sorted_gtf_tabix_object {
my ($self)=@_;
my $gtf_file=$self->options->{'gtf'};
my $gtf_file_no_ext=$self->options->{'f'};
	$log->debug("Using GTF file: ".$gtf_file);
	#check for bgzip and index file...
	if (-e "$gtf_file_no_ext.sorted.gz" && -e "$gtf_file_no_ext.sorted.gz.tbi" ) 
	{
		$log->debug("File $gtf_file_no_ext.sorted.gz and Tabix index exists: skipping this step ....");
	}
	#Create tabix indexed gtf...
	else{
		$log->debug("Creting sorted gtf and tabix index file");
		$log->debug("CR: create file ==>$gtf_file_no_ext<==.sorted.gz and tabix it");
		my ($out,$stderr,$exit) = capture{system("zless $gtf_file | sort -k1,1 -k4,4n | bgzip  >$gtf_file_no_ext.sorted.gz && tabix -p gff $gtf_file_no_ext.sorted.gz")};
		if($exit) {
			$log->logcroak('Unable to create sorted gtf file');
		}
	}
	#create Tabix object...
	$log->debug("Tabix indexed gtf file exits, Creating Tabix object");
	my $tabix_obj = new Tabix(-data => "$gtf_file_no_ext.sorted.gz");
	$log->debug("Tabix object created successfully");
	$self->{'tabix_obj'}=$tabix_obj;
	return;
}

=head2 
checks overlap with bed file
Inputs
=over 2
=item std_chr - reserved to process only user provided chr names [ all non standard chromosomes will be ignored ]
			currently using value from default array declared in this module
=back
=cut
sub _check_tabix_overlap {
	my ($self) = @_;
	my $tabix = $self->tabix_obj;
	my $gtf_record;
	# create outpt bed files names...
	my $unique_regions_bed=$self->options->{'f'}.'.unique_regions.bed';
	my $global_transcript_bed=$self->options->{'f'}.'.global_transcript.bed';
	my $geneboundaries_bed=$self->options->{'f'}.'.geneboundaries.bed';
	my $global_transcript_gene_length_bed=$self->options->{'f'}.'.global_transcript_gene_length.bed';
	my $unique_segment_gene_length_bed=$self->options->{'f'}.'.unique_segment_gene_length.bed';
	# create file handlers to use...
	my ($fh_array)=$self->_create_fh([$unique_regions_bed,$global_transcript_bed,$geneboundaries_bed,$global_transcript_gene_length_bed,$unique_segment_gene_length_bed],1);
	# deal with user defined chromosomes if any...
	my @chr_to_analyse;
	
	if (!$self->options->{'a'}) {
			$log->debug("Warning option -chr or  -a : Chromosomes to analyse not defined");
			$log->debug("using standard list of chromosomes defined for Human [ 1..22,X,Y,MT]");
			$log->debug("To process all chromosmes in a GTF file please provide option -a ");
		} 
	elsif($self->options->{'chr'}) {
		open(my $chr_user, '<', $self->options->{'chr'});
		@chr_to_analyse=<$chr_user>; 
	}
	elsif(defined $self->options->{'a'}) {
			$log->debug("Processing all chromosomes");
	}
	# get chromosome names
	my @chrnames=$tabix->getnames;
	###
	# Querying is ALWAYS half open regardless of underlying file type [ i.e zero start and 1 end -- same as bed file  ]
	###
	foreach my $chr (@chrnames) {
		# use user provided chromosome names
		if($self->options->{'chr'}) {
			next if !{map { $_ => 1 } split('\n', "@chr_to_analyse")}->{"$chr"};
		}
		elsif(!$self->options->{'a'}) {
			next if !{map { $_ => 1 } @std_chr}->{$chr};
		} 
		
		my $res = $tabix->query($chr);
		while(my $record = $tabix->read($res)){
		# the returned data is the ORIGINAL string from the file
			my ($chr,$exon,$start,$stop,$line) = (split '\t', $record) [0,2,3,4,8];
			next if $exon ne 'exon';
			my ($gene)= (split '\"', $line) [1];
			if(defined $gtf_record->{$gene}){
				# substract -1 from start to create bed interval as per BED format specification
				$gtf_record->{$gene}.="$chr\t".($start - 1)."\t$stop\t$gene\n";
			}
			else{
				$gtf_record->{$gene}="$chr\t".($start - 1)."\t$stop\t$gene\n";
				}
		}
		$log->debug(">>>>>>>Calculating global transcript intervals for chromosome:$chr ");
		$self->_process_gtf_lines_per_gene($gtf_record,$chr,@$fh_array[0],@$fh_array[1],@$fh_array[2],@$fh_array[3],@$fh_array[4]);
		$gtf_record=();
		$log->debug(">>>>>>>Done.......chr: $chr");
	}
	$self->_close_fh($fh_array);
	$log->debug(">>>>>>>GTF preprosessing completed successfully for:".$self->options->{'gtf'});
	
}

=head2 _process_gtf_lines_per_gene
get GTF lines per genes to merge and get GLOBAL transcript and unique intervals 
Inputs
=over 2
=item gtf_record -has containing GTF lines per chromosome
=item chr -chromosome number
=item fh_* file handlers to store respective data
=back
=cut

sub _process_gtf_lines_per_gene {
	my ($self,$gtf_record,$chr,$fh_unique_regions,$fh_global_transcript,$fh_geneboundaries,$fh_global_transcript_gene_length,$fh_unique_segment_gene_length)=@_;
	my $count=0;
	my $total=keys %$gtf_record;
	my $tmp_gene_file=$self->options->{'o'}.'/'.'tmp_gene_bed.txt';
	foreach my $gene (sort keys %$gtf_record) {
		$count++;
		if($count % 100 == 0) {
			$log->debug("Completed :$count genes of total:$total genes on chromosme:$chr");
		}
		my($fh)=$self->_create_fh([$tmp_gene_file],1);
		my $fh_str=@$fh[0];
		print $fh_str $gtf_record->{$gene};
		my ($unique_regions,$gt_intervals,$gt_start,$gt_end,$gt_length,$unique_seg_gene_length) = $self->_merge_intervals($gene,$tmp_gene_file);
		print $fh_unique_regions map { $_."\t$gene\n" } (split '\n', $unique_regions); 
		print $fh_global_transcript map { $_."\t$gene\n" } (split '\n', $gt_intervals);
		print $fh_geneboundaries "$chr\t$gt_start\t$gt_end\t$gene\n";
		print $fh_global_transcript_gene_length "$gene\t$gt_length\n";
		print $fh_unique_segment_gene_length "$gene\t$unique_seg_gene_length\n";
		close($fh_str);
	}
	return
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
	my ($self,$file_names,$overwrite)=@_;;
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


=head2 _merge_intervals
merge bed intervals per gene, get global transcript intervals, gt length and gt boundaries
Inputs
=over 2
=item gene - gene name to process
=item tmp_file file containing all intervals for a gene
=back
=cut

sub _merge_intervals {
	my ($self,$gene,$tmp_file)=@_;
	my ($g_start, $g_end, $chr, $start, $stop);
	# merge intervals to get global transcript intervals...
	my ($out,$stderr,$exit) = capture {system("mergeBed -i $tmp_file")};
	if($exit){ 
		$log->debug("Unable to perform mergeBed opertaion for $gene $stderr Trying to sort intervals");
		($out,$stderr,$exit) = capture {system("sortBed -i $tmp_file | mergeBed -i - ")};
	}
	if($exit){ 
		$log->logcroak("Unable to perform mergeBed opertaion for $gene $stderr  after sorting intervals");
	}
	
	my $total_len=0;
	foreach my $line((split("\n", $out))) {
		($chr, $start, $stop) = (split '\t', $line); 
		# Start, Stop 
		# no need to add +1 as bed positions are open ended... 
		$total_len+=($stop-$start);
		push(@$g_start,$start);
		push(@$g_end,$stop);
	}
	#get GT boundaries...
	my $gt_start = min(@$g_start);
	my $gt_end = max(@$g_end);
	
	my($unique_regions,$unique_seg_gene_length) = $self->_get_unique_regions($chr,$gt_start,$gt_end,$gene,$total_len,$out);
	return ($unique_regions,$out,$gt_start,$gt_end,$total_len,$unique_seg_gene_length);
}

=head2 _get_unique_regions
Get unique regions for a gene [subtract any regions that overlaps with one or more other genes]
Inputs
=over 2
=item chr -chromosome number
=item start zero_start start coordinate
=item stop one_end -stop coordinate
=item gene - gene name to process
=item output - output from merge bed
=back
=cut

sub _get_unique_regions {
	my ($self,$chr,$start,$stop,$gene,$total_len,$output)=@_;
	my $overlaped_gene_flag=0;
	my $tmp_sorted_bed=$self->options->{'o'}.'/'.'tmp_gene_sortedBed.txt';
	my $tmp_overlap_file=$self->options->{'o'}.'/'.'tmp_genes_overlap.txt';
	my($fh_overlap)=$self->_create_fh([$tmp_overlap_file],1);
	my $fh_overlap_str=@$fh_overlap[0];
	
	my $tabix=$self->tabix_obj;
	#extract all the GTF lines in a interval
	my $res = $tabix->query($chr,$start,$stop);
		while(my $record = $tabix->read($res)){	
		  my ($chr,$exon,$start,$stop,$line) = (split '\t', $record) [0,2,3,4,8];
			next if $exon ne 'exon';
			my ($tmp_gene)= (split '\"', $line) [1];
			# if GTF line doesn't match with underlying gene then consider this a overlap
			if($gene ne $tmp_gene) {
				print $fh_overlap_str "$chr\t".($start - 1)."\t$stop\t$tmp_gene\n";
				$overlaped_gene_flag=1;
			}
		}
		close($fh_overlap_str);
		# if overlapped genes found do subtract bed 
		if($overlaped_gene_flag){
			my($fh_sorted)=$self->_create_fh([$tmp_sorted_bed],1);
			my $fh_soretd_str=@$fh_sorted[0];
			print $fh_soretd_str $output;
			my($out,$stderr,$exit) = capture {system("subtractBed -a $tmp_sorted_bed -b $tmp_overlap_file | mergeBed -i - ")};
			if($exit){ 
				$log->debug("Unable to perform subtractBed opertaion for $gene $stderr Trying to sort Overlapped interval file file");
				($out,$stderr,$exit) = capture {system("sortBed $tmp_overlap_file | subtractBed -a $tmp_sorted_bed -b - | mergeBed -i - ")};
			}
			if($exit){ 
				$log->logcroak("Unable to perform subtractBed opertaion for $gene $stderr after sorting intervals");
			}
			close($fh_soretd_str);
			my $u_total_len=0;
			foreach my $line((split("\n", $out))) {
				my ($u_chr, $u_start, $u_stop) = (split '\t', $line); 
				# Start, Stop 
				# no need to add +1 as bed positions are open ended...
				$u_total_len+=($u_stop-$u_start);
			}
			
			
			return ($out,$u_total_len);
		}
		# No overlap found return original gene intervals as unique regions
		else {
			return ($output,$total_len);
		}
}


=head2 _close_fh
Close file handlers 
Inputs
=over 2
=item file_names -array of file handlers
=back
=cut
sub _close_fh {
	my ($self,$file_names)=@_;;
	foreach my $fh(@$file_names) {
		close($fh);
	}
}

sub tabix_obj {
shift->{'tabix_obj'};
}


sub options {
shift->{'options'};
}





1;
