#Package- create Global Transcript intervals and Uniques Regions per genes
package PeakRescue::GlobalTranscript; 
use PeakRescue;
our $VERSION = PeakRescue->VERSION;

BEGIN {
  $SIG{__WARN__} = sub {warn $_[0] unless(( $_[0] =~ m/^Subroutine Tabix.* redefined/) || ($_[0] =~ m/^Use of uninitialized value \$buf/))};
};

use strict;

use Tabix;
use File::Path qw(mkpath remove_tree);
use File::Basename;
use File::Spec;
use FindBin qw($Bin);
use List::Util qw(max min);
use Capture::Tiny qw(:all);
use Data::Dumper;
use Log::Log4perl;

use PeakRescue::Base;

Log::Log4perl->init("$Bin/../config/log4perl.gt.conf");
my $log = Log::Log4perl->get_logger(__PACKAGE__);

my @std_chr=(1..22,'X','Y','MT');

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
  $self->{'cfg_path'}=PeakRescue::Base->get_paths;	
	$self->_init($options);
	$self->_create_sorted_gtf_tabix_object();
	$self->_check_tabix_overlap();
	return $self;
}

sub cfg_path {
shift->{'cfg_path'};
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
	# create tmp dir
	mkpath($options->{'o'}.'/'.'tmp_gt');
	$options->{'tmpdir_gtf'}=$options->{'o'}.'/'.'tmp_gt';
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
	my $cmd = "zless $gtf_file | sort -k1,1 -k4,4n | $Bin/bgzip  >$gtf_file_no_ext\_sorted.gz";
	#check for bgzip and index file...
	if (-e "$gtf_file_no_ext\_sorted.gz" && -e "$gtf_file_no_ext\_sorted.gz.tbi" ) 
	{
		$log->debug("File $gtf_file_no_ext.sorted.gz and Tabix index exists: skipping this step ....");
	}
	#Create tabix indexed gtf...
	else{
		$log->debug("Creting sorted gtf and tabix index file");
	  PeakRescue::Base->_run_cmd($cmd);
	  $cmd="$Bin/tabix -p gff $gtf_file_no_ext\_sorted.gz"; 
	  PeakRescue::Base->_run_cmd($cmd);
	}
	#create Tabix object...
	$log->debug("Tabix indexed gtf file exits, Creating Tabix object");
	my $tabix_obj = new Tabix(-data => "$gtf_file_no_ext\_sorted.gz");
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
	my $unique_regions_bed=$self->options->{'f'}.'_unique_regions.bed';
	my $global_transcript_bed=$self->options->{'f'}.'_global_transcript.bed';
	my $geneboundaries_bed=$self->options->{'f'}.'_geneboundaries.bed';
	my $global_transcript_gene_length_bed=$self->options->{'f'}.'_global_transcript_gene_length.tab';
	my $unique_segment_gene_length_bed=$self->options->{'f'}.'_unique_segment_gene_length.tab';
	my $non_overlapping_genes=$self->options->{'f'}.'_non_overlapping_geneboundaries.bed';

	$self->options->{'unique_regions'}=$unique_regions_bed;
	$self->options->{'global_transcript'}=$global_transcript_bed;
	$self->options->{'geneboundaries'}=$geneboundaries_bed;
	$self->options->{'global_transcript_gene_length'}=$global_transcript_gene_length_bed;
	$self->options->{'unique_segment_gene_length'}=$unique_segment_gene_length_bed;
	$self->options->{'non_overlapping_geneboundaries'}=$non_overlapping_genes;
	
	if ( -e $unique_regions_bed) {return;}
	# create file handlers to use...
	my ($fh_array)=PeakRescue::Base->_create_fh([$unique_regions_bed,$global_transcript_bed,$geneboundaries_bed,$global_transcript_gene_length_bed,$unique_segment_gene_length_bed,$non_overlapping_genes],1);
	# add header lines
	$self->_add_header($fh_array);
	if($self->options->{'u'}) {
		$log->debug("Calculating unique regions");
	}else {
		$log->debug("To calculate unique regions please sepcify -u ");
	}
	
	# deal with user defined chromosomes if any...
	my $chr_to_analyse=undef;
	if (!$self->options->{'a'} && !$self->options->{'chr'} ) {
			$log->debug("Warning option -chr or  -a : Chromosomes to analyse not defined");
			$log->debug("using standard list of chromosomes defined for Human [ 1..22,X,Y,MT]");
			$log->debug("To process all chromosmes in a GTF file please provide option -a ");
		} 
	elsif($self->options->{'chr'}) {
		 @$chr_to_analyse=split(',',$self->options->{'chr'});
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
			#next if !{map { $_ => 1 } split('\n', "@chr_to_analyse")}->{"$chr"};
			next if (grep (/^$chr$/, @$chr_to_analyse) ne 1 );
		}
		elsif(!$self->options->{'a'}) {
			next if (grep (/^$chr$/, @std_chr) ne 1);
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
		$self->_process_gtf_lines_per_gene($gtf_record,$chr,@$fh_array[0],@$fh_array[1],@$fh_array[2],@$fh_array[3],@$fh_array[4],@$fh_array[5]);
		$gtf_record=();
		$log->debug(">>>>>>>Done.......chr: $chr");
	}
	PeakRescue::Base::_close_fh($fh_array);
	$log->debug(">>>>>>>GTF preprosessing completed successfully for:".$self->options->{'gtf'});
	#cleanup tmp folder
	#PeakRescue::Base->cleanup_dir($self->options->{'tmpdir_gtf'});
}

=head2 _add_header
Add header lines to each file
=item gtf_record -has containing GTF lines per chromosome
=item chr -chromosome number
=item fh_array file handlers to store respective data
=back
=cut

sub _add_header {
	my ($self,$fh_array)=@_;
	my $fh_unique_regions = @$fh_array[0];
	my $fh_global_transcript = @$fh_array[1];
	my $fh_geneboundaries = @$fh_array[2];
	my $fh_global_transcript_gene_length = @$fh_array[3];
	my $fh_unique_segment_gene_length = @$fh_array[4];
	my $fh_non_overlapping_genes = @$fh_array[5];
	print $fh_unique_regions "#chr\tstart\tend\tgene\n";
  print $fh_global_transcript "#chr\tstart\tend\tgene\n";
  print $fh_geneboundaries "#chr\tstart\tend\tgene\n";
  print $fh_global_transcript_gene_length "#gene\tgt_length\n";
  print $fh_unique_segment_gene_length "#gene\tUnique_seg_length\n";
  print $fh_non_overlapping_genes "#chr\tstart\tend\tgene\n";
 
 1;
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
	my ($self,$gtf_record,$chr,$fh_unique_regions,$fh_global_transcript,$fh_geneboundaries,$fh_global_transcript_gene_length,$fh_unique_segment_gene_length,$fh_non_overlapping_genes)=@_;
	my $count=0;
	my $total=keys %$gtf_record;
	my $tmp_gene_file=$self->options->{'tmpdir_gtf'}.'/tmp_gene_bed.txt';
	foreach my $gene (sort keys %$gtf_record) {
		$count++;
		if($count % 100 == 0) {
			$log->debug("Completed :$count genes of total:$total genes on chromosme:$chr");
		}
		my($fh)=PeakRescue::Base->_create_fh([$tmp_gene_file],1);
		my $fh_str=@$fh[0];
		print $fh_str $gtf_record->{$gene};
		my ($unique_regions,$gt_intervals,$gt_start,$gt_end,$gt_length,$unique_seg_gene_length,$overlaped_gene_flag) = $self->_merge_intervals($gene,$tmp_gene_file);
		if($self->options->{'u'}) {
			print $fh_unique_regions map { $_."\t$gene\n" } (split "\n", $unique_regions);
			print $fh_unique_segment_gene_length "$gene\t$unique_seg_gene_length\n";
			# write non overlapping gene boundaries, useful for coverage calculation
			if(!$overlaped_gene_flag) {
				print $fh_non_overlapping_genes "$chr\t$gt_start\t$gt_end\t$gene\n";
			}
		}
		print $fh_global_transcript map { $_."\t$gene\n" } (split "\n", $gt_intervals);
		print $fh_geneboundaries "$chr\t$gt_start\t$gt_end\t$gene\n";
		print $fh_global_transcript_gene_length "$gene\t$gt_length\n";		
		close($fh_str);
	}
	return
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
	my $cmd="$Bin/mergeBed -i $tmp_file";
	my ($out,$stderr,$exit) = capture {system($cmd)};
	if($exit){ 
		$log->debug("Unable to perform mergeBed opertaion for $gene $stderr Trying to sort intervals");
		$cmd="$Bin/sortBed -i $tmp_file | $Bin/mergeBed -i - ";
	  ($out) = PeakRescue::Base->_run_cmd($cmd);
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
	my $unique_regions=undef;
	my $unique_seg_gene_length=undef;
	my $overlaped_gene_flag=0;
	if($self->options->{'u'}) {
		($unique_regions,$unique_seg_gene_length,$overlaped_gene_flag) = $self->_get_unique_regions($chr,$gt_start,$gt_end,$gene,$total_len,$out);
	}
	return ($unique_regions,$out,$gt_start,$gt_end,$total_len,$unique_seg_gene_length,$overlaped_gene_flag);
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
	my $tmp_sorted_bed=$self->options->{'tmpdir_gtf'}.'/tmp_gene_sortedBed.txt';
	my $tmp_overlap_file=$self->options->{'tmpdir_gtf'}.'/tmp_genes_overlap.txt';
	my($fh_overlap)=PeakRescue::Base->_create_fh([$tmp_overlap_file],1);
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
			my($fh_sorted)=PeakRescue::Base->_create_fh([$tmp_sorted_bed],1);
			my $fh_soretd_str=@$fh_sorted[0];
			print $fh_soretd_str $output;
			my $cmd="$Bin/subtractBed -a $tmp_sorted_bed -b $tmp_overlap_file | $Bin/mergeBed -i - ";
			my($out,$stderr,$exit) = capture {system($cmd)};
			if($exit){ 
				$log->debug("Unable to perform subtractBed opertaion for $gene $stderr Trying to sort Overlapped interval file file");
	      $cmd="$Bin/sortBed $tmp_overlap_file | $Bin/subtractBed -a $tmp_sorted_bed -b - | $Bin/mergeBed -i - ";
				($out) = PeakRescue::Base->_run_cmd($cmd);
			}
			close($fh_soretd_str);
			my $u_total_len=0;
			foreach my $line((split("\n", $out))) {
				my ($u_chr, $u_start, $u_stop) = (split '\t', $line); 
				# Start, Stop 
				# no need to add +1 as bed positions are open ended...
				$u_total_len+=($u_stop-$u_start);
			}
			return ($out,$u_total_len,$overlaped_gene_flag);
		}
		# No overlap found return original gene intervals as unique regions
		else {
			return ($output,$total_len,$overlaped_gene_flag);
		}
}




sub tabix_obj {
shift->{'tabix_obj'};
}


sub options {
shift->{'options'};
}





1;
