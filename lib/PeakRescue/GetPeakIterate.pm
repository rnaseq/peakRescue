package PeakRescue::GetPeakIterate; 
use PeakRescue;
our $VERSION = PeakRescue->VERSION;

use strict;
use Bio::DB::Sam;

use Tabix;
use FindBin qw($Bin);
use List::Util qw(max min);
use File::Path qw(mkpath remove_tree);
use Capture::Tiny qw(:all);
use File::Basename;
use File::Spec;
use Data::Dumper;
use Const::Fast qw(const);
use Log::Log4perl;
Log::Log4perl->init("$Bin/../config/log4perl.gt.conf");
my $log = Log::Log4perl->get_logger(__PACKAGE__);

use PeakRescue::Base;

const my $PROPER_PAIRED => 0x2;
const my $UNMAPPED => 0x4;
const my $REVERSE_STRAND => 0x10;
const my $MATE_REVERSE_STRAND => 0x20;
const my $NOT_PRIMARY_ALIGN => 0x100;
const my $MATE_UNMAPPED => 0x0008;
const my $READ_PAIRED => 0x0001;
const my $FIRST_IN_PAIR => 0x40;
const my $MAX_PILEUP_DEPTH => '1000000';
const my $SUPP_ALIGNMENT => 0x800;
const my $DUP_READ => 0x400;
const my $VENDER_FAIL => 0x200;


sub new {
	my ($class,$options)=@_;
	$log->trace('new');
	my $self={};
	bless $self, $class;
  $self->{'cfg_path'}=PeakRescue::Base->get_paths;	
	$self->_init($options);
	$self->_do_max_peak_calculation();
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
	my ($file_name,$dir_name,$suffix) = fileparse($options->{'bam'},qr/\.[^.]*/);
	$options->{'s'}=$suffix;
	$options->{'d'}=$dir_name;
	$options->{'f'}=$file_name;
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
	mkpath($options->{'o'}.'/'.'tmp_peak');
	$options->{'tmpdir_peak'}=$options->{'o'}.'/'.'tmp_peak';
	$options->{'tmpdir_peak2'}=$options->{'o'}.'/'.'tmp_concatenate_peakouts';
	$options->{'tmpdir_peak3'}=$options->{'o'}.'/'.'tmp_intermediary_peakoutputs';
	$self->{'options'} = $options;
	$log->debug("<<<< Using bam file:".$options->{'bam'});
	$log->debug("<<<< Using bed file:".$options->{'bed'});
	$log->debug("<<<< Using coverage calculation method: ".$options->{'alg'});
	$log->debug("<<<< Using temp dir:".$options->{'tmpdir_peak'});
	return;
}


=head2  _get_bam_object
create biodbsam object for a given bam file
Inputs
=item bam_file - bam file to create object
=over 2
=back
=cut

sub _get_bam_object {
	my ($self,$bam_file)=@_;
		my $bam = Bio::DB::Sam->new(-bam => $bam_file,
															-fasta => $self->options->{'g'},
															-expand_flags => 1,
															-split_splices => 1
															);
		$bam->max_pileup_cnt($MAX_PILEUP_DEPTH); # changes default coverage cap from 8000 to $new_value
		return $bam;
}

=head2  _do_max_peak_calculation
main sub which does peak calculation using other subs
Also stores progress for every peak
Inputs
=over 2
=back
=cut

sub _do_max_peak_calculation { 
	my ($self)=@_;
	my($bam_object)=$self->_get_bam_object($self->options->{'bam'});
	my($tabix)=$self->_get_tabix_object($self->options->{'bed'});
	# -- Get the total number of genes to be processed from the .bed file itself
	my $geneBoundariesBedFile = $self->options->{'bed'};
	my $cmd = "wc -l $geneBoundariesBedFile | awk -F ' ' '{print \$1}'";
	my $totalNbGenesWithFileName = PeakRescue::Base->_run_cmd($cmd);
	my ($totalNbGenes, $tmpVarInputFile) = split / /, $totalNbGenesWithFileName;
	$totalNbGenes = $totalNbGenes-1; ## header line exists
	$log->debug(">>>>>>>>> In peak calculation: check the total nb of genes in " . $geneBoundariesBedFile . ": " . $totalNbGenes . " >>>>>>>>>");
	# --
	my($tabix_gt)=$self->_get_tabix_object($self->options->{'gt'});
	my $peak_out=$self->options->{'o'}.'/'.$self->options->{'f'}.'_peak.txt';
	$self->options->{'sample_peak'}=$peak_out;
	if (-e $self->options->{'sample_peak'} ) { $log->debug("Outfile exists:".$self->options->{'sample_peak'}."Skipping _do_max_peak_calculation step"); return;}
	open(my $peak_fh, '>>', $peak_out);
	# check progress required as this is long process if beaks in between can be restarted from where it left using progress file
	my $progress_file=$self->options->{'o'}.'/'.'progress.txt';
	if( -e $progress_file) { $log->info("Progress file exists, analysis will be resumed, to restart analysis please remove progress file"); }
	# open file in res write and append mode
	open(my $progress_fh, '>>', $progress_file);
	open(my $read_progress, '<', $progress_file);
	# -- Extend the progress file to also report the gene names already processed for each chromosome: chr [1-* relationship] gene pairs.
	my @chr_analysed = [] ;
	my %hash_chr_genes_processed;
	while(<$read_progress>) { 
		my ( $tag_chr, $tag_gene ) = split;
		$tag_gene =~ s/^\s+|\s+$//g;
		push @{ $hash_chr_genes_processed{$tag_chr} }, $tag_gene;
		push @chr_analysed, $tag_chr;
	}
	close($read_progress);
	# -- Check if the peak calculation was completed in the previous runs: i.e. is the number of genes in the hash hash_chr_genes_processed equal to the total nb of genes in boundaris BED
	my $size= 0 ;
	foreach ( values %hash_chr_genes_processed ) { $size += @$_ }
	$log->debug(">>>>>>>>> In peak calculation: check the total nb of genes: " . $totalNbGenes . " versus " . $size . " >>>>>>>>>");
	my $tag_return_peakoutput_filename = 1;
	if ($size != $totalNbGenes) { 
		# end of progress check	
		$log->logcroak("Unable to create tabix object") if (!$tabix);
		# get chromosome names
		my @chrnames=$tabix->getnames;
		# -- add gene_counter
		my $gene_counter = 0;
		foreach my $chr (@chrnames) {
			if (grep (/^$chr$/, @chr_analysed)) { 
				$log->debug(">>>>>>>Calculating peak for chromosome (continue): $chr ");
			} else {
				$log->debug(">>>>>>>Calculating peak for chromosome: $chr ");
			}
			# -- Extract list of genes already processed for this chromosome 
			my @genes_analysed = ();
			foreach my $tag_chr (keys %hash_chr_genes_processed) {
				foreach (@{$hash_chr_genes_processed{$tag_chr}}) {
					if ($tag_chr==$chr) {
						#print "$k:$_\n";
						push @genes_analysed, $_ ; 
					}
				}
			}
			# --
			my ($peak_data, $counter);
			my $res = $tabix->query($chr);
			while(my $record = $tabix->read($res)){		
					chomp;
					my ($chr,$start,$end,$gene)=(split "\t", $record) [0,1,2,3];
					$gene =~ s/^\s+|\s+$//g;
					# -- skip the genes already processed  
					next if (grep (/^$gene$/, @genes_analysed));
					$counter++;
					if($counter % 100 == 0) {
						$log->debug("Calculated peak for chromosome: $chr ==>".$counter.' genes');
						# can be used if hash is running out of memory
						$self->_print_peak($peak_data,$peak_fh);
						$peak_data=();
					}
					my($max_peak)=$self->_get_gene_bam_object($bam_object,$chr,$start,$end,$gene,$tabix_gt);
					$peak_data->{$gene}=$max_peak;
					# -- Record progress for each gene processed
					print $progress_fh $chr."\t".$gene."\n";
					# -- Check if 1K genes processed then exit to empty all created objects - workaround "too many open files" bug..
					$gene_counter++;
					if($gene_counter % 1000 == 0) { 
						$self->_print_peak($peak_data,$peak_fh);
						close($peak_fh);
						if ($tag_return_peakoutput_filename == 1) {
							print($peak_out);
							$tag_return_peakoutput_filename = 0;
						}
						exit;
					}
			}
			# -- print peak data for each chromosome
			$self->_print_peak($peak_data,$peak_fh);
			close($peak_fh);
			if ($tag_return_peakoutput_filename == 1) {
				print($peak_out);
				$tag_return_peakoutput_filename = 0;
			}
			$log->debug(">>>>>>>>> Completed peak calculation for chromosome: $chr ===>".$counter.' gene >>>>>>>>>');
		}
	} ## end of size check on nb of genes processed
	else {
		# -- Create an output file to specify that this step is complete
		# -- => log_peaks_for_all_genes.success in main dir ## not in the output dir as this is not known within the master shell script: run_pr_with_restart.sh
		# my $successOutputPeakCalculation = $self->options->{'o'}."/log_peaks_for_all_genes.success";
		my $successOutputPeakCalculation = "log_peaks_for_all_genes.success";
		my $cmd_success = "touch $successOutputPeakCalculation";
		PeakRescue::Base->_run_cmd($cmd_success);
		$log->debug(">>>>>>>>> In peak calculation: created an output file " . $successOutputPeakCalculation . " to specify that this step is complete  >>>>>>>>>");
		# -- Copy tmp concatenated file to the output directory to proceed with the rest of the pipeline
		my $tmpConcatenatedPeaksFile = $self->options->{'o'}."/tmp_concatenate_peakouts/tmp_peakout.tsv";
		my $cmd = "cp $tmpConcatenatedPeaksFile $peak_out";
		PeakRescue::Base->_run_cmd($cmd);
		$log->debug(">>>>>>>>> In peak calculation: copied " . $tmpConcatenatedPeaksFile . " to " . $peak_out . "  >>>>>>>>>");
	}
	close($peak_fh);
	close($progress_fh);
 # cleanup will happen in run script	
 #PeakRescue::Base->cleanup_dir($self->options->{'tmpdir_peak'});
	
}

=head2  _get_tabix_object
create tabix object for a user provided bed file
Inputs
=over 2
=item bed
=back
=cut

sub _get_tabix_object {
	my ($self,$bed_file)=@_;
	my $tabix_obj;
	my ($file_name,$dir_name,$suffix) = fileparse($bed_file,qr/\.[^.]*/);
  my $tmp_bed = $self->options->{'tmpdir_peak'}."/$file_name\_tabix.bed";
	my $cmd = "$Bin/bedtools sort -i $bed_file | $Bin/bgzip >$tmp_bed.gz";
  PeakRescue::Base->_run_cmd($cmd);
	$cmd="$Bin/tabix -p bed $tmp_bed.gz";	
  PeakRescue::Base->_run_cmd($cmd);
	$tabix_obj = new Tabix(-data => "$tmp_bed.gz");
	$log->debug("Tabix object created successfully for $bed_file");
	return $tabix_obj;
}

=head2  _get_gene_bam_object
fetch reads
Inputs
=over 2
=item bam_object - Bio::DB bam object
=item region - chr:start-end format region to get reads
=item gene - gene name 
=back
=cut

sub _get_gene_bam_object {
my ($self,$bam_object,$chr,$start,$end,$gene,$tabix_gt)=@_;
	my $tmp_gene_file=$self->options->{'tmpdir_peak'}.'/tmp_gene.bam';
	
	# --start: test bug fixing: too many files open...
	my %open_tmp_bam_files = ();
     	if (!defined $open_tmp_bam_files{$gene}) {
    		$open_tmp_bam_files{$gene} = Bio::DB::Bam->open($tmp_gene_file,'w');
     	}
     	#$self->{bam} = $open_tmp_bam_files{$gene};
	# This test replaces the creation of a $bam object to a hash containing the opened filehandler
	# -- end: test bug fixing: too many files open...
	# -- create bio db bam object
	# my $bam = Bio::DB::Bam->open($tmp_gene_file,'w');
	
	#write header
	#$bam->header_write($bam_object->header);
	# $self->{bam}->header_write($bam_object->header);
     	$open_tmp_bam_files{$gene}->header_write($bam_object->header);
	# create temp file
	my $read_flag=undef;
	$bam_object->fetch("$chr:$start-$end", sub {
		my $a = shift;
		my $flags = $a->flag;
		return if $flags & $NOT_PRIMARY_ALIGN;
		return if $flags & $VENDER_FAIL;
		return if $flags & $UNMAPPED;
		return if $flags & $DUP_READ;
		#return if $flags & $SUPP_ALIGNMENT;
		# exact match with XF tag value
		if ($a->get_tag_values('XF') eq $gene) {
			#$bam->write1($a->{'align'});
			#$self->{bam}->write1($a->{'align'});
			$open_tmp_bam_files{$gene}->write1($a->{'align'});
			$read_flag=1;
		}elsif($a->aux=~ m/XF:Z:$gene/){
			#$bam->write1($a->{'align'});
			#$self->{bam}->write1($a->{'align'});
			$open_tmp_bam_files{$gene}->write1($a->{'align'});
			$read_flag=1;
		}
	});
	#undef $bam;
	#delete $bam;
	#undef $self->{bam};
	#delete $self->{bam};
	undef $open_tmp_bam_files{$gene};
	# read flag is undefined 
	if (!$read_flag) {
		#print "No reads for :  $gene\n";
		return "0";
	} 
	Bio::DB::Bam->index_build($tmp_gene_file);

	# name sorted bam required for clipOver
	
	#----------biodbsam------------------- Option1
	if($self->options->{'alg'} eq 'biodbsam') {
	  	# create tmp file
		my ($gene_bam)=$self->_get_bam_object($tmp_gene_file);	
		my ($coverage) = $gene_bam->features(-type => 'coverage', -seq_id => $chr, -start => $start, -end => $end);
		return max($coverage->coverage);
	}
	#------------ClipOver------------------Option2
	elsif($self->options->{'alg'} eq 'clipover') {
	  	# create tmp file
		my $tmp_gene_sorted=$self->options->{'tmpdir_peak'}.'/tmp_gene_name_sorted';
		# read name sorted bam file
		Bio::DB::Bam->sort_core(1,$tmp_gene_file,$tmp_gene_sorted);
		my($clipped_bam)=$self->_run_clipOver("$tmp_gene_sorted.bam");
		# -- Workaround: start: too many open files...
		# my ($gene_bam)=$self->_get_bam_object($clipped_bam);	
		my %open_gene_bam = ();
		if (!defined $open_gene_bam{$gene}) {
			$open_gene_bam{$gene} = $self->_get_bam_object($clipped_bam);
		}
		# -- Workaround: end: too many open files...
		# my ($coverage) = $gene_bam->features(-type => 'coverage', -seq_id => $chr, -start => $start, -end => $end);
		my ($coverage) = $open_gene_bam{$gene}->features(-type => 'coverage', -seq_id => $chr, -start => $start, -end => $end);
		#delete $gene_bam;
		#delete $clipped_bam;
		#delete $tmp_gene_sorted;
		undef $open_gene_bam{$gene};
		return max($coverage->coverage);
	}
	#-------------mpileup------------------Option3
	elsif ($self->options->{'alg'} eq 'mpileup') {	
  	#need samtools 1.1 or above which takes care of overlapping read pairs... 
		my $cmd= $self->cfg_path->{'SAMTOOLS_1'}."/samtools mpileup $tmp_gene_file -d $MAX_PILEUP_DEPTH -A -f ".$self->options->{'g'}. " -r $chr:$start-$end --no-BAQ ";
		my($out)=PeakRescue::Base->_run_cmd($cmd);
		my ($max)=$self->_parse_pileup($out);
		return $max;
	}
	#------------ gatk ----------Option4
	if($self->options->{'alg'} eq 'gatk') {
		my ($max)=$self->_run_gatk($tmp_gene_file,$chr,$start,$end,$tabix_gt);
		return $max;
	}
	
}

=head2  _run_clipOver
run bamutils clipover to clip overlapping read pairs to calculate coverage peak for a gene
Inputs
=over 2
=item gene_bam bam file containing gene specific reads
=back
=cut

sub _run_clipOver {
	my($self,$gene_bam)=@_;
	my $tmp_gene_clipped=$self->options->{'tmpdir_peak'}.'/tmp_gene_clipped';
	my $cmd="$Bin/samtools view -h $gene_bam | $Bin/bam clipOverlap --readName --in - --out $tmp_gene_clipped.bam";
  PeakRescue::Base->_run_cmd($cmd);	
	# create coordinate sorted bam ...
	Bio::DB::Bam->sort_core(0,$tmp_gene_clipped.'.bam',$tmp_gene_clipped.'_coord');
	Bio::DB::Bam->index_build($tmp_gene_clipped.'_coord.bam');
	return $tmp_gene_clipped.'_coord.bam';
}

=head2  _run_gatk
run gatk and get peak coverage
Inputs
=over 2
=item gene_bam bam file containing gene specific reads
=back
=cut

sub _run_gatk {
	my($self,$gene_bam,$chr,$start,$end,$tabix_gt)=@_;
	#my($gt_int_file)=$self->_get_intervals($tabix_gt,$chr,$start,$end);
	my $tmp_cov=$self->options->{'tmpdir_peak'}.'/tmp.coverage';
	my $cmd = "java -Xmx2G -jar ".$self->cfg_path->{'GATK'}."/GenomeAnalysisTK.jar ".
	"-T DepthOfCoverage ".
	" -R ".$self->options->{'g'}.
	" -I $gene_bam ".
	" -o $tmp_cov ".
	" --countType COUNT_FRAGMENTS".
	" -U ALL ".
	" -L $chr:$start-$end".
	" --downsampling_type NONE".
	" -omitIntervals".
	" -l OFF".
	" -omitLocusTable".
	" -omitSampleSummary".
	" --outputFormat table ";
	
 PeakRescue::Base->_run_cmd($cmd);
 my($max_cov)=$self->_parse_gatk($tmp_cov);
 
 return $max_cov;
 
}

=head2  _get_intervals
get global transcript intervals
Inputs
=over 2
=item tabix_gt - tabix object containing global transcript intrevals
=item chr,start,end - gene boundaries
=back
=cut

sub _get_intervals{
	my($self,$tabix,$chr,$start,$end)=@_;
	my $tmp_gt_file=$self->options->{'tmpdir_peak'}.'/tmp.list';
	open(my $gt_fh,'>',$tmp_gt_file);
	my $res = $tabix->query($chr,$start,$end);
	while(my $record = $tabix->read($res)){		
		print $gt_fh $record;
	}
	return $tmp_gt_file;
}

=head2  _parse_gatk
parse gatk coverage file and retuen max coverage peak
Inputs
=over 2
=item tmp_cov - gatk coverage output for given gene interval
=back
=cut
sub _parse_gatk {
	my($self,$tmp_cov)=@_;
#Locus   Total_Depth     Average_Depth_sample    Depth_for_flux
#22:45898118     0       0.00    0
#22:45898119     0       0.00    0
	open(my $cov_fh,'<',$tmp_cov)|| $log->logcroak("unable to open $!");
	my $max;
	while (<$cov_fh>) {
		my ($cov) = (split "\t", $_)[2];
		$max = $cov if $cov > $max;
	}
	return $max;
}

# parse pilpeup output
# currently not in use
sub _parse_pileup {
	my ($self,$pileup_out)=@_;	
	my $max;
	foreach my $line((split("\n", $pileup_out))) 
	{
		my ($cov)=(split "\t" , $line) [3];
		$max = $cov if $cov > $max;
	}
	$max;
}

=head2  _print_peak
print peak data to file
Inputs
=over 2
=item hash - has string peak data 
=item fh - file handler to write data
=back
=cut

sub _print_peak {
	my ($self,$hash,$fh)=@_;
	foreach my $key (sort keys %$hash) {
		if(defined($hash->{$key}))
		{	
			print $fh "$key\t".$hash->{$key}."\n";	
		}
	}
}


sub options {
shift->{'options'};
}



1;
