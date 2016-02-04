package PeakRescue::RunPeakRescue; 
use PeakRescue;
our $VERSION = PeakRescue->VERSION;

use Bio::DB::Sam;
use strict;
use FindBin qw($Bin);
use List::Util qw(max min);
use File::Path qw(mkpath remove_tree);
use Capture::Tiny qw(:all);
use File::Basename;
use File::Spec;
use Data::Dumper;
use Log::Log4perl;
Log::Log4perl->init("$Bin/../config/log4perl.gt.conf");
my $log = Log::Log4perl->get_logger(__PACKAGE__);

use PeakRescue::Base;
use PeakRescue::GlobalTranscript;
use PeakRescue::GetPeak;


sub new {
	my ($class,$options)=@_;
	$log->trace('new');
	my $self={};
	bless $self, $class;
	$self->{'cfg_path'}=PeakRescue::Base->get_paths;
	$self->_init($options);
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
	mkpath($options->{'o'}.'/'.'tmp_runPeakRescue');
	$options->{'tmpdir_pipeline'}=$options->{'o'}.'/'.'tmp_runPeakRescue';
	if(!defined $options->{'st'}) {
	 $options->{'st'}='no';
	}	
	$self->{'options'} = $options;
	
	$log->debug("Using bam file:".$options->{'bam'});
	$log->debug("Using gtf file:".$options->{'gtf'});
	$log->debug("Using genome file:".$options->{'g'});
	$log->debug("Results will be written to :".$options->{'o'});
	$log->debug("Using temp dir:".$options->{'tmpdir_pipeline'});
	return;
}


sub cfg_path {
 shift->{'cfg_path'};
}

sub run_pipeline {
	my($self)=@_;
	$self->_run_htseq;
	$self->_run_htseq_disambiguate;
  $self->_process_sam;
  $self->_process_gtf;
  $self->_get_peak;
	#$log->debug("Exit for testing");
	#return 1;
  my ($flag)= $self->_runPeakrescue;
	if($flag) {
		$self->_process_output;
		$log->info("PeakRescue pipeline completed successfully");
		$log->info("Process log is written in peakrescue.log");
		$log->info("Results are stored in folder ".$self->options->{'o'});
		$log->info("Cleaning tmp data, ignore errors occured during cleanup as it might be due to file permissions");
		sleep(5);
		PeakRescue::Base->cleanup_dir($self->options->{'tmpdir_pipeline'});
		PeakRescue::Base->cleanup_dir($self->options->{'tmpdir_peak'});
		PeakRescue::Base->cleanup_dir($self->options->{'tmpdir_gtf'});
		$log->info("Cleanup completed successfully");
  }else {
	
		$log->debug("Analysis completed successfully");
	}
  return 1;
}


=head2 _run_htseq
run htseq modified version which outputs sam file with XF:Z flag and counts data containing multimapped reads ?? to be confirmed
Inputs
=over 2
=item gtf - GTF file
=item bam - bam file from splice aware aligner 
=back
=cut


sub _run_htseq {
	my($self)=@_;
	# store output ...
	$self->options->{'htseq_sam'}=$self->options->{'tmpdir_pipeline'}.'/'.$self->options->{'f'}.'_htseq.sam';
	$self->options->{'htseq_count'}=$self->options->{'tmpdir_pipeline'}.'/'.$self->options->{'f'}.'_htseq_count.out';
	if (-e $self->options->{'htseq_count'} ) { $log->debug("Outfile exists:".$self->options->{'htseq_count'}." Skipping <<< _run_htseq >>> step"); return;}
	# requires read name sorted sam file...
	my $cmd = "$Bin/samtools sort -on ".$self->options->{'bam'}.' '.$self->options->{'tmpdir_pipeline'}.'/tmpsort >'.$self->options->{'tmpdir_pipeline'}.'/tmpsort.bam';
	if(! -s $self->options->{'tmpdir_pipeline'}.'/tmpsort.bam'){
		PeakRescue::Base->_run_cmd($cmd);
	}
  else{
	$log->debug("File exists :".$self->options->{'tmpdir_pipeline'}.'/tmpsort.bam');
	}
	  $cmd = "$Bin/samtools view ".$self->options->{'tmpdir_pipeline'}."/tmpsort.bam | ".
		"python ".
		" $Bin/HTSeq-0.5.3p3_peakRescue/HTSeq/scripts/count_peakRescue_step1.py ".
			" --mode=union ".
			" --stranded=".$self->options->{'st'}.
			" --samout=".$self->options->{'htseq_sam'}.
			" --type=exon ".
			" --idattr=gene_id ". 
			"- ".
			$self->options->{'gtf'}.
			" >".$self->options->{'htseq_count'};
	# run command
	if(! -s $self->options->{'htseq_sam'}){
		PeakRescue::Base->_run_cmd($cmd);
	}
  else{
		$log->debug("File exists :".$self->options->{'htseq_sam'});
	}
	return 1;
}

=head2 _run_htseq_disambiguate
run htseq modified version to outputs sam file with additional XF:Z flag for disambiguated reads  
and counts data containing unique disambiguated reads , ambiguous reads and multimapped reads 
Inputs
=over 2
=item gtf - GTF file
=item htseq_sam - sam file from previous step
=back
=cut

sub _run_htseq_disambiguate {
	my($self)=@_;
	$self->options->{'disambiguated_sam'} =$self->options->{'tmpdir_pipeline'}.'/'.$self->options->{'f'}.'_disambiguated.sam';
	$self->options->{'disambiguated_count'}=$self->options->{'tmpdir_pipeline'}.'/'.$self->options->{'f'}.'_disambiguated_count.out';
	$self->options->{'multimapped_rngn'}=$self->options->{'tmpdir_pipeline'}.'/'.$self->options->{'f'}.'_multimapped_readname_gene_name.out';
	$self->options->{'amb_rngn'}=$self->options->{'tmpdir_pipeline'}.'/'.$self->options->{'f'}.'_ambiguous_readname_gene_name.out';
	if (-e $self->options->{'disambiguated_count'} ) { $log->debug("Outfile exists:".$self->options->{'disambiguated_count'}." Skipping <<< _run_htseq_disambiguate >>> step"); return;}
	
  my $cmd = "grep -P \"ambiguous|alignment_not_unique\" ".$self->options->{'htseq_sam'}.
		" | python ".
		" $Bin/HTSeq-0.5.3p3_peakRescue/HTSeq/scripts/count_peakRescue_step2.py ".
			" --mode=union ".
			" --stranded=".$self->options->{'st'}.
			" --samout=".$self->options->{'disambiguated_sam'}.
			" --type=exon ".
			"--idattr=gene_id ". 
			"- ".
			$self->options->{'gtf'}." ".
			$self->options->{'multimapped_rngn'}." ".
			$self->options->{'amb_rngn'}." ".
			" >".$self->options->{'disambiguated_count'};
	# run command
	PeakRescue::Base->_run_cmd($cmd);
	
}


=head2 _process_sam
create Karyotypic sorted bam file using picard required for GATK pcoverage calculation
Inputs
=over 2
=back
=cut

sub _process_sam {
	my ($self)=@_;
  my $tmp_combined_bam=$self->options->{'tmpdir_pipeline'}.'/'.$self->options->{'f'}.'_combined_sorted.bam';
	$self->options->{'kayrotypic'}=$self->options->{'tmpdir_pipeline'}.'/'.$self->options->{'f'}.'_kayrotypic.bam';
		
	if (-e $self->options->{'kayrotypic'} ) { $log->debug("Outfile exists:".$self->options->{'kayrotypic'}." Skipping <<< _process_sam >>> step"); return;}
  
	# add disambiguated reads containing additional XF:Z tags to original sam with updated 
 	my $cmd = "$Bin/samtools view -H ".
 	   $self->options->{'bam'}.
 	  " | cat -  ".$self->options->{'disambiguated_sam'}. " ".$self->options->{'htseq_sam'}." | "."$Bin/samtools view -bS -| ".
		"$Bin/samtools sort -o - ".$self->options->{'tmpdir_pipeline'}."/tmpsort_2 >$tmp_combined_bam";  # took 19 min to create 1.6GB file
   
    if (! -e $tmp_combined_bam ) {
    	PeakRescue::Base->_run_cmd($cmd);
    	Bio::DB::Bam->index_build($tmp_combined_bam);
    }
  	
    
    #$cmd="samtools index $tmp_combined_bam";
   	#PeakRescue::Base->_run_cmd($cmd);
 
 # creates bam to Karyotipic sorted bam
		$cmd = "java -Xmx2G -jar ".$self->cfg_path->{'PICARD'}."/picard.jar ReorderSam ".
						"I= $tmp_combined_bam ".
						"O= ".$self->options->{'kayrotypic'}.
						" REFERENCE= ".$self->options->{'g'}.
						" ALLOW_CONTIG_LENGTH_DISCORDANCE=true".
						" ALLOW_INCOMPLETE_DICT_CONCORDANCE=true";
		
		if ( ! -e $self->options->{'kayrotypic'} && $self->options->{'alg'} eq 'gatk') {								
			PeakRescue::Base->_run_cmd($cmd);
			Bio::DB::Bam->index_build($self->options->{'kayrotypic'});
		}
		else {
			$self->options->{'kayrotypic'}=$tmp_combined_bam;
		}
		#$cmd="samtools index ".$self->options->{'kayrotypic'};
		#PeakRescue::Base->_run_cmd($cmd);
		
}

=head2 _process_gtf
process gtf file using GlobalTranscript package
Inputs
=over 2
=back
=cut

sub _process_gtf {
	my ($self)=@_;
	$self->options->{'u'}=1;
	$self->options->{'a'}=1;
	my $gt=PeakRescue::GlobalTranscript->new($self->options);
	undef $gt;
}

=head2 _get_peak
get peak 
Inputs
=over 2
=back
=cut

sub _get_peak {
	my ($self)=@_;
	my $getPeak_options;
	$getPeak_options->{'bed'}=$self->options->{'geneboundaries'};
  $getPeak_options->{'bam'}=$self->options->{'kayrotypic'};
  $getPeak_options->{'g'}=$self->options->{'g'};
  $getPeak_options->{'gt'}=$self->options->{'global_transcript'};
  $getPeak_options->{'o'}=$self->options->{'o'};
  $getPeak_options->{'alg'}=$self->options->{'alg'};
 	my $peak=PeakRescue::GetPeak->new($getPeak_options);
  $self->options->{'peak_file'}=$peak->options->{'sample_peak'};
  undef $peak;
}

=head2 _runPeakrescue
run peakrescue  
Inputs
=over 2
=back
=cut

sub _runPeakrescue {
	my ($self)=@_;
	$self->options->{'final_output'}= $self->options->{'o'}.'/peakRescueFinalCount.out';	
	if (-e $self->options->{'final_output'} ) { $log->debug("Outfile exists:".$self->options->{'final_output'}." Skipping <<< _runPeakrescue >>> step"); return 0;}
	# Dev: if we would like to analyse all reads in single pass
  	#my $combinedReadNameGeneName = $self->options->{'tmpdir_pipeline'}.'/tmpCombinedRNGN.tab';
	#my $cmd_cat = "cat ".$self->options->{'amb_rngn'}." ".$self->options->{'multimapped_rngn'}." > $combinedReadNameGeneName";
	#PeakRescue::Base->_run_cmd($cmd_cat);
	# -- Assign multimapped reads' proportions to genes
	$log->info("Running probabilistic assignment of multimapped reads... ");
	$self->_runReadToGeneAssignment($self->options->{'multimapped_rngn'},'multimappers');
	$log->info("Probabilistic assignment of multimapped reads: completed.");
	# -- Assign ambiguous uniquely mapped reads' proportions to genes
	$log->info("Running probabilistic assignment of ambiguous uniquely mapped reads...");
	$self->_runReadToGeneAssignment($self->options->{'amb_rngn'},'ambiguous_unique');
	$log->info("Probabilistic assignment of ambiguous uniquely mapped reads: completed.");
	# Dev: if we would like to analyse all reads in single pass
	#$self->_runReadToGeneAssignment($combinedReadNameGeneName,'combined');
}

=head2 _runReadToGeneAssignment
runs python script _runReadToGeneAssignment
Inputs
=over 2
=back
=cut

sub  _runReadToGeneAssignment {
	my($self,$ReadNameGeneNameFile,$tag) = @_;
	$self->options->{$tag}=$self->options->{'o'}.'/results_peakrescue_readtype_'.$tag.'_all_genes.out';

	if( -e $self->options->{$tag}) { return;}
	my $cmd= "python $Bin/readToGeneAssignment.py -d ".
	$self->options->{'o'}.
	" -p ".$self->options->{'peak_file'}.
	" -m $ReadNameGeneNameFile".
	" -l ".$self->options->{'global_transcript_gene_length'}.
	" -t $tag";
	PeakRescue::Base->_run_cmd($cmd);	
}

=head2 _process_output
process output files to create final output file with FPKM 
Inputs
=over 2
=back
=cut

sub _process_output {
	my ($self)=@_;
	# merge files in pairs using ensid as common key 
	# instead of using paste merge files by ensid key
	
	my $cmd = "$Bin/mergeFiles.pl -f1 ".$self->options->{'htseq_count'}.
	" -f2 ".$self->options->{'disambiguated_count'}.
	" -c1 1 -c2 1 -o ".$self->options->{'tmpdir_pipeline'}.'/tmp_htseq_count_n_disambiguated_count.out';  
	
	PeakRescue::Base->_run_cmd($cmd);
	
	$cmd = "$Bin/mergeFiles.pl -f1 ".$self->options->{'ambiguous_unique'}.
	" -f2 ".$self->options->{'multimappers'}.
	" -c1 1 -c2 1 -o ".$self->options->{'tmpdir_pipeline'}.'/tmp_ambiguous_unique_n_multimappers.out';
	
	PeakRescue::Base->_run_cmd($cmd);
	
	$cmd = "$Bin/mergeFiles.pl -f1 ".$self->options->{'tmpdir_pipeline'}.'/tmp_htseq_count_n_disambiguated_count.out'.
	" -f2 ".$self->options->{'tmpdir_pipeline'}.'/tmp_ambiguous_unique_n_multimappers.out'.
	" -c1 1 -c2 1 -o ".$self->options->{'tmpdir_pipeline'}.'/tmp_count_n_contributions.out';
	
	PeakRescue::Base->_run_cmd($cmd);
	# add length
	$cmd = "$Bin/mergeFiles.pl -f1 ".$self->options->{'tmpdir_pipeline'}.'/tmp_count_n_contributions.out'.
	" -f2 ".$self->options->{'global_transcript_gene_length'}.
	" -c1 1 -c2 1 -o ".$self->options->{'tmpdir_pipeline'}.'/tmp_count_n_contributions_n_length.out';
	
	PeakRescue::Base->_run_cmd($cmd);
	
  #gene   uc   gene    uc_d    mm_tr   amb_tr  gene    amb_p   gene    mm_p    ensid   gene_len_gt 
  #0      1      2       3      4        5       6       7       8       9       10      11             
  
  my $total_read_count;
  my $all_data;
  open (my $fh_data, '<', $self->options->{'tmpdir_pipeline'}.'/tmp_count_n_contributions_n_length.out');
  while (<$fh_data>) {
    chomp;
    my($gene, $uc, $uc_d, $mm_tr, $amb_tr, $amb_p, $mm_p, $length)=(split "\t", $_)[0,1,3,4,5,7,9,11];
    my $final_count=($uc+$uc_d+$amb_p+$mm_p);
    $total_read_count+=$final_count;
    my $line;
    push(@$line,$uc,$uc_d,$amb_tr,$amb_p,$mm_tr,$mm_p,$final_count,$length);
    $all_data->{$gene}=$line;
  }

  eval{$total_read_count= int($total_read_count + $total_read_count/abs($total_read_count*2));};
 
  open(my $fh_final,'>', $self->options->{'final_output'});
	#fpkm loop 
	print $fh_final 
	"Gene\t".
	"uniqueCount\t".
	"uniqueDisambiguatedCount\t".
	"ambiguousUniqueToRescue\t".
	"ambiguousUniqueProportion\t".
	"multimapperToRescue\t".
	"multimapperProportion\t".
	"totalCount\t".
	"globalTranscriptLength\t".
	"fpkmTotalCount\n";
	
  foreach my $gene (keys %$all_data) {
    my $fpkm = ((10**9) *  @{$all_data->{$gene}}[6] ) / ($total_read_count * @{$all_data->{$gene}}[7] );
    print $fh_final "$gene\t".join("\t",@{$all_data->{$gene}})."\t".sprintf("%.3f",$fpkm)."\n";
  }
  
	close($fh_final,$fh_data);

}



sub options {
shift->{'options'};
}



1;
