#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Spec;
use Getopt::Long;
use Pod::Usage qw(pod2usage);
use Data::Dumper;

# Author: sb43@sanger.ac.uk

my ($options) = option_builder();

open(my $fh1, $options->{'f1'}) || warn "print $!";
open(my $fh2, $options->{'f2'}) || warn "print $!";

my ($f1,$d1,$s1) = fileparse($options->{'f1'},qr/\.[^.]*/);
my ($f2,$d2,$s2) = fileparse($options->{'f2'},qr/\.[^.]*/);

my $col_f1=$options->{'c1'};
my $col_f2=$options->{'c2'};
chomp $col_f1;
chomp $col_f2;

my $fh_common;
# create file handle for filed to write [ common columns and different columns
if ($options->{'o'}) { 
	open($fh_common, '>', $options->{'o'});
}else{	
	open($fh_common, '>', $f1.'_'.$f2.'common'.$s1);
}

# create hash of both the files....

my ($file1_hash)=get_file_hash($fh1,$col_f1,$options->{'h'},$f1);
my ($file2_hash)=get_file_hash($fh2,$col_f2,$options->{'h'},$f2);

#combine files
combine_hash($file1_hash,$file2_hash,$fh_common,$options->{'h'});

# sub to combine file data...
sub combine_hash {
	my($hash_f1,$hash_f2,$fh_combined,$header)=@_;
	my $count=0;
	# Write common rows file...
	if($header){
		print $fh_common $hash_f1->{'header'}."\t".$hash_f2->{'header'}."\n";
	}
	#print Dumper $hash_f1;
	foreach my $key(keys %$hash_f1) {
		next if $key eq 'header';
		if ($hash_f2->{$key}) {
		  print $fh_combined $hash_f1->{$key}."\t".$hash_f2->{$key}."\n";
		  $count++;
		}
	}	
	print "Found $count matching lines\n"; 
}

# create hash from file
sub get_file_hash {
	my ($fh,$columns,$header,$f)=@_;
	my $file_hash;
	my $count=0;	
	while (my $line=<$fh>) {	
		$count++;
		$line=~s/\n//g;
		my ($col_key) = get_col_key($line,$columns);
		if($header && $count == 1) {
			$file_hash->{'header'}="$line";
			print "Columns to compare from file: $f: $col_key\n"; 
			next;
		} 
		$file_hash->{$col_key}="$line";
	}

return $file_hash


}


# get common key from both files
sub get_col_key {
	my ($line,$columns)=@_;
	my $col_key;
                my @col=split(',',$columns);
                foreach my $col_num(@col) {
                        my($field)=(split "\t", $line)[$col_num - 1];
                        $col_key.='_'.$field;
                }	
return $col_key;	
}



1;


#########################################

sub option_builder {
	my ($factory) = @_;

	my %opts;

	&GetOptions (
					'h|help'    => \$opts{'h'},
					'f1|file1=s' => \$opts{'f1'},
					'f2|file2=s'  => \$opts{'f2'},
					'c1|col1=s'  => \$opts{'c1'},
					'c2|col2=s'  => \$opts{'c2'},
					'o|outFile=s'  => \$opts{'o'},
					'hd|header=s'  => \$opts{'hd'},
					'v|version'  => \$opts{'v'},
	);

	pod2usage(q{'-f1' file1 must be specified.}) unless(defined $opts{'f1'}) ;
	pod2usage(q{'-f2' file2 must be specified.}) unless(defined $opts{'f2'}) ;
	return \%opts;
}

__END__
=head1 NAME

mergeFiles.pl - creates file with matching columns from both files    

=head1 SYNOPSIS

mergeFiles.pl [-h] -f1 -f2 -c1 -c2 -o  [-h -v ]

  Required Options (project must be defined):
    --help          (-h)  This message and format of input file
     One or more of the following:
     
    --file1  (-f1) file1 
    --file2  (-f2) file2 
    --col1  (-c1) comma separated column numbers fom file1  
    --col2  (-c2) comma separated columns numbers fom file2 
    --outFile  (-o) outfile name
    --header  (-hd) [1|0] print header line in output 

=cut


