#!/bin/bash 

# -------------------------------------------------------------------------
# -- peakRescue: Get input parameters and check their validity
# -------------------------------------------------------------------------
## Argument = -b $bamFile -g $gtfFile -f $fastaFile -a $algoParameter -o $outputDirectory

usage()
{
cat << EOF
usage: $0 options

This script runs peakRescue with the following input parameters:

OPTIONS:
   -h      Show this message
   -b      Input BAM file
   -g      Input GTF file
   -f      Input FASTA file
   -a      Algorithm used to calculate read coverage
   -o      Output directory

EOF
}

bamFile=
gtfFile=
fastaFile=
algoParameter=
outputDirectory=
while getopts “hb:g:f:a:o:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         b)
             bamFile=$OPTARG
             ;;
         g)
             gtfFile=$OPTARG
             ;;
         f)
             fastaFile=$OPTARG
             ;;
         a)
             algoParameter=$OPTARG
             ;;
         o)
             outputDirectory=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

# Check that any of the passed parameters are not equal to an empty string 
if [[ -z $bamFile ]] || [[ -z $gtfFile ]] || [[ -z $fastaFile ]] || [[ -z $algoParameter ]] || [[ -z $outputDirectory ]]
then
     usage
     exit 1
fi

# Check that the specified input files actually exist
if [ ! -f $bamFile ]; then
    echo "File $bamFile not found."
    usage
    exit 1
fi

if [ ! -f $gtfFile ]; then
    echo "File $gtfFile not found."
    usage
    exit 1
fi

if [ ! -f $fastaFile ]; then
    echo "File $fastaFile not found."
    usage
    exit 1
fi


# -------------------------------------------------------------------------
# -- Run peakRescue pipeline
# -------------------------------------------------------------------------
peakOutputFile=`/usr/bin/env perl -I ./lib/perl5/ ./bin/runPeakRescue.pl -bam $bamFile -gtf $gtfFile -g $fastaFile -alg $algoParameter -o $outputDirectory` 

tmpConcatenatePeakOutputs="$outputDirectory/tmp_concatenate_peakouts"
mkdir -p $tmpConcatenatePeakOutputs
cat $peakOutputFile >> $tmpConcatenatePeakOutputs/tmp_peakout.tsv

temp_dir="$outputDirectory/tmp_intermediary_peakoutputs"
mkdir -p $temp_dir
mv $peakOutputFile $temp_dir

