#!/bin/bash 

# -------------------------------------------------------------------------
# -- MASTER SCRIPT: Run peakRescue pipeline 
# -- Implementation for servers with limited nb of open files 
# -- (i.e. ulimit -n <100K) 
# -- Iterate over sets of 1K genes in the peak calculation - exit to avoid
# -- the 'Too many open file' error and re-start program where left.
# -------------------------------------------------------------------------

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
# -- Iterate over sets of 1K genes in the peak calculation - exit to avoid
# -- the 'Too many open file' error and re-start program where left.
# -------------------------------------------------------------------------
while [ 1 ]; do
	if [ -e $PWD/log_peaks_for_all_genes.success ]; then                                                                                                                                                                                                                                                                       
		echo -n "Finished processing peak calculation for all genes. "
		echo -n "Resume pipeline to rescue/assign reads to genes..."
		/bin/bash run_peakrescue.sh -b $bamFile -g $gtfFile -f $fastaFile -a $algoParameter -o $outputDirectory
		exit
	else 
		/bin/bash run_peakrescue.sh -b $bamFile -g $gtfFile -f $fastaFile -a $algoParameter -o $outputDirectory
	fi
done

