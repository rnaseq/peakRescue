LICENCE
=======

PeakRescue is freely available under a GNU Public License.

PeakRescue
===========

PeakRescue is a robust fragment counting method for RNA-seq data.

This tool takes BAM files produced using any splice aware aligner e.g, TopHat, STAR etc., and  produces fragment count data per gene. 

---

##Dependencies

Python 2.5+ (< 3) - Please note peakRescue was tested with python 2.6.6 and 2.7.10.

Perl 5.14.2+

Tested on various Unix platforms [ Not tested on Mac OS X].

#### Note

Please read the Licensing agreement for respective tools before downloading them for commercial use.

#### Manual installation

###### Python package numpy

peakRescue was tested with numpy version 1.6.2:

https://pypi.python.org/pypi/numpy/1.6.2

###### Python package cython 

peakRescue was tested with cython version 0.22:

https://pypi.python.org/pypi/Cython/0.22

###### Please install the following tools before running setup.sh script and make sure that the path is set in config/peakrescue.ini file [optional if using 'clipover' as coverage calculation option]                                                                                                                             
                                                                                                                                                                                                                                                                                                                               
GATK - https://github.com/broadgsa/gatk-protected [ required only if user wants to use algorithm option: gatk for coverage calculation ]                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                                                               
https://github.com/broadgsa/gatk-protected/archive/2.8.tar.gz                                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                                                                               
Picard - https://github.com/broadinstitute/picard [ required only if user wants to use algorithm option: gatk for coverage calculation ]                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                                                               
https://github.com/broadinstitute/picard/releases/download/1.131/picard-tools-1.131.zip                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                               
Samtools - http://samtools.sourceforge.net [ required only if user wants to use algorithm option: mpielup for coverage calculation ]                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                                                               
https://github.com/samtools/samtools/archive/1.2.tar.gz     

#### Automated installation using setup.sh script (see section "Installing PeakRescue" below)

Bio::DB::Sam - http://search.cpan.org/~lds/Bio-SamTools/lib/Bio/DB/Sam.pm 

http://search.cpan.org/CPAN/authors/id/L/LD/LDS/Bio-SamTools-1.42.tar.gz

BamUtil - http://genome.sph.umich.edu/wiki/BamUtils 

http://genome.sph.umich.edu/w/images/7/70/BamUtilLibStatGen.1.0.13.tgz

Tabix - http://samtools.sourceforge.net/tabix.shtml 

http://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2

BedTools - http://bedtools.readthedocs.org/en/latest/

https://bedtools.googlecode.com/files/BEDTools.v2.17.0.tar.gz

#### Update .bashrc with the python version used to install numpy, cython

Edit your ~/.bashrc and add the following line - replacing 'full_path_to_your_python_installation_dir' to your own python path: alias python='full_path_to_your_python_installation_dir/python'

Save ~/.bashrc & quit.

Update environment variables with the following command: source ~/.bashrc. 

Check the use of the correct python path and version: which python && python -V

#### Download PeakRescue

wget https://github.com/rnaseq/peakRescue/archive/x.x.x.tar.gz

tar -xvzf x.x.x.tar.gz

##### for development only 

git clone https://github.com/rnaseq/peakRescue.git

# Installing PeakRescue

cd peakRescue-x.x.x/

/bin/bash ./setup.sh  \<path_to_install_dir\>

# Running PeakRescue

perl ./bin/runPeakRescue.pl -bam datasets/chr21.bam -gtf datasets/chr21.gtf.gz -g datasets/chr21.fa -alg clipover -o output

#### choosing algorithm [-alg] option                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                                                               
clipover - uses bamutils clipOverlap option to merge overlapping reads [ fast and gives simialr results to gatk on test data set ]                                                                                                                                                                                             
                                                                                                                                                                                                                                                                                                                               
gatk - uses gatk DepthOfCoverage [ slow but accurate coverage for overlapping read pairs from same fragment ]                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                                                                               
mpileup - uses samtools v1.1 and above [ faster than gatk - propelry calculates coverage for overlapping reads ]   

