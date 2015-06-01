LICENCE
=======

PeakRescue is freely available under a GNU Public License.

PeakRescue
===========

PeakRescue is robust fragment counting method for RNA-seq data.

This tool takes BAM files produced using any splice aware aligner e.g, TopHat, STAR etc., and  produces fragment count data per gene. 

---

##Dependencies

Python2.5+

Perl 5.14.2+


#### Note

Please read the Licensesing agreement for respective tools before downloading them for commercial use.

#### Manual installation

Please install following tools before running setup.sh script and make sure that the path is set in config/peakrescue.ini file 

GATK - https://github.com/broadgsa/gatk-protected [ required only if user want to use algorithm option - gatk for coverage calculation ] 

https://github.com/broadgsa/gatk-protected/archive/2.8.tar.gz

Picard - https://github.com/broadinstitute/picard [ required only if user want to use algorithm option - gatk for coverage calculation ] 

https://github.com/broadinstitute/picard/releases/download/1.131/picard-tools-1.131.zip

Samtools - http://samtools.sourceforge.net [ required only if user want to use algorithm option - mpielup for coverage calculation ]

https://github.com/samtools/samtools/archive/1.2.tar.gz

#### Automated installation using setup.sh script

Bio::DB::Sam - http://search.cpan.org/~lds/Bio-SamTools/lib/Bio/DB/Sam.pm 

http://search.cpan.org/CPAN/authors/id/L/LD/LDS/Bio-SamTools-1.41.tar.gz

BamUtil - http://genome.sph.umich.edu/wiki/BamUtils 

http://genome.sph.umich.edu/w/images/7/70/BamUtilLibStatGen.1.0.13.tgz

Tabix - http://samtools.sourceforge.net/tabix.shtml 

http://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2

BedTools - http://bedtools.readthedocs.org/en/latest/

https://bedtools.googlecode.com/files/BEDTools.v2.17.0.tar.gz

---

# Installing PeakRescue

cd peakRescue
./setup.sh  /path_to_install_dir

# Running PeakRescue

bin/runPeakRescue.pl -bam datasets/chr21.bam -gtf datasets/chr21.gtf.gz -g datasets/chr21.fa -alg clipover -o output

#### choosing algorithm [-alg] option

gatk - uses gatk DepthOfCoverage [ slow but accurate coverage for overlapping read peairs from same fragment ]
clipover - uses bamutils clipOverlap option to merge overlapping reads [ fast and gives simialr results to gatk on test data set ]
mpileup - uses samtools v1.1 and above [ faster than gatk - propelry calculates coverage for overlapping reads ]
biodbsam - uses Bio::Db:Sam coverage method [ fastest but calculates overlapping read pairs twice] 

