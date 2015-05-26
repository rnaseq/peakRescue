LICENCE
=======

PeakRescue is freely available under a GNU Public License.

PeakRescue
===========

PeakRescue is the robust fragment counting method for RNA-seq data.

This tool takes BAM files produced using any splice aware aligner e.g, TopHat, STAR etc., and  produces fragment count data per gene. 

---

##Dependencies
Some of the code in PeakRescue package has dependencies on following utility perl modules.



Please install following tools before running this script and make sure that the path is set either in config/peakrescue.ini file or PATH environment variable 

###Note
Please read the Licensesing agreement for respective tools before downloading them for commercial use.

GATK - https://github.com/broadgsa/gatk-protected [ required only if you select algorithm option - gatk for coverage calculation ] 

https://github.com/broadgsa/gatk-protected/archive/2.8.tar.gz

Picard - https://github.com/broadinstitute/picard [ required only if you select algorithm option - gatk for coverage calculation ] 

https://github.com/broadinstitute/picard/archive/1.109.tar.gz

Bio::DB::Sam - http://search.cpan.org/~lds/Bio-SamTools/lib/Bio/DB/Sam.pm

http://search.cpan.org/CPAN/authors/id/L/LD/LDS/Bio-SamTools-1.41.tar.gz

BamUtil - http://genome.sph.umich.edu/wiki/BamUtils [ required only if you select algorithm option - clipover for coverage calculation

http://genome.sph.umich.edu/w/images/7/70/BamUtilLibStatGen.1.0.13.tgz

Samtools - http://samtools.sourceforge.net

https://github.com/samtools/samtools/archive/1.2.tar.gz

Tabix - https://github.com/samtools/tabix


Python2.7+

Perl 5.16+

BedTools - http://bedtools.readthedocs.org/en/latest/

https://bedtools.googlecode.com/files/BEDTools.v2.17.0.tar.gz

---

# Running PeakRescue

perl scripts/perl/bin/runPeakRescue.pl -bam datasets/chr21.bam -gtf datasets/chr21.gtf.gz -g datasets/chr21.fa -alg clipover -o results_dir
