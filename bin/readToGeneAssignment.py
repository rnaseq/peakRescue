#!/usr/bin/env python

import os
import sys
import optparse
import warnings
import traceback
import readToGeneAssignmentWithCython as r2g

"""
PeakRescue workflow: Assignment of ambiguously mapped reads to associated genes.

"""

#####################################################################################################################
## -- FUNCTION MAIN: 
## -- COLLECTS PASSED PARAMETERS & CALLS 'peakrescue_probabilistic_assignment' 
## -- USING THE readToGeneAssignmentWithCython module.
#####################################################################################################################
def main():
   
   optParser = optparse.OptionParser( 
      
      usage = "%prog [options] -d outdir -p peak_gene -m readname_genename_mappings -l gene_length_filename -r read_type ",
      
      description=
         "This script takes as input an output directory (outdir), " +
         "a list of peak contribution for each gene (peak_gene), " + 
         "a read name-to-gene names mappings (readname_genename_mappings), a gene (global " +
         "transcript) length table and a read type [either ambiguous or multimappers]. " +
         "It does the assignment of ambiguously mapped reads to associated genes based on " +
         "each given gene's highest per-base read coverage (peak) reported for each gene.")

   optParser.add_option( "-d", "--outputDirectory", type="string", dest="outdir",
      default = "", help = "Output directory " )

   optParser.add_option( "-p", "--peakInputFile", type="string", dest="peak_filename",
      default = "", help = "Input file containing the highest per-base read coverage (peak) " + 
      "for each gene (one gene per row) " )

   optParser.add_option( "-m", "--mappingsReadsToGenes", type="string", dest="mappings_reads2genes_filename",
      default = "", help = "Input file containing the mappings of ambiguous reads to their genes " )

   optParser.add_option( "-l", "--geneLengthFilename", type="string", dest="gene_length_filename",
      default = "", help = "Input file containing the length of each gene's global transcript. " )
     
   optParser.add_option( "-t", "--readType", type="choice", dest="readtype",
      choices = ( "ambiguous_unique", "multimappers" ), 
      default = "multimappers", help = "This option specifies the read type handled in the read name-to-gene names mappings " +
         "(choices: ambiguous_unique, multimappers; default: multimappers)" )
            
   optParser.add_option( "-v", "--verbose", action="store_true", dest="verbose",
      help = "suppress progress report and warnings" )

   if len( sys.argv ) == 1:
      print "pass here"
      optParser.print_help()
      sys.exit( 1 )

   (opts, args) = optParser.parse_args()
   
   if (not opts.outdir) or (not opts.peak_filename) or (not opts.mappings_reads2genes_filename) or (not opts.gene_length_filename) or (not opts.readtype):
      sys.stderr.write( sys.argv[0] + ": Error - Please see the list of required parameters to provide on the command line \n" )
      sys.stderr.write( "Call with '-h' to get usage information.\n" )
      sys.exit( 1 )
      
   try:
      r2g.peakrescue_probabilistic_assignment(opts.outdir, opts.peak_filename, opts.mappings_reads2genes_filename, opts.gene_length_filename, opts.readtype)
   except:
      sys.stderr.write( "Error: %s\n" % str( sys.exc_info()[1] ) )
      sys.stderr.write( "[Exception type: %s, raised in %s:%d]\n" % 
         ( sys.exc_info()[1].__class__.__name__, 
           os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
           traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
      sys.exit( 1 )


if __name__ == '__main__':
	main()
	
