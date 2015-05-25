import sys
import optparse
import itertools
import warnings
import traceback
import os, os.path

import HTSeq

########################################################################################
## -- Global variables
########################################################################################
debug = 0
#debug = 1
tag_report_instances_same_multiread_on_same_gene = 1
########################################################################################

class UnknownChrom( Exception ):
   pass


def is_read_in_gene_interval(readObject, features, dict_read_name_genes_names,ambiguous_tag,exons):
	"""
	Check whether a given read is in a gene interval.
	This function is called to check any ambiguously mapped read in a fragment (e.g. multimapper and ambiguous unique).
	
	Input:
		+ readObject: read object (either r[0] or r[1])
		+ features: gene intervals (see HTSeq documentation)
		+ dict_read_name_genes_names: mappings of ambiguous reads-to-genes 
		+ ambiguous_tag: tag specifying the type of ambiguous read (e.g. ambiguous unique or multimapped)
		  @todo: remove this 'ambiguous_tag' parameter not required as input.
		+ exons: contains object with exon intervals (see HTSeq documentation and count.py)
		
	Output:
		+ Binary (0|1): 1: read is in gene interval; 0: read is not present in gene interval.
		+ fs_genes: two variables now containing same information as fs_genes (for historical development reasons) 
		+ fs_exons: variable contains same information as fs_genes (for historical development reasons) 
		  @todo: fix method interface to remove 'fs_exons' as it contains same info as fs_genes redundant information.          
		+ dict_read_name_genes_names: hash tale containing the ambiguous reads-to-genes mappings
		+ ambiguous_tag: binary variable set to 1 if the read maps on more than one gene (0 if it maps on a single gene). 
		
	"""
	fs_exons = set()
	fs_genes = set()
	fs_gene_N= set()
	ambiguous_tag=0

	# -- Check if read is ambiguous i.e overlaps on two genes
	genes_ge2_tag=0
	fs_exons_temp= set()
	for iv3_temp, fs_exon_temp in features[ readObject.iv ].steps():
		fs_exons_temp = fs_exons_temp.union( fs_exon_temp )
		if len(fs_exons_temp) > 1:
			genes_ge2_tag=1
			break
	# -- Check if read interval overlaps with gene if not return 0
	if  len(fs_exons_temp) == 0:
		return (0, fs_genes, fs_exons, dict_read_name_genes_names, ambiguous_tag)
	
	threshold_split_region = 1
	tag_split_read = 0
	read_subregions = []
	regions_iv = []
	N_subregions = []
	N_size = 1	

	# -- If read overlaps on two genes then check if it is due to split mapping 
	if genes_ge2_tag:
		for cigar_region in readObject.cigar:
			if cigar_region.type == "N" and int(cigar_region.size) >= threshold_split_region:
				N_subregions = [ cigar_region.ref_iv ]
				N_subregion_iv = N_subregions[0]
				for iv_N, fs_exon_N in features[ N_subregion_iv ].steps():
					fs_gene_N.update(fs_exon_N)
				tag_split_read = 1
			# -- Check number of genes in M regions of CIGAR
			if cigar_region.type == "M":
				read_subregions.append(cigar_region.ref_iv)
		if tag_split_read:
				regions_iv = read_subregions
				fs_all_genes = set()
				for cigar_region in readObject.cigar:
					if cigar_region.type == "M":
						for read_m_interval, fs_exon in ( exons[cigar_region.ref_iv].steps() ):
					     		for exon in list(fs_exon):
								fs_all_genes = fs_all_genes.union( set([exon.name]) )
								if debug:
									print "gene---%s" %(fs_all_genes)	     		
								if ((cigar_region.ref_iv.start == exon.iv.start) or (cigar_region.ref_iv.end == exon.iv.end)):
									fs_genes = fs_genes.union( set([exon.name]) )
						if len(fs_genes) == 1:
							return (1, fs_genes, fs_exon_temp, dict_read_name_genes_names,ambiguous_tag)
						else:
							# Re-initialise fs_genes for each match position in CIGAR of split read
							fs_genes = set()
				ambiguous_tag = 1
				if debug:
					print "READ: AMBIGUOUS - flag=0"
				return (0, fs_all_genes ,fs_exon_temp, dict_read_name_genes_names,ambiguous_tag)
		else:
				## Read is not a split read and the read maps on 2 genes on the same location (i.e. ambiguous) 
				ambiguous_tag=1
				if debug:
					print "READ: AMBIGUOUS (condition no split) - flag=0"
				return (0,fs_exon_temp ,fs_exon_temp, dict_read_name_genes_names,ambiguous_tag)
	else:
		if debug:
			print "DEBUG:: single gene"
		regions_iv = [readObject.iv]
	# --
	fs_exons = set()
	fs_gene_unique=set()
	read_name = readObject.read.name
	results_subregions  = []
	unique_gene_found = 0
	for rObj_iv in regions_iv:
		for iv3, fs_exon in features[ rObj_iv ].steps():
			if len(fs_exon) == 1:
				fs_gene_unique = fs_gene_unique.union( fs_exon )		
				unique_gene_found=1
			fs_exons = fs_exons.union( fs_exon )
			if unique_gene_found:
				fs_genes = fs_gene_unique
			else:
				fs_genes=fs_exons
	if len(fs_gene_unique) == 1 and len(fs_gene_unique - fs_exons) == 0:
		# Read is in gene interval and maps on single gene
		# >>> a=set(["A","B"])  #>>> u=set(["A"]) #>>> u-a #set([]) #>>> a-u #set(['B'])	
		return (1, fs_gene_unique, fs_gene_unique, dict_read_name_genes_names,ambiguous_tag)

	elif len(fs_gene_unique) > 1:
		ambiguous_tag=1
		return (0, fs_genes, fs_exons, dict_read_name_genes_names,ambiguous_tag)

	if len(fs_genes) != 0:
		ambiguous_tag=1
	return (0, fs_genes, fs_exons, dict_read_name_genes_names,ambiguous_tag)


def _populate_read_name_gene_name(dict_read_name_genes_names, fs_genes, read_name, tag_report_instances_same_multiread_on_same_gene=1):
	"""
	Store the set of genes 'fs_genes' a multiread instance maps to in a hash table (dict_read_name_genes_names) with the multiread 
	name as key ('read_name' argument). 
	If the multiread (read_name) was already present in the hash table, then add the new set of genes (on which the current 
	instance of the multiread maps to) to the existing genes (already stored for previously encountered instances of the same multiread).
	
	If multiple (>=2) instances of the same multiread (i.e. same 'read_name') map on the same gene, two options are available to store
	these instances:
		+ if tag_report_instances_same_multiread_on_same_gene is set to 1 ## default option.
			=> Report as many occurrences of a given gene as there are instances of a given multiread mapping to it.
		+ if tag_report_instances_same_multiread_on_same_gene is set to 0:
			=> Report *only* once the given gene regardless of how many instances of a given multiread map to it.

	Input:
		+ dict_read_name_genes_names: hash table containing read name-to-gene names mappings.
		+ fs_genes: set of genes
		+ read_name: a given read name
		+ tag_report_instances_same_multiread_on_same_gene (default to 1): see description in text above.
	
	Output:
		+ dict_read_name_genes_names: hash table containing read name-to-gene names mappings.

	N.B. Same rational applies for the ambiguous uniquely mapped reads.
	
	"""
	if read_name not in dict_read_name_genes_names.keys(): 
		## First instance of a given multiread encountered
		if tag_report_instances_same_multiread_on_same_gene == 1:
			dict_read_name_genes_names[read_name] = list(fs_genes)
		else:
			dict_read_name_genes_names[read_name] = fs_genes
	else:
		## Additional instances of a given multiread encountered
		if tag_report_instances_same_multiread_on_same_gene == 1:
			dict_read_name_genes_names[read_name].extend( list(fs_genes) )
		else:
			genes_prev_instances_of_same_multiread = dict_read_name_genes_names[read_name]
			dict_read_name_genes_names[read_name] = genes_prev_instances_of_same_multiread.union( fs_genes )
	return dict_read_name_genes_names


def add_unique_counts_per_feature(dict_gene_unique_counts, fs):
	"""
	Store the count of non-ambiguous unique reads for each gene in a hash table.
	 
	"""
	feature_name1 = list(fs)[0] 
	if feature_name1 in dict_gene_unique_counts.keys():
		dict_gene_unique_counts[feature_name1] += 1 
	elif feature_name1 not in dict_gene_unique_counts.keys():
		dict_gene_unique_counts[feature_name1] = 1 
	return (dict_gene_unique_counts)


def add_non_unique_counts_per_feature(fs_genes, dict_nonunique):
	"""
	Store the count of multimapped reads for each gene in a hash table.
	 
	"""
	for gene_feature in list(fs_genes):
		dict_nonunique[ gene_feature ] += 1
	return (dict_nonunique)

def add_unique_counts_per_feature_ambiguous(fs_genes, dict_unique_ambiguous):
	for gene_feature in list(fs_genes):
		dict_unique_ambiguous[ gene_feature ] += 1
	return (dict_unique_ambiguous)

def initialise_counts_per_feature(dict_features, feature_name):
	"""
	Initialise read count (set to 0) per gene feature.
			
	"""
	dict_features[feature_name] = 0
	return dict_features

def initalize_read_name_and_interval(read_object0, read_object1):
	"""
	Store previous fragment information (intervals for each read (r0 r1))
	To avoid counting twice reads that map at same location with distinct CIGARs.

	"""
	temp_interval0 = None
	temp_interval1 = None
	if (read_object0 is not None and read_object1 is None):
		temp_interval0 = read_object0.iv
		temp_read_name = read_object0.read.name
	elif (read_object0 is None and read_object1 is not None):
		temp_interval1 = read_object1.iv
		temp_read_name = read_object1.read.name
	elif (read_object0 is not None and read_object1 is not None):
		temp_interval0 = read_object0.iv
		temp_interval1 = read_object1.iv
		temp_read_name = read_object0.read.name
	return (temp_read_name, str(temp_interval0), str(temp_interval1))

def invert_strand( iv ):
   iv2 = iv.copy()
   if iv2.strand == "+":
      iv2.strand = "-"
   elif iv2.strand == "-":
      iv2.strand = "+"
   else:
      raise ValueError, "Illegal strand"
   return iv2

def _add_non_unique_counts_per_feature_based_on_read_interval_and_readname(readObject, temp_interval_read, temp_read_name, fs_genes, dict_nonunique, dict_read_name_genes_names):
	"""
	Store multimapped read count for each gene.
	
	Function called in 'count_reads_in_features' and in the following condition (multimapper reads):
	               if (( r[0] is not None and r[0].optional_field( "NH" ) > 1 ) or \
                     ( r[1] is not None and r[1].optional_field( "NH" ) > 1 )):
	Input:
		+ readObject: HTSeq read object (See HTSeq documentation)
		+ temp_interval_read: read interval 
		+ temp_read_name: read name
		+ fs_genes: set of genes 
		+ dict_nonunique: hash table storing mutlimapped read count per gene (keys) 
		+ dict_read_name_genes_names: hash table storing the read name to gene names mappings
		
	Output:
		+ dict_nonunique: hash table storing mutlimapped read count per gene (keys) 		
		+ flag_aln_not_unique
		+ dict_read_name_genes_names: hash table storing the read name to gene names mappings
	
	N.B. Parameter: tag_report_instances_same_multiread_on_same_gene is a global variable.
	
	"""
	flag_aln_not_unique = 0
	if ( temp_interval_read != str(readObject.iv)): 
		(dict_nonunique)= add_non_unique_counts_per_feature(fs_genes, dict_nonunique)
		dict_read_name_genes_names = _populate_read_name_gene_name(dict_read_name_genes_names, fs_genes, readObject.read.name, tag_report_instances_same_multiread_on_same_gene)
		flag_aln_not_unique = 1
	else:
		if (temp_read_name != readObject.read.name):
			(dict_nonunique)= add_non_unique_counts_per_feature(fs_genes, dict_nonunique)
			dict_read_name_genes_names = _populate_read_name_gene_name(dict_read_name_genes_names, fs_genes, readObject.read.name, tag_report_instances_same_multiread_on_same_gene)
			flag_aln_not_unique = 1
	return (dict_nonunique, flag_aln_not_unique, dict_read_name_genes_names)


############################################################################################################################################################
## -- MAIN FUNCTION: COUNT READS PER GENE
############################################################################################################################################################
def count_reads_in_features( sam_filename, gff_filename, stranded, overlap_mode, feature_type, id_attribute, quiet, minaqual, samout, \
							 filename_read_names_gene_names,filename_read_names_gene_names_amb_unique):
   """
	Main function to count reads in features i.e. genes. 
	
	Input:
		+ sam_filename: Input alignment with all the ambiguously mapped reads
		+ gff_filename: GTF containing all genes for a given species
		+ stranded: specify whether data are stranded - see -s option
		+ overlap_mode: mode to handle reads overlapping more than one feature (e.g. union) - 
		  See -m option: choices = ( "union", "intersection-strict", "intersection-nonempty")
		+ feature_type: see -t option
		+ id_attribute: see -i option
		+ quiet: see -q option
		+ minaqual: see -a option 
		+ samout: SAM output file storing disambiguated reads (see -o option).
		+ filename_read_names_gene_names: filename for the output file containing the mappings readName to geneNames for multimapped reads
		+ filename_read_names_gene_names_amb_unique: filename for the output file containing the mappings readName to geneNames for ambiguously mapped reads
      
	Output:
		+ Writes readName to geneName outputs.
		+ Writes SAM output file for ddisambiguated uniquely mapped reads.
		+ Writes to stdout the genes and their read counts with read count for distinct read type: non-ambiguous unique, multimapped and ambiguous unique. 
		  This output redirected and stored to an output file in main peakRescue pipeline. 
		  This output is used in the later stage of the peakRescue pipeline to rescue the reads present in the readName to genNames mappings.
	
   """
   # Output filhandles for readName to geneNames mappings
   fh_read_names_gene_names = open(filename_read_names_gene_names, 'w')
   fh_read_names_gene_names_amb_unique = open(filename_read_names_gene_names_amb_unique, 'w')
   
   def write_to_samout( r, assignment ):
      if samoutfile is None:
         return
      if not pe_mode:
         r = (r,)
      for read in r:
         if read is not None:
            samoutfile.write( read.original_sam_line.rstrip() + 
               "\tXF:Z:" + assignment + "\n" )
   if quiet:
      warnings.filterwarnings( action="ignore", module="HTSeq" ) 
      
   if samout != "":
      samoutfile = open( samout, "w" )
   else:
      samoutfile = None
      
   features = HTSeq.GenomicArrayOfSets( "auto", stranded != "no" )     
   ## Hash table to store unique reads per exon (if modified GTF)
   counts = {}
   ## Hash table to store original non unique reads per gene (without 
   dict_nonunique = {}
   ## Hash table to store all unique reads as per original GTF
   dict_gene_unique_counts = {}
   ## hast table to store ambigouous read count for unique reads...
   dict_gene_unique_counts_ambiguous = {}
   ## Hash table to store all non-unique reads including shared reads 
   ## (either split reads or read pair matching on two distinct exons, same gene)
   dict_gene_nonunique_counts = {}
   ## Hash to store the non-unique read-names as key and genes names as values (fragments)
   dict_read_name_genes_names = {}
   ## Hash to store the non-unique read-names as key and genes names as values (fragments) including instances of a given multimapped read on same gene
   dict_read_name_genes_names_final = {} 
   dict_read_name_genes_names_ambiguous = {}
   ## @todo: tag_gff - parameter to be removed - only deal with gene level information 
   ## tag_gff: type to specify whether it contains gene or exons information 
   tag_gff = "gene_gff" 
   # Try to open samfile and fail early in case it is not there
   if sam_filename != "-":
      open( sam_filename ).close() 
      
   gff = HTSeq.GFF_Reader( gff_filename )   
   exons = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
   
   i = 0
   try:
      for f in gff:
         if f.type == feature_type:
	    exons[ f.iv ] += f # added to get exon interval data
            try:
               feature_id = f.attr[ id_attribute ]
            except KeyError:
               sys.exit( "Feature %s does not contain a '%s' attribute" % 
                  ( f.name, id_attribute ) )
            if stranded != "no" and f.iv.strand == ".":
               sys.exit( "Feature %s at %s does not have strand information but you are "
                  "running htseq-count in stranded mode. Use '--stranded=no'." % 
                  ( f.name, f.iv ) )
            features[ f.iv ] += feature_id
            counts[ f.attr[ id_attribute ] ] = 0
	    # -- Initialisation 
	    feature_name = f.attr[ id_attribute ]
	    # -- Added tag_gff for GFF type
	    if tag_gff == "gene_gff":
		# Original GTF (genes) 
		dict_nonunique = initialise_counts_per_feature(dict_nonunique, feature_name)
		dict_gene_unique_counts = initialise_counts_per_feature(dict_gene_unique_counts, feature_name)
		dict_gene_nonunique_counts = initialise_counts_per_feature(dict_gene_nonunique_counts, feature_name)
		dict_gene_unique_counts_ambiguous = initialise_counts_per_feature(dict_gene_unique_counts_ambiguous, feature_name)
         i += 1
         if i % 100000 == 0 and not quiet:
            sys.stderr.write( "%d GFF lines processed.\n" % i )
   except:
      sys.stderr.write( "Error occured in %s.\n" % gff.get_line_number_string() )
      raise
      
   if not quiet:
      sys.stderr.write( "%d GFF lines processed.\n" % i )
      
   if len( counts ) == 0 and not quiet:
      sys.stderr.write( "Warning: No features of type '%s' found.\n" % feature_type )
   
   try:
      if sam_filename != "-":
         read_seq = HTSeq.SAM_Reader( sam_filename )
         first_read = iter(read_seq).next()
      else:
         read_seq = iter( HTSeq.SAM_Reader( sys.stdin ) )
         first_read = read_seq.next()
         read_seq = itertools.chain( [ first_read ], read_seq )
      pe_mode = first_read.paired_end
      #pe_mode = 1 ## Added by us
   except:
      sys.stderr.write( "Error occured when reading first line of sam file.\n" )
      raise

   ###################################################################################################   
   try:
      if pe_mode:
         read_seq_pe_file = read_seq
         read_seq = HTSeq.pair_SAM_alignments( read_seq )
      empty = 0
      ambiguous = 0
      ambiguous_tag=0
      notaligned = 0
      lowqual = 0
      nonunique = 0
      nonunique_nonamb_to_be_rescued = 0
      temp_read_name="NA"
      previous_read_name="NA"
      temp_interval_r0="NA"
      temp_interval_r1="NA"
      counter_fragment = 0	
      flag_result = 0
      i = 0   
      pe_mode_for_SE = 0
      ## -- Added pe_mode on for SE files so that multireads reads will be accounted for
      if not pe_mode: # real SE
      	pe_mode_for_SE = 1 #
      	read_seq_pe_file = read_seq
      	pe_mode=1
      ## -- End
      index_fragment = 0
      for r in read_seq:
         prev_index_fragment = index_fragment
	 tag_nonunique_NH = 0
	 tag_overlapping_genes = 0
	 flag_aln_not_unique = 0 #
	 flag_ambiguous = 0 #
	 #-- LOOP OVER ALL READS IN INPUT BAM FILE
	 if pe_mode_for_SE:
	 	r = (r, None)
      	 counter_fragment += 1	
         i += 1
         if not pe_mode:
	    # -- SINGLE_END mode
            if not r.aligned:
               notaligned += 1
               #write_to_samout( r, "not_aligned" )
               continue
            try:
               if r.optional_field( "NH" ) > 1:
		  # --- Rescue multimappers in singel-end mode
                  #write_to_samout( r, "alignment_not_unique" )
                  #nonunique += 1
                  continue
            except KeyError:
               pass
            if r.aQual < minaqual:
               lowqual += 1
               #write_to_samout( r, "too_low_aQual" )
               continue
            if stranded != "reverse":
               iv_seq = ( co.ref_iv for co in r.cigar if co.type == "M" )
            else:
               iv_seq = ( invert_strand( co.ref_iv ) for co in r.cigar if co.type == "M" )            
         else:
	    # -- PAIRED-END
            if r[0] is not None and r[0].aligned:
               if stranded != "reverse":
                  iv_seq = ( co.ref_iv for co in r[0].cigar if co.type == "M" )
               else:
                  iv_seq = ( invert_strand( co.ref_iv ) for co in r[0].cigar if co.type == "M" )
            else:
               iv_seq = tuple()
            if r[1] is not None and r[1].aligned:            
               if stranded != "reverse":
                  iv_seq = itertools.chain( iv_seq, 
                     ( invert_strand( co.ref_iv ) for co in r[1].cigar if co.type == "M" ) )
               else:
                  iv_seq = itertools.chain( iv_seq, 
                     ( co.ref_iv for co in r[1].cigar if co.type == "M" ) )
            else:
               if ( r[0] is None ) or not ( r[0].aligned ):
                  #write_to_samout( r, "not_aligned" )
                  notaligned += 1
                  continue         
            try:
               if (( r[0] is not None and r[0].optional_field( "NH" ) > 1 ) or \
                     ( r[1] is not None and r[1].optional_field( "NH" ) > 1 )):
	       	  tag_nonunique_NH = 1
               	  if ( r[0] is not None and r[1] is None ):
			result, fs_genes, fs_exons,dict_read_name_genes_names,ambiguous_tag = is_read_in_gene_interval(r[0], features,dict_read_name_genes_names,ambiguous_tag, exons)
			if result:
				flag_result = 1
				(dict_nonunique, flag_aln_not_unique, dict_read_name_genes_names) = _add_non_unique_counts_per_feature_based_on_read_interval_and_readname(r[0], \
												temp_interval_r0, temp_read_name, fs_genes, dict_nonunique, dict_read_name_genes_names)
			else:
				if len(fs_genes) != 0:
					(dict_nonunique, flag_aln_not_unique, dict_read_name_genes_names) = _add_non_unique_counts_per_feature_based_on_read_interval_and_readname(r[0], \
													temp_interval_r0, temp_read_name, fs_genes, dict_nonunique, dict_read_name_genes_names)
               	  if ( r[0] is None and r[1] is not None ):
			result, fs_genes, fs_exons,dict_read_name_genes_names,ambiguous_tag = is_read_in_gene_interval(r[1], features,dict_read_name_genes_names,ambiguous_tag,exons)
			if result:
				flag_result = 1
				(dict_nonunique, flag_aln_not_unique, dict_read_name_genes_names) = _add_non_unique_counts_per_feature_based_on_read_interval_and_readname(r[1], \
												temp_interval_r1, temp_read_name, fs_genes, dict_nonunique, dict_read_name_genes_names)
			else:
				if len(fs_genes) != 0:
					(dict_nonunique, flag_aln_not_unique, dict_read_name_genes_names) = _add_non_unique_counts_per_feature_based_on_read_interval_and_readname(r[1], \
													temp_interval_r1, temp_read_name, fs_genes, dict_nonunique, dict_read_name_genes_names)
               	  if ( r[0] is not None and r[1] is not None ):
			result1, fs_genes1, fs_exons1,dict_read_name_genes_names,ambiguous_tag = is_read_in_gene_interval(r[0], features,dict_read_name_genes_names,ambiguous_tag,exons)
			result2, fs_genes2, fs_exons2,dict_read_name_genes_names,ambiguous_tag = is_read_in_gene_interval(r[1], features,dict_read_name_genes_names,ambiguous_tag,exons)

		        if len(fs_genes1.intersection(fs_genes2)) > 0:
				fs_genes = fs_genes1.intersection(fs_genes2)
		        elif len(fs_genes1.intersection(fs_genes2))==0:
				fs_genes = fs_genes1.union(fs_genes2)

			if result1 and not result2:
				flag_result = 1
				(dict_nonunique, flag_aln_not_unique, dict_read_name_genes_names) = _add_non_unique_counts_per_feature_based_on_read_interval_and_readname(r[0], \
												temp_interval_r0, temp_read_name, fs_genes, dict_nonunique, dict_read_name_genes_names)
			elif result2 and not result1:
				flag_result = 1
				(dict_nonunique, flag_aln_not_unique, dict_read_name_genes_names) = _add_non_unique_counts_per_feature_based_on_read_interval_and_readname(r[1], \
												temp_interval_r1, temp_read_name, fs_genes, dict_nonunique, dict_read_name_genes_names)
			else:
				if len(fs_genes1) != 0 or len(fs_genes2) != 0:
					flag_result = 1
					if ( ( ((temp_interval_r0 != str(r[0].iv)) or (temp_interval_r1 != str(r[1].iv))) or (temp_read_name != r[0].read.name) ) ):
						(dict_nonunique)= add_non_unique_counts_per_feature(fs_genes, dict_nonunique)
						dict_read_name_genes_names = _populate_read_name_gene_name(dict_read_name_genes_names, fs_genes, r[0].read.name, tag_report_instances_same_multiread_on_same_gene)
						flag_aln_not_unique = 1
                  #write_to_samout( r, "alignment_not_unique" )
	          nonunique += 1

		  if flag_result:
			
			  if r[0] is not None and r[1] is None:		
				non_uniq_read_name = r[0].read.name
			  elif r[0] is None and r[1] is not None:		
				non_uniq_read_name = r[1].read.name
			  elif r[0] is not None and r[1] is not None:		
				non_uniq_read_name= r[0].read.name
			  non_uniq_read_name2 = dict_read_name_genes_names.keys()[0]
			  if flag_aln_not_unique:
				nonunique_nonamb_to_be_rescued += 1
	          	  # -- Re-initialise hash
			  # previous_read_name: read which falls into at least one gene interval
			  # tmp_read_name: the previous read in the bam file
			  # BAM is sorted by read name hence each multimapper will be arranged one after another 
			  if previous_read_name == "NA":
				previous_read_name = non_uniq_read_name

		  	  if non_uniq_read_name != previous_read_name:
				if previous_read_name in dict_read_name_genes_names.keys():
					fs_genes_names = dict_read_name_genes_names[previous_read_name]
					fh_read_names_gene_names.write("%s\t%s\n" % (previous_read_name, "\t".join(list(fs_genes_names)) ))
				previous_read_name = non_uniq_read_name
				tmp_dict = {}
				if non_uniq_read_name in dict_read_name_genes_names.keys():
					#print "non_uniq_read_name IN dict_read_name_genes_names.keys()"
					tmp_dict[non_uniq_read_name] = dict_read_name_genes_names[non_uniq_read_name]
				dict_read_name_genes_names.clear() # only one read stored
				dict_read_name_genes_names = tmp_dict	

		  flag_result = 0
		  flag_aln_not_unique = 0 #
		  (temp_read_name, temp_interval_r0, temp_interval_r1) = initalize_read_name_and_interval(r[0], r[1]) 
		  continue
            # except KeyError:
            except KeyError:
               pass
            if ( r[0] and r[0].aQual < minaqual ) or ( r[1] and r[1].aQual < minaqual ):
               lowqual += 1
               #write_to_samout( r, "too_low_aQual" )
               continue         
          
         try:
	    # --
            if overlap_mode == "union":
               fs = set()
               for iv in iv_seq: # interval from bam file for each fragment
                  if iv.chrom not in features.chrom_vectors:
                     raise UnknownChrom
                  for iv2, fs2 in features[ iv ].steps():
		     	#if debug:
				#print "****Unique_feature %s and feature_interval %s" %(fs2,iv2)	
		        fs = fs.union( fs2 )
			
            elif overlap_mode == "intersection-strict" or overlap_mode == "intersection-nonempty":
               fs = None
               for iv in iv_seq:
                  if iv.chrom not in features.chrom_vectors:
                     raise UnknownChrom
                  for iv2, fs2 in features[ iv ].steps():
                     if len(fs2) > 0 or overlap_mode == "intersection-strict":
                        if fs is None:
                           fs = fs2.copy()
                        else:
                           fs = fs.intersection( fs2 )
            else:
               sys.exit( "Illegal overlap mode." )

	    fs_genes = fs
            if fs_genes is None or len( fs_genes ) == 0:
               #write_to_samout( r, "no_feature" )
               empty += 1
		# ambiguous read count and/or one of the read pair mapping on different gene (potential gene fusion events)...
		# elif len( fs ) > 1:
            elif len( fs_genes ) > 1:
	       ###############################################################
	       ## AMBIGUOUS UNIQUE
	       ###############################################################
	       is_disambiguated = 0
	       if not tag_nonunique_NH:
                  if ( r[0] is not None and r[1] is None ):
			       result, fs_genes, fs_exons,dict_read_name_genes_names_ambiguous, ambiguous_tag = is_read_in_gene_interval(r[0], features, dict_read_name_genes_names_ambiguous, ambiguous_tag,exons)
			       if result:
		       			(dict_gene_unique_counts) = add_unique_counts_per_feature(dict_gene_unique_counts, fs_genes)
	       				is_disambiguated = 1
			       if ambiguous_tag:
			       		(dict_gene_unique_counts_ambiguous) = add_unique_counts_per_feature_ambiguous(fs_genes, dict_gene_unique_counts_ambiguous)
					flag_ambiguous = 1
					# write in the file ambiguous read name gene name data...
					fh_read_names_gene_names_amb_unique.write("%s\t%s\n" % (r[0].read.name, "\t".join(list(fs_genes)) ))
                  if ( r[0] is None and r[1] is not None ):
			       result, fs_genes, fs_exons, dict_read_name_genes_names_ambiguous, ambiguous_tag = is_read_in_gene_interval(r[1], features, dict_read_name_genes_names_ambiguous, ambiguous_tag,exons)
			       if result:
		       			(dict_gene_unique_counts) = add_unique_counts_per_feature(dict_gene_unique_counts, fs_genes)
	       				is_disambiguated = 1
			       if ambiguous_tag:
			       		(dict_gene_unique_counts_ambiguous) = add_unique_counts_per_feature_ambiguous(fs_genes, dict_gene_unique_counts_ambiguous)
					flag_ambiguous = 1
					fh_read_names_gene_names_amb_unique.write("%s\t%s\n" % (r[1].read.name, "\t".join(list(fs_genes)) ))
                  if ( r[0] is not None and r[1] is not None ):
			       result1, fs_genes1, fs_exons, dict_read_name_genes_names_ambiguous, ambiguous_tag1 = is_read_in_gene_interval(r[0], features, dict_read_name_genes_names_ambiguous, ambiguous_tag,exons)
			       result2, fs_genes2, fs_exons, dict_read_name_genes_names_ambiguous, ambiguous_tag2 = is_read_in_gene_interval(r[1], features, dict_read_name_genes_names_ambiguous, ambiguous_tag,exons)
			       if debug:
			       		print "IN UNIQUE DISAMBIGUATION -->r[0].read.name=%s\t%s\t%s\t%s\t%s\n" % (r[0].read.name,result1, result2, fs_genes1, fs_genes2)
			       if len(fs_genes1.intersection(fs_genes2))==1:
					fs_genes = fs_genes1.intersection(fs_genes2)
					(dict_gene_unique_counts) = add_unique_counts_per_feature(dict_gene_unique_counts, fs_genes)
	       				is_disambiguated = 1
			       elif len(fs_genes1.intersection(fs_genes2)) > 1:
					fs_genes = fs_genes1.intersection(fs_genes2)
					(dict_gene_unique_counts_ambiguous) = add_unique_counts_per_feature_ambiguous(fs_genes, dict_gene_unique_counts_ambiguous)
					flag_ambiguous = 1
					fh_read_names_gene_names_amb_unique.write("%s\t%s\n" % (r[0].read.name, "\t".join(list(fs_genes)) ))
			       elif len(fs_genes1.intersection(fs_genes2))==0:
					fs_genes = fs_genes1.union(fs_genes2)
					if (fs_genes1 == set([]) or fs_genes2 == set([])) and len(fs_genes) == 1: 					
						## Disambiguate the uniquely mapped to the single gene it maps on
						(dict_gene_unique_counts) = add_unique_counts_per_feature(dict_gene_unique_counts, fs_genes)
	       					is_disambiguated = 1
					elif (fs_genes1 != set([]) or fs_genes2 != set([])):
						## Add fragment to the RN-GN for ambiguous uniquely mapped based on 
						## union of both fs_genes (fs_genes1 & fs_genes2) > 1
						(dict_gene_unique_counts_ambiguous) = add_unique_counts_per_feature_ambiguous(fs_genes, dict_gene_unique_counts_ambiguous)
						flag_ambiguous = 1
						fh_read_names_gene_names_amb_unique.write("%s\t%s\n" % (r[0].read.name, "\t".join(list(fs_genes)) ))

	       if flag_ambiguous:
			ambiguous += 1
			#write_to_samout( r, "ambiguous[" + '+'.join( fs ) + "]" )
               if is_disambiguated:
			write_to_samout( r, list(fs_genes)[0] )
            else:
	       if debug:
		       #print "DEBUG::CR:: len(fs) <-> 1:: fs = %s" %fs
			pass
               write_to_samout( r, list(fs)[0] )

               rr2 = r[0] if r[0] is not None else r[1]

	       if not tag_nonunique_NH:
			(dict_gene_unique_counts) = add_unique_counts_per_feature(dict_gene_unique_counts, fs_genes)
			
         except UnknownChrom:
            if not pe_mode:
               rr = r 
            else: 
               rr = r[0] if r[0] is not None else r[1]
            if not quiet:
               sys.stderr.write( ( "Warning: Skipping read '%s', because chromosome " +
                  "'%s', to which it has been aligned, did not appear in the GFF file.\n" ) % 
                  ( rr.read.name, iv.chrom ) )

         if i % 100000 == 0 and not quiet:
            sys.stderr.write( "%d sam %s processed.\n" % ( i, "lines " if not pe_mode else "line pairs" ) )

	 flag_ambiguous = 0 ## re-initialise....
	 index_fragment += 1
      #########################
      # This is to store the last read/fragment since it will no pass in previous condition:
      # => if non_uniq_read_name != previous_read_name:
      # -- At same level as the for loop (outside of the for loop) - column: 7
      #fh_read_names_gene_names.close()
      if dict_read_name_genes_names.keys() != []:
	#print "dict_read_name_genes_names passing"
	non_uniq_read_name = dict_read_name_genes_names.keys()[0]
	fs_genes_names = dict_read_name_genes_names[non_uniq_read_name]
	fh_read_names_gene_names.write("%s\t%s\n" % (non_uniq_read_name, "\t".join(list(fs_genes_names)) ))
      # -- 
      fh_read_names_gene_names.close() 
      fh_read_names_gene_names_amb_unique.close()
   ###################################################################################################   
   #except UnboundLocalError:
   except AttributeError:
   #except:
      if not pe_mode:
         sys.stderr.write( "Error occured in %s.\n" % read_seq.get_line_number_string() )
      else:
         sys.stderr.write( "Error occured in %s.\n" % read_seq_pe_file.get_line_number_string() )
      raise

   if not quiet:
      sys.stderr.write( "%d sam %s processed.\n" % ( i, "lines " if not pe_mode else "line pairs" ) )
         
   if samoutfile is not None:
      samoutfile.close()

   if tag_gff == "gene_gff":
	   tuples_genenames_exontag = [(fn, fn) for fn in dict_gene_unique_counts.keys()]
   tuples_genenames_exontag.sort()

   previous_gene_name = "NA"

   for gene_name, fn in tuples_genenames_exontag:
	gene_name = gene_name.strip()
	fn = fn.strip()
	
   	if tag_gff == "gene_gff": #
		if gene_name in dict_gene_unique_counts.keys():
			print "%s\t%i\t%i\t%s" % ( fn, dict_gene_unique_counts[gene_name], dict_nonunique[gene_name],dict_gene_unique_counts_ambiguous[gene_name] )
		else:
			# -- No non-unique reads for that gene_name
			print "%s\t%i\t%i\t%i" % ( fn, dict_gene_unique_counts[gene_name], 0,dict_gene_unique_counts_ambiguous[gene_name] )
		
	# -- Re-initialise gene name
	previous_gene_name = gene_name
			
   print "no_feature\t%d" % empty
   print "ambiguous\t%d" % ambiguous
   print "too_low_aQual\t%d" % lowqual
   print "not_aligned\t%d" % notaligned
   print "alignment_not_unique\t%d" % nonunique
   print "nonunique_nonamb_to_be_rescued:\t%d"  % nonunique_nonamb_to_be_rescued
   ## -- Close output file for non_unique read names --> [genes..]
   #fh_read_names_gene_names.close() 
   #fh_save_gene_unique_reads_total_rescued.close()


############################################################################################################################################################
## -- EXTRACT PASSED PARAMETERS & RUN MAIN FUNCTION: COUNT READS PER GENE
############################################################################################################################################################
def main():
   
   optParser = optparse.OptionParser( 
      
      usage = "%prog [options] sam_file gff_file -m <mode> -s <strandedness> -a <minaqual> -t <featureType> -i <gene_id> -o <samout>",
      
      description=
         "This script takes an alignment file in SAM format and a " +
         "feature file in GFF format and calculates for each feature " +
         "the number of reads mapping to it. See " +
         "http://www-huber.embl.de/users/anders/HTSeq/doc/count.html for details.",
         
      epilog = 
         "Written by Simon Anders (sanders@fs.tum.de), European Molecular Biology " +
         "Laboratory (EMBL). (c) 2010. Released under the terms of the GNU General " +
         "Public License v3. Part of the 'HTSeq' framework, version %s.\n "% HTSeq.__version__  +
         "Extended version written by the peakRescue team to handle ambiguously mapped reads." )
         
   optParser.add_option( "-m", "--mode", type="choice", dest="mode",
      choices = ( "union", "intersection-strict", "intersection-nonempty" ), 
      default = "union", help = "mode to handle reads overlapping more than one feature" +
         "(choices: union, intersection-strict, intersection-nonempty; default: union)" )
         
   optParser.add_option( "-s", "--stranded", type="choice", dest="stranded",
      choices = ( "yes", "no", "reverse" ), default = "yes",
      help = "whether the data is from a strand-specific assay. Specify 'yes', " +
         "'no', or 'reverse' (default: yes). " +
         "'reverse' means 'yes' with reversed strand interpretation" )
      
   optParser.add_option( "-a", "--minaqual", type="int", dest="minaqual",
      default = 0,
      help = "skip all reads with alignment quality lower than the given " +
         "minimum value (default: 0)" )
      
   optParser.add_option( "-t", "--type", type="string", dest="featuretype",
      default = "exon", help = "feature type (3rd column in GFF file) to be used, " +
         "all features of other type are ignored (default, suitable for Ensembl " +
         "GTF files: exon)" )
         
   optParser.add_option( "-i", "--idattr", type="string", dest="idattr",
      default = "gene_id", help = "GFF attribute to be used as feature ID (default, " +
      "suitable for Ensembl GTF files: gene_id)" )

   optParser.add_option( "-o", "--samout", type="string", dest="samout",
      default = "", help = "write out all SAM alignment records into an output " +
      "SAM file called SAMOUT, annotating each line with its feature assignment " +
      "(as an optional field with tag 'XF')" )

   optParser.add_option( "-q", "--quiet", action="store_true", dest="quiet",
      help = "suppress progress report and warnings" )

   if len( sys.argv ) == 1:
      optParser.print_help()
      sys.exit(1)

   (opts, args) = optParser.parse_args()

   if len( args ) != 4:
      sys.stderr.write( sys.argv[0] + ": Error: Please provide (originally two) now 4 arguments.\n" )
      sys.stderr.write( "  Call with '-h' to get usage information.\n" )
      sys.exit( 1 )

   warnings.showwarning = my_showwarning
   try:
      filename_read_names_gene_names = args[2]
      filename_read_names_gene_names_amb_unique = args[3]
	
   except:
      sys.stderr.write( "Error: %s\n" % str( sys.exc_info()[1] ) )
      sys.stderr.write( "[Exception type: %s, raised in %s:%d]\n" % 
         ( sys.exc_info()[1].__class__.__name__, 
           os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
           traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
      sys.exit( 1 )

   # -- Added outside the try-except block to get a more meaningful error msg
   count_reads_in_features( args[0], args[1], opts.stranded, 
   	opts.mode, opts.featuretype, opts.idattr, opts.quiet, opts.minaqual,
        opts.samout, filename_read_names_gene_names,filename_read_names_gene_names_amb_unique)

def my_showwarning( message, category, filename, lineno = None, line = None ):
   sys.stderr.write( "Warning: %s\n" % message )

if __name__ == "__main__":
   main()

