import sys, optparse, itertools, warnings, traceback, os.path

import HTSeq

# -- Global variable for debug # added CR
debug = 0
#debug = 1

class UnknownChrom( Exception ):
   pass

def weighting_non_unique_reads_with_max_exon_count(fh_read_names_gene_names, dict_gene_exons_unique_counts):
	"""
	Weighting the non-unique read counts per gene based on the proportion of unique reads mapped to each gene 
	(based on max number of unique reads over all exons of a gene)

	"""
	#print "--- WEIGHTING METHOD:  dict_gene_exons_unique_counts.keys()[0] = **%s**" % dict_gene_exons_unique_counts.keys()[0]
	line = fh_read_names_gene_names.readline()
	#dict_gene_max_of_unique_counts_over_all_exons = {}
	dict_gene_non_unique_proportions = {}
	dict_gene_non_unique_uniform_proportions = {}
	while line:
		dict_gene_max_of_unique_counts_over_all_exons = {}
		#dict_gene_non_unique_proportions = {}
		# -- 
		fields = line.split("\t")
		read_name = fields[0]
		list_gene_names = fields[1:]
		#print "read_gene_names-->%s" %(line) 
		sum_max_unique_counts_exons_all_genes = 0
		for gene_name in list_gene_names:
			gene_name = gene_name.strip()
			if gene_name in dict_gene_exons_unique_counts.keys():
				#print  "Example gene: %s\tExons: %s" % (gene_name, dict_gene_exons_unique_counts[gene_name].keys())
				exon_name = dict_gene_exons_unique_counts[gene_name].keys()[0]
				#print  "Example [dict_gene_exons_unique_counts[gene_name][exon] = : %s" % (dict_gene_exons_unique_counts[gene_name][exon_name])
				list_exons_names = dict_gene_exons_unique_counts[gene_name].keys()
				max_unique_count_per_gene = max([int(dict_gene_exons_unique_counts[gene_name][exon]) for exon in list_exons_names])
				sum_max_unique_counts_exons_all_genes += max_unique_count_per_gene 
				dict_gene_max_of_unique_counts_over_all_exons[gene_name] = max_unique_count_per_gene
				#dict_gene_non_unique_proportions[gene_name] = 0
				#print "--- gene_name: %s \tMAX: %i\tSUM: %i" % (gene_name, max_unique_count_per_gene, sum_max_unique_counts_exons_all_genes)
		if dict_gene_max_of_unique_counts_over_all_exons != {}:
			for gene_name in dict_gene_max_of_unique_counts_over_all_exons.keys():
				#print "gene_name-->%s\t sum: %s" % ( gene_name, sum_max_unique_counts_exons_all_genes)
				proportion = float(dict_gene_max_of_unique_counts_over_all_exons[gene_name]) / sum_max_unique_counts_exons_all_genes
				#print "proportion %1.2f" %proportion
				if gene_name in dict_gene_non_unique_proportions.keys():
					dict_gene_non_unique_proportions[gene_name] += proportion
					dict_gene_non_unique_uniform_proportions[gene_name] += 1./len(list_gene_names)
				else:
					dict_gene_non_unique_proportions[gene_name] = proportion
					dict_gene_non_unique_uniform_proportions[gene_name] = 1./len(list_gene_names)
		line = fh_read_names_gene_names.readline()
	return (dict_gene_non_unique_proportions, dict_gene_non_unique_uniform_proportions)


def set_exons_to_set_genes(fs_exons):
	fs_genes = set()
	
	if not (fs_exons is None or len(fs_exons) == 0):
		for fs_exon in list(fs_exons):
			#print "Debug: fs_exon_loop1_exonset = ", fs_exon
			fs_gene_name = fs_exon.split("_")[1]
			fs_gene=set([fs_gene_name])
			#print "Debug: fs_gene_name = ", fs_gene_name
			fs_genes = fs_genes.union(fs_gene)
			
	#print "Debug: fs_exon_loop_genes = ", fs_genes
	return fs_genes

def is_read_in_gene_interval(readObject, features, dict_read_name_genes_names):
	"""
	Method called in cases of non-unique fragments based on NH tag.
	fs_genes contains sets of exons.
	Read parameter is a read (r[0] or r[1])...

	Re-initialisation of the hash dict_read_name_genes_names outside method
	at the end of condition over multi-mapped reads (> NH:i:1) to avoid
	memory issues regarding recording of many reads (keys).
	
	if the condition is true then return 1 else return 0 # int(start) >= int(iv_start) and int(end) <= int(iv_end) and chr == iv_chr
	 
	 condition1-->r[0] is not None and r[1] is None
	 readname SRR001362.3657620 readObject.iv chr18:[75161453,75161709)/- ----Debug: fs_interval chr18:[75161453,75161454)/. fs_exons --> set(['ENSMUSG00000062328'])
	 readname SRR001362.3657620 readObject.iv chr18:[75161453,75161709)/- ----Debug: fs_interval chr18:[75161454,75161533)/. fs_exons --> set(['ENSMUSG00000062328'])
	 readname SRR001362.3657620 readObject.iv chr18:[75161453,75161709)/- ----Debug: fs_interval chr18:[75161533,75161598)/. fs_exons --> set(['ENSMUSG00000064844', 'ENSMUSG00000062328'])
	 readname SRR001362.3657620 readObject.iv chr18:[75161453,75161709)/- ----Debug: fs_interval chr18:[75161598,75161677)/. fs_exons --> set(['ENSMUSG00000064844', 'ENSMUSG00000062328'])
	 readname SRR001362.3657620 readObject.iv chr18:[75161453,75161709)/- ----Debug: fs_interval chr18:[75161677,75161709)/. fs_exons --> set(['ENSMUSG00000064844', 'ENSMUSG00000062328'])

	"""
	if debug:
		print "=>DEBUG1:: readName=%s\tReadLocation:: chrom=%s\tstart=%s\tend=%s" % (readObject.read.name, readObject.iv.chrom, readObject.iv.start, readObject.iv.end)
	fs_exons = set()
	fs_genes = set()
	##### return (0, fs_genes, fs_exons, dict_read_name_genes_names)
	#Check if read is ambiguous i.e overlaps on two genes ::SB
	genes_ge2_tag=0
	fs_exons_temp= set()
	for iv3_temp, fs_exon_temp in features[ readObject.iv ].steps():
		fs_exons_temp = fs_exons_temp.union( fs_exon_temp )
		if len(fs_exons_temp) > 1:
			genes_ge2_tag=1
			if debug:
				print " genes_ge2 fs_exons", fs_exons_temp
	
			break
	# check if read intervaloverlapswith gene if not return 0
	if  len(fs_exons_temp) == 0:
		if debug:
			print "No gene in read interval - exit iv3_temp %s fs_exon_temp %s" %(iv3_temp,fs_exon_temp)
		return (0, fs_genes, fs_exons, dict_read_name_genes_names)
		#else:
		#	print " genes_lt2 fs_exons", fs_exons_temp
	
	#print "DEBUG::CR:: dir(readObject) = %s" % dir(readObject)
	#print "=>DEBUG:: readName=%s\tdir(readObj)=%s\tdir(readObj.iv)=%s" % (readObject.read.name, dir(readObject), dir(readObject.iv))
	if debug:
		print "=>DEBUG2:: readName=%s\tReadLocation:: chrom=%s\tstart=%s\tend=%s" % (readObject.read.name, readObject.iv.chrom, readObject.iv.start, readObject.iv.end)
	#print "start=%s\tstart_as=%s" % (readObject.iv.start, readObject.iv.start_as_pos)
	#print "end=%s\tend_as=%s" % (readObject.iv.end, readObject.iv.end_as_pos)
	threshold_split_region = 50
	readObject.cigar.reverse()
	tag_split_read = 0
	read_subregions = []
	regions_iv = []
	#if read overlaps on two genes then check if it is due to split mapping 
	if genes_ge2_tag:
				#chek if the region between split ends of a read maps on only one gene.. then assign read to gene G1  (eg.,
				#----====---------====---- G1
				#      <-----N------>
				#---===================---- G2
				# if ends of the read maps on signle gene and N maps on one gene assign read to G1
				# -====--------------====---- G1
				#    <-------N-------->	   Read
				#--------=========----------- G2
				#--intron , ==exon

		for cigar_region in readObject.cigar:
			#print "==>DEBUG:: CigarType=%s\tCigarSize=%s\tQeryFrom=%s\tQueryTo=%s\t%s\n" % ( cigar_region.type, cigar_region.size, cigar_region.query_from, cigar_region.query_to, dir(readObject.cigar))
			if cigar_region.type == "N" and int(cigar_region.size) > threshold_split_region:
				#print "Split read..."
				#print "==DEBUG:: "
				tag_split_read = 1
				

			if cigar_region.type == "M":
				subregion_start = readObject.iv.start + cigar_region.query_from
				subregion_end = subregion_start + cigar_region.query_to
				subregion_iv = HTSeq.GenomicInterval( readObject.iv.chrom, subregion_start, subregion_end, readObject.iv.strand )
				#print "SB redObject.cigar--:%s read_name-->%s\tsubregion_iv: %s" %(readObject.cigar,readObject.read.name, subregion_iv)
				read_subregions.append(subregion_iv)
		if tag_split_read:
				regions_iv = read_subregions
		else:
				## It's not a split read and the read maps on 2 genes on the same location (i.e. ambiguous) 
				return (0, fs_genes, fs_exons, dict_read_name_genes_names)
				#regions_iv = [readObject.iv]
	else:
		if debug:
			print "DEBUG:: single gene"
		regions_iv = [readObject.iv]
	#print "DEBUG::CR:: readObject.cigar = **%s**" % readObject.cigar
	#print "DEBUG::CR:: str(readObject.cigar) = **%s**" % str(readObject.cigar)
	#readObject.cigar.reverse()
	#obj= readObject.cigar.pop()
	#print "DEBUG::CR::***obj=%s\tType=%s\tDIR=%s***" % (obj, type(obj), dir(obj))
	#print "DEBUG::CR::***obj.type=%s\tobj.size=%s\n" % (obj.type, obj.size)
	fs_exons = set()
	#SB counter added to check the number of gene intervals the read overlaps 
	counter=0
	read_name = readObject.read.name
	## -- CR: start commented
	##for iv2, fs2 in features[ iv ].steps():
	#(chr,start_end)=str(readObject.iv).split(":")
	#(start, end_strand)= start_end.split(",")
	#start = start[1:]; end = end_strand.split(")")[0]
	## -- CR: end commented
	##print "chr %s--start %s--end %s" %(chr,start,end)
	results_subregions  = []
	for rObj_iv in regions_iv: ## added_CR
		(chr,start_end)=str(rObj_iv).split(":") ## added_CR
		(start, end_strand)= start_end.split(",") ## added_CR
		start = start[1:]; end = end_strand.split(")")[0] ## added_CR
		# for iv3, fs_exon in features[ readObject.iv ].steps():
		for iv3, fs_exon in features[ rObj_iv ].steps(): ## added_CR
			counter +=1
			(iv_chr, iv_start_end)=str(iv3).split(":")
			(iv_start, iv_end_strand)= iv_start_end.split(",")
			iv_start = iv_start[1:]; iv_end = iv_end_strand.split(")")[0]
			#print "counter %s iv_chr %s--iv_start %s--iv_end %s" %(counter, iv_chr,iv_start,iv_end)
			#if int(start) >= int(iv_start) and int(end) <= int(iv_end) and chr == iv_chr: 	
			#chr18:[75161453,75161709)/+
			#print  "int(start) %s >= int(iv_start) %s  int(end) %s <= int(iv_end) %s "%(start,iv_start,end,iv_end)
			fs_exons = fs_exons.union( fs_exon )
			#print "readname_test %s readObject.iv %s ----Debug: fs_interval %s fs_exon ### %s fs_exons-->%s "%(read_name,readObject.iv,iv3, fs_exon,fs_exons)
			#print " readname %s readObject.iv %s ----Debug: fs_interval %s fs_exon ### %s fs_exons-->%s "%(read_name, rObj_iv, iv3, fs_exon,fs_exons)
			# -- @TODO: ADD A TAG TO SPECIFY WHICH GFF TYPE IS USED (e.g. GENE OR EXON)
			# -- HERE WE COMMENT THE FOLLOWING LINE AS GENE GFF IS USED AND ADD LINE TO ASSIGN fs_exons to fs_genes
			# fs_genes = set_exons_to_set_genes(fs_exons)
			fs_genes = fs_exons
			if debug:
				print "fs_genes-->%s" %(fs_genes)
			if not ((fs_genes is None or len( fs_genes ) == 0 ) and \
				(fs_exons is None or len(fs_exons) == 0)):
				#print " ----Debug: fs_genes = %s", fs_genes
				#print " ----Debug: fs_exons = %s", fs_exons
				if read_name in dict_read_name_genes_names.keys(): 
					# Second read (r[1]) of a given fragment is encountered 
					fs_genes_first_read = dict_read_name_genes_names[read_name]
					dict_read_name_genes_names[read_name] = fs_genes_first_read.union(fs_genes)
				else:
					# First read (r[0]) of a given fragment is encountered 
						dict_read_name_genes_names[read_name] = fs_genes
						#print "dict_with_gene_name---> %s\t%s\n" % (read_name, dict_read_name_genes_names[read_name])
			#SB if condition at the end of for loop to avoid null values in fs_genes, fs_exons, dict_read_name_genes_names
			if genes_ge2_tag:
				if int(start) >= int(iv_start) and int(end) <= int(iv_end) and chr == iv_chr and genes_ge2_tag: 	
					# return (1, fs_genes, fs_exons, dict_read_name_genes_names)
					results_subregions.append((1, fs_genes, fs_exons, dict_read_name_genes_names))
			else:
				results_subregions.append((1, fs_genes, fs_exons, dict_read_name_genes_names))
		#print "return value %s fs_genes", fs_genes
	#return (0, fs_genes, fs_exons, dict_read_name_genes_names)
		#print "readname_test %s readObject.iv %s ----Debug: fs_interval %s fs_exon ### %s fs_exons-->%s "%(read_name,readObject.iv,iv3, fs_exon,fs_exons)
	if results_subregions==[]:
		#results_subregions.append((0, fs_genes, fs_exons, dict_read_name_genes_names))
		#return (0, results_subregions)
		return (0, fs_genes, fs_exons, dict_read_name_genes_names)
	else:
		if tag_split_read: #(and gene_ge2_tag)
			list_subregions_in_genes = [tag_subregion_in_intervals for tag_subregion_in_intervals, fs_genes, fs_exons, dict_read_name_genes_names in results_subregions if tag_subregion_in_intervals==1]
			if debug:	
				print "len(list_subregions_in_genes)=%s len(regions_iv)=%s len(fs_genes)= %s" %(len(list_subregions_in_genes), len(regions_iv), len(fs_genes))
			if len(list_subregions_in_genes) == len(regions_iv):
				#return (1, results_subregions)
				if debug:
					print "DEBUG::CR:: subregions are all present in gene interval!"
				return (1, fs_genes, fs_exons, dict_read_name_genes_names)
		#else:
		#	#return (0, results_subregions)
		#	return (0, fs_genes, fs_exons, dict_read_name_genes_names)
	
	if len(fs_genes) > 1:
		return (0, fs_genes, fs_exons, dict_read_name_genes_names)
	elif len(fs_genes) == 1:
		return (1, fs_genes, fs_exons, dict_read_name_genes_names)
	# new condition added to take care of ambiguous reads mapping on two genes
	else:
		return (0, fs_genes, fs_exons, dict_read_name_genes_names)

def count_and_rescue_unique_reads(tag_gff, counts, dict_gene_unique_counts, fs, dict_gene_exons_unique_counts):
	"""
	Count unique reads and rescue unique reads from multimapped reads which map uniquely to the 
	transcriptome as opposed to those reads mappped to other non-coding 
	regions of the genome.
	count_and_rescue_unique_read_to_transcriptome_from_multimapped_reads(dict_gene_unique_counts, fs)
	
	Tophat with option transcriptome-only will map reads to trnascriptome so NH:i:* (multimapping) would 
	be with respect to transcriptome. No need to follow the rescue step.
	
	Fragment is unique (map on unique location)

	"""
	if debug:
		#print "AAA_DEBUG::CR::IN_count_and_rescue_unique_reads:: tag_gff=%s\tcounts=%s\tdict_gene_unique_counts=%s\tfs=%s"%(tag_gff, counts, dict_gene_unique_counts, fs)
		pass
	#if len( fs ) >= 1: # exon
		#print "debug: ", dir(fs)
	# @TODO ADD TAG FOR GFF TYPE !!!!!!!!!!!!!!!!
	if tag_gff == "exon_gff":
		feature_name1 = list(fs)[0].split("_")[1] # gene name (same gene for all exons)
		for fs_exon in list(fs):
			counts[fs_exon] += 1
			# -- Hash of hash Gene --> exons counts
			if feature_name1 in dict_gene_exons_unique_counts.keys():
				if fs_exon in dict_gene_exons_unique_counts[feature_name1].keys(): 
					dict_gene_exons_unique_counts[feature_name1][fs_exon] += 1
				else: 
					dict_gene_exons_unique_counts[feature_name1][fs_exon] = 1
			else:
				dict_gene_exons_unique_counts[feature_name1] = {}
				dict_gene_exons_unique_counts[feature_name1][fs_exon] = 1
	elif tag_gff == "gene_gff":
		feature_name1 = list(fs)[0] # IF TAG GFF TYPE IS GENE TYPE
		#print "DEBUG::GENE_GFF TYPE:: feature_name1 = %s" % feature_name1
	
	# feature_name1 = list(fs)[0].split("_")[1] # gene name (same gene for all exons)
	#feature_name2 = list(fs)[1].split("_")[1] # 
	#if feature_name1 == feature_name2:
		# split read over two exons on same gene
		 #counts[ list(fs)[0] ] += 1 # exon
		 #counts[ list(fs)[1] ] += 1 # 
		 # print "split reads --> %s" % r[0].read.name
	if feature_name1 in dict_gene_unique_counts.keys():
		dict_gene_unique_counts[feature_name1] += 1 
		if debug: ##  feature_name1=="ENSMUSG00000069392"
			#print "ZZZZ_1DEBUG::feature_name1 = %s" % feature_name1
			pass
	elif feature_name1 not in dict_gene_unique_counts.keys():
		dict_gene_unique_counts[feature_name1] = 1 
		if debug: ## feature_name1=="ENSMUSG00000069392":
			#print "ZZZZ_2DEBUG::feature_name1 = %s" % feature_name1
			pass
	#elif len( fs ) == 1:
	#	counts[ list(fs)[0] ] += 1 # exon
	#	feature_name = list(fs)[0].split("_")[1] # gene name
	#	if feature_name in dict_gene_unique_counts.keys():
	#		dict_gene_unique_counts[feature_name] += 1 
	#	elif feature_name not in dict_gene_unique_counts.keys():
	#		dict_gene_unique_counts[feature_name] = 1 
	if debug:
		#print "BBB_DEBUG::CR::IN_count_and_rescue_unique_reads:: tag_gff=%s\tcounts=%s\tdict_gene_unique_counts=%s\tfs=%s"%(tag_gff, counts, dict_gene_unique_counts, fs)
		pass
	return (dict_gene_unique_counts, counts, dict_gene_exons_unique_counts)


def add_non_unique_counts_per_feature(fs_genes,fs_exons, dict_nonunique, dict_non_unique_exons):
	#print "I am in add_non_unique_counts_per_feature"
	#print " ADD----Debug: fs_genes = **%s**" % fs_genes
	#print " ADD----Debug: fs_exons = **%s**" % fs_exons
	#print "START OF METHOD add_non_unique_counts_per_feature"
	#print "fs_genes = %s" % fs_genes
	#print "fs_exons = %s" % fs_exons
	for gene_feature in list(fs_genes):
		if debug :
			print "B4 ----Debug: gene_feature = **%s**\tdict_nonunique[] = %i" %(gene_feature, dict_nonunique[ gene_feature ])
		dict_nonunique[ gene_feature ] += 1
		if debug :
			print "AFTER ----Debug: gene_feature = **%s**\tdict_nonunique[] = %i" %(gene_feature, dict_nonunique[ gene_feature ])
	for exon_feature in list(fs_exons):
		#print " ADD2 ----Debug: exon_feature = **%s**" % exon_feature
		dict_non_unique_exons[ exon_feature ] += 1
	#print "gene_feature= %s" % gene_feature
	#print "END OF METHOD add_non_unique_counts_per_feature"
	return (dict_nonunique, dict_non_unique_exons)

def initialise_counts_per_feature(dict_features, feature_name):
	"""
	Initialisation: to do: optimisation to avoid re-assigning count to same gene (with multiple exons)
			
	"""
	dict_features[feature_name] = 0
	return dict_features

def initalize_read_name_and_interval(read_object0, read_object1):
	temp_interval0 = None
	temp_interval1 = None
	if (read_object0 is not None and read_object1 is None):
		temp_interval0 = read_object0.iv
		temp_read_name = read_object0.read.name
		pass
	elif (read_object0 is None and read_object1 is not None):
		temp_interval1 = read_object1.iv
		temp_read_name = read_object1.read.name
		pass
	elif (read_object0 is not None and read_object1 is not None):
		temp_interval0 = read_object0.iv
		temp_interval1 = read_object1.iv
		temp_read_name = read_object0.read.name
		pass
	return (temp_read_name, temp_interval0, temp_interval1)

def invert_strand( iv ):
   iv2 = iv.copy()
   if iv2.strand == "+":
      iv2.strand = "-"
   elif iv2.strand == "-":
      iv2.strand = "+"
   else:
      raise ValueError, "Illegal strand"
   return iv2

def count_reads_in_features( sam_filename, gff_filename, outfile_table_genenames_unique_counts_total_proportion_rescued, stranded, overlap_mode, feature_type, id_attribute, quiet, minaqual, samout, filename_read_names_gene_names):
   # -- Commented as we commented the writing in that file further down      
   fh_read_names_gene_names = open(filename_read_names_gene_names, 'w')
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
   ##SB flag to set debug mode on 
   ## CR: replaced the variable 'debug' with a global variable (of same name) for use outside the scope of this method
   # debug=1 ## CR: please see new global variable set above (Line 6)
   ## Hash table to store unique reads per exon (if modified GTF)
   counts = {}
   ## Hash table to store original non unique reads per gene (without 
   dict_nonunique = {}
   ## Hash table to store non unique reads per exon (if modified GTF)
   dict_non_unique_exons = {}
   ## Hash table to store all unique reads as per original GTF
   dict_gene_unique_counts = {}
   ## Hash table to store all non-unique reads including shared reads 
   ## (either split reads or read pair matching on two distinct exons, same gene)
   dict_gene_nonunique_counts = {}
   ## Hash of hash table to store exon counts for each gene
   dict_gene_exons_unique_counts = {}
   ## Hash to store the non-unique read-names as key and genes names as velues (fragments)
   dict_read_name_genes_names = {}
   ## GFF type to specify whether it contains GENE or EXONS based information 
   tag_gff = "gene_gff" #"exon_gff" ## SEE REMAINING LOCATIONS TO ADD TAG @TODO
   # Try to open samfile to fail early in case it is not there
   if sam_filename != "-":
      open( sam_filename ).close() 
      
   gff = HTSeq.GFF_Reader( gff_filename )   
   i = 0
   try:
      for f in gff:
         if f.type == feature_type:
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
	    ##added by CR
	    # dict_nonunique[ f.attr[ id_attribute ] ] = 0
	    # -- Initialisation 
	    feature_name = f.attr[ id_attribute ]
	    # -- Added tag_gff for GFF type
	    #if "_" in feature_name:
	    if tag_gff == "exon_gff":
		# GTF contains modified exon names with underscore before gene name 120_ENSSG*****0001
		# Exon
		dict_non_unique_exons = initialise_counts_per_feature(dict_non_unique_exons, feature_name)
		# Gene
		dict_nonunique = initialise_counts_per_feature(dict_nonunique, feature_name.split("_")[1])
		dict_gene_unique_counts = initialise_counts_per_feature(dict_gene_unique_counts, feature_name.split("_")[1])
		dict_gene_nonunique_counts = initialise_counts_per_feature(dict_gene_nonunique_counts, feature_name.split("_")[1])
	    elif tag_gff == "gene_gff":
		# Original GTF (genes) 
		dict_non_unique_exons = initialise_counts_per_feature(dict_non_unique_exons, feature_name)
		dict_nonunique = initialise_counts_per_feature(dict_nonunique, feature_name)
		dict_gene_unique_counts = initialise_counts_per_feature(dict_gene_unique_counts, feature_name)
		dict_gene_nonunique_counts = initialise_counts_per_feature(dict_gene_nonunique_counts, feature_name)
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

   try:
      if pe_mode:
         read_seq_pe_file = read_seq
         read_seq = HTSeq.pair_SAM_alignments( read_seq )
      empty = 0
      ambiguous = 0
      notaligned = 0
      lowqual = 0
      nonunique = 0
      #added by SB
      temp_read_name="NA"
      previous_read_name="NA"
      temp_interval_r0="NA"
      temp_interval_r1="NA"
      ## added by CR	
      nonunique2 = 0
      counter_fragment = 0	
      flag_result = 0
      #added by SB
      i = 0   
      pe_mode_for_SE = 0
      ## -- Added pe_mode on for SE files so that MU reads will be accounted for
      if not pe_mode: # real SE
      	pe_mode_for_SE = 1 ## Added by us
      	read_seq_pe_file = read_seq
      	pe_mode=1
      ## -- End
      for r in read_seq:
	 #print "LOOP OVER READS"
	 if pe_mode_for_SE:
	 	r = (r, None)
      	 counter_fragment += 1	
         i += 1
         if not pe_mode:
	    # -- SINGLE_END mode
            if not r.aligned:
               notaligned += 1
               write_to_samout( r, "not_aligned" )
               continue
            try:
               if r.optional_field( "NH" ) > 1:
		  # --- Rescue multimappers in singel-end mode
                  write_to_samout( r, "alignment_not_unique" )
                  nonunique += 1
                  continue
            except KeyError:
               pass
            if r.aQual < minaqual:
               lowqual += 1
               write_to_samout( r, "too_low_aQual" )
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
                  write_to_samout( r, "not_aligned" )
                  notaligned += 1
                  continue         
            try:
	       tag_nonunique_NH = 0
               if (( r[0] is not None and r[0].optional_field( "NH" ) > 1 ) or \
                     ( r[1] is not None and r[1].optional_field( "NH" ) > 1 )):
	       	  tag_nonunique_NH = 1
                  nonunique += 1
		  #print "nonunique = %i" % nonunique
		  #print "Multimapper case... (NH>1)"
               	  if ( r[0] is not None and r[1] is None ):
			if debug:
				print "condition1-->r[0] is not None and r[1] is None"
			result, fs_genes, fs_exons,dict_read_name_genes_names = is_read_in_gene_interval(r[0], features,dict_read_name_genes_names)
			if debug:
				print "result flag after first NH:1 %s" % result ## CR: corrected format %s -> % result 
			if result:	
				flag_result = 1
				if debug:
					print "passed condition if ( r[0] is not None and r[1] is None ): %s --flagresult-->%s " %(result,flag_result)
					# condition checks if the read is PCR duplicated if starts at same location SB
					print "readName=%s\ttemp_interval_r0 %s , r[0].iv -->%s" %(r[0].read.name, temp_interval_r0, r[0].iv)
				if ( temp_interval_r0 is not r[0].iv):
					#print "up passed ( temp_interval_r0 is not r[0].iv)--flag_result--%s" %(flag_result)
					(dict_nonunique, dict_non_unique_exons) = add_non_unique_counts_per_feature(fs_genes, fs_exons, dict_nonunique, dict_non_unique_exons )
					#print "down passed ( temp_interval_r0 is not r[0].iv)--flag_result--%s" %(flag_result)
					#print "add_non_unique_counts_per_feature passed ( temp_interval_r0 is not r[0].iv)--flag_result--%s" %(flag_result)
				else:
					if (temp_read_name != r[0].read.name):
						(dict_nonunique, dict_non_unique_exons) = add_non_unique_counts_per_feature(fs_genes, fs_exons, dict_nonunique, dict_non_unique_exons )
					#print "passed (temp_read_name != r[0].read.name)--flag_result--%s" %(flag_result)
		  if debug:
			print " flag value after condition -->1 flag_result--> %s" %(flag_result)
				
               	  if ( r[0] is None and r[1] is not None ):
			if debug:
				print "(condition2--> r[0] is None and r[1] is not None )"
			result, fs_genes, fs_exons,dict_read_name_genes_names = is_read_in_gene_interval(r[1], features,dict_read_name_genes_names)
			if result:
				flag_result = 1
				if debug:
					print "read in gene interval r[0] is None and r[1] is not None", result
 				if ( temp_interval_r1 is not r[1].iv):
					(dict_nonunique, dict_non_unique_exons)= add_non_unique_counts_per_feature(fs_genes, fs_exons, dict_nonunique, dict_non_unique_exons )
				else:
					if (temp_read_name != r[1].read.name):
						(dict_nonunique, dict_non_unique_exons)= add_non_unique_counts_per_feature(fs_genes, fs_exons, dict_nonunique, dict_non_unique_exons )
		  	if debug:
				print " flag value after condition -->2 flag_result --->%s" %(flag_result)
					
               	  if ( r[0] is not None and r[1] is not None ):
			if debug:
				print "(condition3--> r[0] is not None and r[1] is not None )"
			result1, fs_genes1, fs_exons1,dict_read_name_genes_names = is_read_in_gene_interval(r[0], features,dict_read_name_genes_names)
			result2, fs_genes2, fs_exons2,dict_read_name_genes_names = is_read_in_gene_interval(r[1], features,dict_read_name_genes_names)

			if result1 and not result2:
				flag_result = 1
				if debug:
					print "read in gene interval r[0] is not None and r[1] is not None", result1 
				if ( temp_interval_r0 is not r[0].iv):
					(dict_nonunique, dict_non_unique_exons)= add_non_unique_counts_per_feature(fs_genes1, fs_exons1, dict_nonunique, dict_non_unique_exons )
				else:
					if (temp_read_name != r[0].read.name):
						(dict_nonunique, dict_non_unique_exons)= add_non_unique_counts_per_feature(fs_genes1, fs_exons1, dict_nonunique, dict_non_unique_exons )

			elif result2 and not result1:
				flag_result = 1
				
				if ( temp_interval_r1 is not r[1].iv):
					(dict_nonunique, dict_non_unique_exons)= add_non_unique_counts_per_feature(fs_genes2, fs_exons2, dict_nonunique, dict_non_unique_exons )
				else:
					if (temp_read_name != r[1].read.name):
						(dict_nonunique, dict_non_unique_exons)= add_non_unique_counts_per_feature(fs_genes2, fs_exons2, dict_nonunique, dict_non_unique_exons )
			elif result1 and result2:
				flag_result = 1
				if ((temp_interval_r0 is not r[0].iv ) and \
					( temp_interval_r1 is not r[1].iv) ):
					if list(fs_genes1)[0] !=  list(fs_genes2)[0]:
						#TO count read1 and read2 to respective genes..
						(dict_nonunique, dict_non_unique_exons)= add_non_unique_counts_per_feature(fs_genes1, fs_exons1, dict_nonunique, dict_non_unique_exons )
						(dict_nonunique, dict_non_unique_exons)= add_non_unique_counts_per_feature(fs_genes2, fs_exons2, dict_nonunique, dict_non_unique_exons )
					else:
						#TO count fragment (read1 and read2) to respective gene
						(dict_nonunique, dict_non_unique_exons)= add_non_unique_counts_per_feature(fs_genes1, fs_exons1, dict_nonunique, dict_non_unique_exons )

				else: 
					# Interval is the same: we count only if read name (fragment name) is different and 
					# both reads are mapping on the same gene. 
					# 
					if (temp_read_name != r[0].read.name):					
						if list(fs_genes1)[0] !=  list(fs_genes2)[0]:
						# To avoid the reads mapping overlapping interavals of two genes...i.e ambiguous
							pass
						else:
							(dict_nonunique, dict_non_unique_exons)= add_non_unique_counts_per_feature(fs_genes1, fs_exons1, dict_nonunique, dict_non_unique_exons )
		        if debug:
				print " flag value after condition -->3 passed flag_restlt --->%s" %(flag_result)				
		  		print "read is multimapper %s" %(flag_result) # added SB

                  write_to_samout( r, "alignment_not_unique" )
		  if flag_result:
			  if debug:
				  #print "read is multimapper %s" %(flag_result) # added SB
				  pass
			  #print "DEBUG: KEYS NUMBER: ", len(dict_read_name_genes_names)
			  non_uniq_read_name2 = dict_read_name_genes_names.keys()[0]
			  if r[0] is not None and r[1] is None:		
				non_uniq_read_name = r[0].read.name
			  elif r[0] is None and r[1] is not None:		
				non_uniq_read_name = r[1].read.name
			  elif r[0] is not None and r[1] is not None:		
				non_uniq_read_name= r[0].read.name
			  #print "--> non_uniq_read_name = ", non_uniq_read_name 
	          	  # -- Re-initialise hash
			  # previous_read_name: read which falls into at least one gene interval
			  # tmp_read_name: the previous read in the bam file
			  # BAM is sorted by read name hence each multimapper will be arranged one after another 
			  if previous_read_name == "NA":
				previous_read_name = non_uniq_read_name
		  	  if non_uniq_read_name != previous_read_name:
				fs_genes_names = dict_read_name_genes_names[previous_read_name]
				fh_read_names_gene_names.write("%s\t%s\n" % (previous_read_name, "\t".join(list(fs_genes_names)) ))
				previous_read_name = non_uniq_read_name
				tmp_dict = {}
				tmp_dict[non_uniq_read_name] = dict_read_name_genes_names[non_uniq_read_name]
				dict_read_name_genes_names.clear() # only one read stored
				dict_read_name_genes_names = tmp_dict	
		  flag_result = 0
		  (temp_read_name, temp_interval_r0, temp_interval_r1) = initalize_read_name_and_interval(r[0], r[1]) 
		  continue
            except KeyError:
               pass
            if ( r[0] and r[0].aQual < minaqual ) or ( r[1] and r[1].aQual < minaqual ):
               lowqual += 1
               write_to_samout( r, "too_low_aQual" )
               continue         
          
         try:
	    # --
            if overlap_mode == "union":
               fs = set()
               for iv in iv_seq: # interval from bam file for each fragment
                  if iv.chrom not in features.chrom_vectors:
                     raise UnknownChrom
                  for iv2, fs2 in features[ iv ].steps():
		     	if debug:
				print "****Unique_feature %s and feature_interval %s" %(fs2,iv2)	
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
	    # --
	    # fs is fs_exon
	    # -- @TODO: ADD A TAG TO SPECIFY WHICH GFF TYPE IS USED (e.g. GENE OR EXON)
	    # -- HERE WE COMMENT THE FOLLOWING LINE AS GENE GFF IS USED AND ADD LINE TO ASSIGN fs_exons to fs
	    # fs_genes = set_exons_to_set_genes(fs)
	    fs_genes = fs
            if fs_genes is None or len( fs_genes ) == 0:
               write_to_samout( r, "no_feature" )
               empty += 1
		# ambiguous read count and/or one of the read pair mapping on different gene (potential gene fusion events)...
		# elif len( fs ) > 1:
            elif len( fs_genes ) > 1:
	       """ for temp_exon in list(fs):
		
		       if (temp_exon == "470_ENSSSCG00000010836"):
				print "%s--%s" % (r[0].read.name, fs_genes)
	       """ 
	       if debug:	
			print "fs_genes:%s readname:%s" %(fs_genes,r[0].read.name) # added to check ambiguous reads
	       write_to_samout( r, "ambiguous[" + '+'.join( fs ) + "]" )
               ambiguous += 1
            else:
	       if debug:
		       #print "DEBUG::CR:: len(fs) <-> 1:: fs = %s" %fs
			pass
               write_to_samout( r, list(fs)[0] )
		
	       # -- See method
	       for temp_exon in list(fs):
		
		       if (temp_exon == "470_ENSSSCG00000010836"):
				pass		#print "%s" % (r[0].read.name)
	       # -- @NEW_UPDATE: REMOVE CALL TO METHOD count_and_rescue_unique_reads AS THIS WAS USED TO COUNT UNIQUE READS PER EXON
	       # -- ADDED CR: check if in list of 14 non-unique reads
               rr2 = r[0] if r[0] is not None else r[1]
	       if debug:

		       #print "ZZZZ_DEBUG::CR:: before call to count_and_rescue_unique_reads: readName = %s\ttag_nonunique_NH= %s" % (rr2.read.name, tag_nonunique_NH)
			pass
	       if not tag_nonunique_NH:
		       (dict_gene_unique_counts, counts, dict_gene_exons_unique_counts) = count_and_rescue_unique_reads(tag_gff, counts, dict_gene_unique_counts, fs, dict_gene_exons_unique_counts)
		       if debug:
			       # -- ADDED CR: check if read in list of 14 non-unique reads
			       # rr2 = r[0] if r[0] is not None else r[1]
			       if rr2.read.name in ["SRR001362.16709855","SRR001362.18249582","SRR001362.26237172","SRR001362.4561888","SRR001362.15281800",\
						    "SRR001362.9224444","SRR001362.3827859","SRR001362.18294547","SRR001362.17122208","SRR001362.23633791",\
						    "SRR001362.3657620","SRR001362.7580262","SRR001362.15547456","SRR001362.27524161"]:
				       print "DEBUG::CR::NU_reads:: fs = %s" % fs
				       print "DEBUG::CR::NU_reads:: dict_gene_unique_counts = %s" % (dict_gene_unique_counts)
	       # -- END ADDED CR
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

      # -- At same level as the for loop (outside of the for loop) - column: 7
      #fh_read_names_gene_names.close()
      if dict_read_name_genes_names.keys() != []:
	#print "dict_read_name_genes_names passing"
	non_uniq_read_name = dict_read_name_genes_names.keys()[0]
	fs_genes_names = dict_read_name_genes_names[non_uniq_read_name]
	fh_read_names_gene_names.write("%s\t%s\n" % (non_uniq_read_name, "\t".join(list(fs_genes_names)) ))
	#dict_read_name_genes_names.clear()
      # -- 

      fh_read_names_gene_names.close() 
      # fh_read_names_gene_names = open("non_unique_read_names_gene_names.tsv", 'r')	
      # fh_read_names_gene_names_4reading = open(outfile_non_unique_read_names_gene_names, 'r')	
      fh_read_names_gene_names_4reading = open(filename_read_names_gene_names, 'r')	
      # -- Find proportions of non-unique reads per genes
      # -- NEW_UPDATE: COMMENT LINE AS PROPORTIONS ARE CALCULATED IN SCRIPT_ASSIGN_*.py
      # dict_gene_non_unique_proportions, dict_gene_non_unique_uniform_proportions = weighting_non_unique_reads_with_max_exon_count(fh_read_names_gene_names_4reading, dict_gene_exons_unique_counts)
      fh_read_names_gene_names_4reading.close() 
      
   except:
      if not pe_mode:
         sys.stderr.write( "Error occured in %s.\n" % read_seq.get_line_number_string() )
      else:
         sys.stderr.write( "Error occured in %s.\n" % read_seq_pe_file.get_line_number_string() )
      raise

   if not quiet:
      sys.stderr.write( "%d sam %s processed.\n" % ( i, "lines " if not pe_mode else "line pairs" ) )
         
   if samoutfile is not None:
      samoutfile.close()
  # print "dict_nonunique.keys() = ", dict_nonunique.keys()
   #print "dict_non_unique_exons.keys() = ", dict_non_unique_exons.keys()
   #print "Gene\tUnique_reads(exon)\tNon_unique_reads(exon)\tTotal\tUnique(Gene)\tNon_Unique_(Gene)\tNon-Unique_proportions\t Total(Non-Unique_proportions+Unique_count_per_gene)\t Non_Unique_Uniform_Proportions(Gene) \tTotal(Unique(Gene)+ Non_Unique_Uniform_Proportions(Gene)))\n"
   print "Gene\tUnique(Gene)\tNon_Unique_(Gene)\tTotal(Unique(Gene)+ Non_Unique(Gene)\n"
   if tag_gff == "exon_gff":
	   tuples_genenames_exontag = [(fn.split("_")[1], fn) for fn in counts.keys()]
   elif tag_gff == "gene_gff":
	   tuples_genenames_exontag = [(fn, fn) for fn in dict_gene_unique_counts.keys()]
   tuples_genenames_exontag.sort()
   #for fn in sorted( counts.keys() ):
   #	gene_name = fn.split("_")[1].strip()	

   # fh_save_gene_unique_reads_total_rescued = open("table_genenames_unique_counts_total_proportion_rescued.tsv", 'w')
   fh_save_gene_unique_reads_total_rescued = open(outfile_table_genenames_unique_counts_total_proportion_rescued, 'w')
   fh_save_gene_unique_reads_total_rescued.write("Gene\tUnique(Gene)\tTotal(Unique(Gene)+Non_Unique_proportions)\t Total(Unique(Gene)+ Non_Unique_Uniform_Proportions(Gene)) \n")
   previous_gene_name = "NA"

   for gene_name, fn in tuples_genenames_exontag:
	gene_name = gene_name.strip()
	fn = fn.strip()
#	if gene_name not in dict_nonunique.keys():
#		dict_nonunique[gene_name] = 0
#	if gene_name not in dict_gene_unique_counts.keys():
#		dict_gene_unique_counts[gene_name] = 0
#	if fn not in dict_non_unique_exons.keys():
#		dict_non_unique_exons[fn] = 0
	#print "%s\t%d\t%d\t%d\t%i\t%i\n" % ( fn, counts[fn] ,dict_non_unique_exons[fn], counts[fn] + dict_non_unique_exons[fn], dict_gene_unique_counts[gene_name], dict_nonunique[gene_name] )
	if tag_gff == "exon_gff":	
		if gene_name in dict_gene_non_unique_proportions.keys():
			print "%s\t%d\t%d\t%d\t%i\t%i\t%1.2f\t%1.2f\t%i\t%i\n" % ( fn, counts[fn] ,dict_non_unique_exons[fn], counts[fn] + dict_non_unique_exons[fn], \
										   dict_gene_unique_counts[gene_name], dict_nonunique[gene_name], dict_gene_non_unique_proportions[gene_name], \
										   dict_gene_unique_counts[gene_name] + dict_gene_non_unique_proportions[gene_name], \
										   dict_gene_non_unique_uniform_proportions[gene_name], \
										   dict_gene_unique_counts[gene_name] + dict_gene_non_unique_uniform_proportions[gene_name] )
			if gene_name != previous_gene_name and gene_name != "NA":
				line2store = "%s\t%i\t%i\t%i\n" % ( gene_name, dict_gene_unique_counts[gene_name], \
								    dict_gene_unique_counts[gene_name] + dict_gene_non_unique_proportions[gene_name], \
								    dict_gene_non_unique_uniform_proportions[gene_name] + dict_gene_unique_counts[gene_name] )
				fh_save_gene_unique_reads_total_rescued.write(line2store)
		
		else:
			# -- No non-unique reads for that gene_name
			print "%s\t%d\t%d\t%d\t%i\t%i\t%i\t%i\t%i\t%i\n" % ( fn, counts[fn] ,dict_non_unique_exons[fn], counts[fn] + dict_non_unique_exons[fn], \
									     dict_gene_unique_counts[gene_name], dict_nonunique[gene_name], 0 , dict_gene_unique_counts[gene_name], \
									     0, 0 + dict_gene_unique_counts[gene_name]) 
			if gene_name != previous_gene_name and gene_name != "NA":
				line2store = "%s\t%i\t%i\t%i\n" % ( gene_name, dict_gene_unique_counts[gene_name], dict_gene_unique_counts[gene_name], 0 + dict_gene_unique_counts[gene_name])
				fh_save_gene_unique_reads_total_rescued.write(line2store)
   	elif tag_gff == "gene_gff": ##@TODO
		# counts{} contains count for exons and is depricated ....
		if gene_name in dict_gene_unique_counts.keys():
			print "%s\t%i\t%i\t%i\n" % ( fn, counts[fn], dict_gene_unique_counts[gene_name], dict_nonunique[gene_name] )
		else:
			# -- No non-unique reads for that gene_name
			print "%s\t%i\t%i\t%i\n" % ( fn, counts[fn], dict_gene_unique_counts[gene_name], 0 )
		
	# -- Re-initialise gene name
	previous_gene_name = gene_name
	
	#	print "%s\t%d\t%d\t%d\n" % ( fn, counts[fn] , dict_nonunique[fn], counts[fn] + dict_nonunique[fn])
		
		##print "*%s\t%d" % (fn, dict_nonunique[fn])
   print "no_feature\t%d" % empty
   print "ambiguous\t%d" % ambiguous
   print "too_low_aQual\t%d" % lowqual
   print "not_aligned\t%d" % notaligned
   print "alignment_not_unique\t%d" % nonunique
   ## -- Close output file for non_unique read names --> [genes..]
   #fh_read_names_gene_names.close() 
   fh_save_gene_unique_reads_total_rescued.close()
      
def main():
   
   optParser = optparse.OptionParser( 
      
      usage = "%prog [options] sam_file gff_file",
      
      description=
         "This script takes an alignment file in SAM format and a " +
         "feature file in GFF format and calculates for each feature " +
         "the number of reads mapping to it. See " +
         "http://www-huber.embl.de/users/anders/HTSeq/doc/count.html for details.",
         
      epilog = 
         "Written by Simon Anders (sanders@fs.tum.de), European Molecular Biology " +
         "Laboratory (EMBL). (c) 2010. Released under the terms of the GNU General " +
         "Public License v3. Part of the 'HTSeq' framework, version %s." % HTSeq.__version__ )
         
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
# 
#   if len( args ) != 2:
   if len( args ) != 4:
      sys.stderr.write( sys.argv[0] + ": Error: Please provide (originally two) now 4 arguments.\n" )
      sys.stderr.write( "  Call with '-h' to get usage information.\n" )
      sys.exit( 1 )
      
   warnings.showwarning = my_showwarning
   try:
      # fh_read_names_gene_names = open("non_unique_read_names_gene_names.tsv", 'w')
      #fh_read_names_gene_names = open(args[3], 'w')
      filename_read_names_gene_names = args[3]
      #count_reads_in_features( args[0], args[1], args[2], opts.stranded, 
      #   opts.mode, opts.featuretype, opts.idattr, opts.quiet, opts.minaqual,
      #   opts.samout,fh_read_names_gene_names )
      count_reads_in_features( args[0], args[1], args[2], opts.stranded, 
         opts.mode, opts.featuretype, opts.idattr, opts.quiet, opts.minaqual,
         opts.samout, filename_read_names_gene_names)
      # fh_read_names_gene_names.close()

   except:
      sys.stderr.write( "Error: %s\n" % str( sys.exc_info()[1] ) )
      sys.stderr.write( "[Exception type: %s, raised in %s:%d]\n" % 
         ( sys.exc_info()[1].__class__.__name__, 
           os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
           traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
      sys.exit( 1 )

def my_showwarning( message, category, filename, lineno = None, line = None ):
   sys.stderr.write( "Warning: %s\n" % message )

if __name__ == "__main__":
   main()

