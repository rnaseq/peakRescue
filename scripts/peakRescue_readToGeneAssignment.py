#!/usr/bin/python

import os
import sys
import optparse
import warnings
import traceback

"""
PeakRescue workflow: Assignment of ambiguously mapped reads to associated genes.

"""


#####################################################################################################################
## -- Get the set of genes which are present in all contributions considered: 
## -- In peakRescue only one contribution is required i.e. the Peak coverage.
#####################################################################################################################

def get_gene_names_intersection(list_dict_gene_parameter_values):
	"""
	Get the intersection of all gene names in all dictionaries
	to make sure that all the genes in the final list have a peak parameter value.
	
	"""
	list_intersection_genenames_for_contribution = []
	tag_first_iteration = 1
	for dict_gene_parameter_values in list_dict_gene_parameter_values:
		if tag_first_iteration == 1:
			list_intersection_genenames_for_contribution = dict_gene_parameter_values.keys()
			s3 = set(list_intersection_genenames_for_contribution)
			tag_first_iteration += 1
		else:	
			tmp_list = dict_gene_parameter_values.keys()
			s2 = set(tmp_list)
			s3 = s2.intersection(s3)
	return list(s3)


#####################################################################################################################
## -- Probabilistic assignment of reads to genes
#####################################################################################################################

def _get_gene_parameter_values(tag_first_gene, list_parameter_values_for_given_read, gene_name, list_dict_gene_parameter_values):
	"""
	Get the sum of all peak values across all the genes an ambiguous read maps to. 
	This corresponds to the denominator in the probabilistic assignment calculation.
	This is calculated for every contribution type considered (e.g. peak). 
	In the current implementation of peakRescue only one contribution is required i.e. Peak coverage contribution.
	The index here corresponds to the contribution considered - e.g. index=0 for peak contribution.
	This generic function was originally written in order to allow for more than one contribution to be used if required. 
	
	Input: 
		+ tag_first_gene: binary (0/1) value to indicate whether this is the first gene encountered
		+ list_dict_gene_parameter_values: list containing dictionaries of parameter values 
		  As many dictionaries as there are contributions considered
		+ gene_name: a gene present in the read name-to- gene name mapping for a given ambiguously mapped read.
		+ list_dict_gene_parameter_values: a list containing as many dictionaries as there are contributions considered.
		  Each dictionary contains the peak value (dict values) for all genes (dict keys).
		  
	Output:
		+ list_parameter_values_for_given_read: 
		A list of sum of contributions across all genes a given ambiguously mapped read maps to. 
		  
	"""
	if tag_first_gene:
		for dict_gene_parameter_values in list_dict_gene_parameter_values:
				list_parameter_values_for_given_read.append(dict_gene_parameter_values[gene_name])
	else:
		index = 0
		for dict_gene_parameter_values in list_dict_gene_parameter_values:
			if gene_name in dict_gene_parameter_values.keys():
				list_parameter_values_for_given_read[index] += dict_gene_parameter_values[gene_name]
			index += 1
	return list_parameter_values_for_given_read


def _get_proportions(list_dict_proportions, list_parameter_values_for_given_read, gene_name, list_dict_gene_parameter_values):
	"""
	Calculates the proportion of ambiguously mapped reads to be assigned to a given gene. Proportions are summed across all 
	the ambiguously mapped reads that map to the given gene (gene_name) with a function call from within a loop over all 
	ambiguously mapped reads.	

	Input: 
		+ list_dict_proportions: list of dictionaries (one per contribution considered - here a single contribution: peak coverage) 
		  containing all the ambiguously mapped read proportions assigned to each gene
		+ list_parameter_values_for_given_read: a list (as many items as contributions considered) of sum of all the peak values across 
		  all the genes a given ambiguously mapped read maps to. 
  		+ gene_name: a gene present in the read name-to- gene name mapping for a given ambiguously mapped read.
		+ list_dict_gene_parameter_values: list containing dictionaries of parameter values. As many dictionaries as there are contributions considered
		  
	Output:
		+ list_dict_proportions: list of dictionaries (one per contribution considered - here a single contribution: peak coverage) 
		  containing all the ambiguously mapped read proportions assigned to each gene
		  		  
	"""
	for i in range(0, len(list_dict_gene_parameter_values)):
		dict_proportions = list_dict_proportions[i] 
		sum_parameter_values = list_parameter_values_for_given_read[i]
		dict_parameter_values = list_dict_gene_parameter_values[i]
		if sum_parameter_values != 0:
			proportion = float(dict_parameter_values[gene_name]) / sum_parameter_values
		else:
			proportion = 0
		if gene_name in dict_proportions.keys():
			dict_proportions[gene_name] += proportion
		else:
			dict_proportions[gene_name] = proportion
		list_dict_proportions[i] = dict_proportions

	return list_dict_proportions


def weighting_ambiguously_mapped_reads(fh_read_names_gene_names, list_dict_gene_parameter_values):
	"""
	Probabilistic assignment of ambiguous reads (e.g. multimapped or ambiguous uniquely mapped) 
	to each gene they map to based on each gene's peak coverage.
	Peak coverage is based on the uniquely mapped reads including the disambiguated uniquely mapped ones.
	
	Methods output  a list of dictionaries of all contributions considered ('list_dict_proportions') 
	over all ambiguously mapped reads for each gene (key=gene name)

	Input: 
		+ fh_read_names_gene_names: filehandle to the ambiguously mapped read name-to-gene name mappings
		+ list_dict_gene_parameter_values: list containing dictionaries of parameter values. 
		  As many dictionaries as there are contributions considered
		  
	Output:
		+ list_dict_proportions: list of dictionaries (one per contribution considered - here a single contribution: peak coverage) 
		  containing all the ambiguously mapped read proportions assigned to each gene

	"""
	# -- Get the intersection of all gene names in all dictionaries
	list_intersection_genenames_for_contribution =  get_gene_names_intersection(list_dict_gene_parameter_values)

	line = fh_read_names_gene_names.readline()
	if line.startswith("EnsemblID") or line.startswith("ensembl"):
		line = fh_read_names_gene_names.readline()
	
	dict_gene_non_unique_proportions = {}
	dict_gene_non_unique_uniform_proportions = {}
	
	# -- Initialise dictionaries to store proportions for each contribution
	list_dict_proportions = []
	for i in range(0, len(list_dict_gene_parameter_values)):
		dict_proportions = {}
		list_dict_proportions.append(dict_proportions)	
		
	while line:
		fields = line.split("\t")
		read_name = fields[0]
		list_gene_names = fields[1:]
		list_gene_names = [g.strip() for g in list_gene_names]
		sum_max_unique_counts_exons_all_genes = 0
		list_parameter_values_for_given_read = []
		## Calculate the denominator of all contributions by looping over all genes associated with an ambiguously mapped read
		tag_first_gene = 1
		for gene_name in list_gene_names:
			if gene_name in list_intersection_genenames_for_contribution:
				list_parameter_values_for_given_read = _get_gene_parameter_values(tag_first_gene, list_parameter_values_for_given_read, gene_name, list_dict_gene_parameter_values)
				tag_first_gene = 0 
		for gene_name in list_gene_names:
			if gene_name in list_intersection_genenames_for_contribution:
				list_dict_proportions = _get_proportions(list_dict_proportions, list_parameter_values_for_given_read, gene_name, list_dict_gene_parameter_values)
		# -- 
		line = fh_read_names_gene_names.readline()
	return list_dict_proportions


#####################################################################################################################
## -- Processing files 
## -- Extract dict of values from input data & 
## -- Create output file to store read proportions, FPKM values &
## -- Create the final combined output file with all read count from non-ambiguous unique and ambiguously mapped rescued
#####################################################################################################################

def get_dict_gene_name_value(fh_input, tag_type="not_htseq"):
	"""
	Create dictionary for a given parameter and store associated value.

	Input:
		+ fh_input contains a .tsv file with 2 columns: 
			1/ gene name
			2/ parameter value (e.g. peak value)
			3/ tag_type: specify if file is from htseq original - then we need to skip header if any and last few summary lines (ambiguous etc...)

	Output:
		+ Returns 'dict_gene_parameter_value': 
		  Keys: gene names and values: parameter value (cf. 2nd col).

	"""
	dict_gene_parameter_value = {}
	lines = fh_input.readlines()
	if tag_type == "htseq_original":
		lines = lines[:-5] 
	if tag_type == "proportions":
		lines = lines[1:] 
	for line in lines:
		fields = line.split("\t")
		fields = [f.strip() for f in fields]
		ensid=fields[0]; parameter_value = float(fields[1])
		dict_gene_parameter_value[ensid] = parameter_value
	return dict_gene_parameter_value


def save_ambiguous_proportions_per_gene_as_tsv(fh_out, list_dict_gene_proportions):
	"""
	Save ambiguously mapped proportions for each gene into a .tsv file.

	Input:
		+ fh_out: filehandle to store the tab separated (.tsv) output file.
		  Column 1: gene
		  Column 2: ambiguously mapped read proportions for the given gene
		+ list_dict_proportions: list of dictionaries (one per contribution considered - here a single contribution: peak coverage) 
		  containing all the ambiguously mapped read proportions assigned to each gene

	Returns: None

	"""
	# Prints all proportions for each gene - The order in the list_dict_gene_proportions defines the order in the output.
	# Proportion means contribution
	list_genes = list_dict_gene_proportions[0].keys()
	for gene_name in list_genes:
		string_contributions = ""
		sum_contribution = 0
		for i in range(0, len(list_dict_gene_proportions)):
			string_contributions += "\t%1.5f" % list_dict_gene_proportions[i][gene_name]
			sum_contribution += list_dict_gene_proportions[i][gene_name]
		if len(list_dict_gene_proportions) >1:
			# Generic function which include cases where more than one contribution would be considered for ambiguous read rescue
			normalised_sum = float(sum_contribution) / float(len(list_dict_gene_proportions))
			fh_out.write("%s%s\t%1.5f\n" % (gene_name, string_contributions, normalised_sum))
		else:
			# PeakRescue only relies on the peak coverage (single contribution) for ambiguous read rescue
			fh_out.write("%s%s\n" % (gene_name, string_contributions))


def get_fpkm_from_counts_per_gene(gene_id, count, tag_all_mappable_reads, dict_gene_length):
	"""
	Convert the total read count per gene into FPKM values.

	Input:
		+ gene_id: gene ID
		+ count: total read count for the given gene ID
		+ tag_all_mappable_reads: total library size
		+ dict_gene_length: global transcript gene length (gene: keys)

	Output:
		+ fpkm_value: FPKM value for the given gene ID.

	"""
    #tag_all_mappable_reads = sum([dict_gene_counts[g] for g in dict_gene_counts.keys()])
    #fpkm_value = float((10**9) * dict_gene_counts[gene_id]) / (tag_all_mappable_reads * dict_gene_length[gene_id])
	if gene_id in dict_gene_length.keys():
		fpkm_value = float((10**9) * count) / (tag_all_mappable_reads * dict_gene_length[gene_id])
	else:
		print "Error: no gene_id %s in dict_gene_length" % gene_id
		fpkm_value = "NA"
        return fpkm_value

#####################################################################################################################
## -- Generate the final output with read counts and FPKM for each gene
## -- @ToDo: Check if we still need this function - whether we combine the output files here or in the main pipeline
#####################################################################################################################
def store_combined_outputs_htseq_and_rescued_assigned_mm_proportions(fh_final_out, fh_htseq, fh_htseq_original, fh_mm_proportions, fh_amb_unique_proportions, dict_gene_unique_peak, dict_gt_genelength):
	"""
	Create final output file with all read count with non-ambiguous unique count as well as all the ambiguously mapped read proportions for each gene.
	Add the final FPKM value for each gene - i.e. normalised the total read count assigned to each gene using FPKM metric with gene length based on the global transcript of a gene.

	Input:
		+ fh_final_out: output file opened in 'w' mode
		+ fh_htseq: HTSeq -union- mode output: containing the non-ambiguous uniquely mapped read count per gene
		+ fh_htseq_original: HTSeq -union- mode output: containing the non-ambiguous uniquely mapped read count per gene
		  @comment: curren implementation requires twice filehandle on HTSeq as it is accessed and read twice once to extract the counts 
		  (get_dict_gene_name_value) and once to loop over lines. @todo: to be modified to reduce redundant input entries.
		+ fh_mm_proportions: fh (open in 'r') on the proportions of multimappers 
		+ fh_amb_unique_proportions: fh (open in 'r') on the proportions of ambigous uniquely mapped
		+ dict_gene_unique_peak: Dict containing peak values for each gene (key)
		+ dict_gene_length: global transcript gene length (gene: keys).

	Returns: None
	Writes data to file opened as fh_final_out.

	"""
	# -- Extract unique (HTSeq original) counts
	dict_gene_htseq_original = get_dict_gene_name_value(fh_htseq_original, "htseq_original")
	# -- Extract all ambiguous uniquely mapped reads' proportions per gene
	#dict_gene_amb_unique_prop =  _extract_ambiguously_mapped_reads_proportions_per_gene(fh_amb_unique_proportions)
	dict_gene_amb_unique_prop = get_dict_gene_name_value(fh_amb_unique_proportions, "proportions")
	# -- Extract all multimapped reads' proportions per gene
	dict_gene_mm_prop = get_dict_gene_name_value(fh_mm_proportions, "proportions")
	# -- Store gene name and HTseq values and ambiguously ampped reads (both read types) proportions (if gene not present in dict proportions - report 0).
	dict_gene_all_values = {}
	# -- Combine HTSeq output data with ambiguously mapped (2 types) reads' contributions 
	lines_htseq = fh_htseq.readlines()
	for line in lines_htseq[1:-6]: ## avoids header and footer lines
		#if line.startswith("ENS"):
		fields = line.split("\t"); fields = [f.strip() for f in fields]
		ensid=fields[0]; unique_disambiguated = int(fields[1]); ambiguous_c = float(fields[3]); mm_c = float(fields[2])
		# -- Add contributions: amb_unique and mm
		if ensid in dict_gene_amb_unique_prop.keys() and ensid in dict_gene_mm_prop.keys():
			dict_gene_all_values[ensid] = [unique_disambiguated, ambiguous_c, dict_gene_amb_unique_prop[ensid], mm_c, dict_gene_mm_prop[ensid]]
		# -- Add contributions: amb_unique only (report 0 for mm)
		elif ensid in dict_gene_amb_unique_prop.keys() and ensid not in dict_gene_mm_prop.keys():
			dict_gene_all_values[ensid] = [unique_disambiguated, ambiguous_c, dict_gene_amb_unique_prop[ensid], mm_c, 0]
		# -- Add contributions:  mm only (report 0 for amb_unique)
		elif ensid not in dict_gene_amb_unique_prop.keys() and ensid in dict_gene_mm_prop.keys():
			dict_gene_all_values[ensid] = [unique_disambiguated, ambiguous_c, 0, mm_c, dict_gene_mm_prop[ensid]]
		# -- Report 0 for both contributions if gene absent from both dicts.
		else:
			dict_gene_all_values[ensid] = [unique_disambiguated, ambiguous_c, 0, mm_c, 0]
	## -- get sum of counts for RPKM calculation
	tag_all_mappable_reads_peakrescue = sum( [ dict_gene_htseq_original[ensid] + dict_gene_all_values[ensid][0] + dict_gene_all_values[ensid][2] + dict_gene_all_values[ensid][4] \
																		for ensid in dict_gene_all_values.keys() ] ) 
	tag_all_mappable_reads_htseqoriginal = sum( [ dict_gene_htseq_original[ensid] for ensid in dict_gene_all_values.keys() ] )
	
	## -- Write combined output to file
	for ensid in dict_gene_all_values.keys():
		# -- Peak information
		if ensid in dict_gene_unique_peak.keys():
			tag_peak = dict_gene_unique_peak[ensid]
		else:
			tag_peak = 0
		htseq_orig_unique = dict_gene_htseq_original[ensid]
		unique_disambiguated = dict_gene_all_values[ensid][0] 
		ambiguous_c = dict_gene_all_values[ensid][1] 
		amb_unique_prop = dict_gene_all_values[ensid][2] 
		mm_c = dict_gene_all_values[ensid][3] 
		mm_prop = dict_gene_all_values[ensid][4] 
		final_count = htseq_orig_unique + unique_disambiguated + amb_unique_prop + mm_prop
		fpkm_peakrescue_final_counts = get_fpkm_from_counts_per_gene(ensid, final_count, tag_all_mappable_reads_peakrescue , dict_gt_genelength)
		fpkm_htseqoriginal = get_fpkm_from_counts_per_gene(ensid, htseq_orig_unique, tag_all_mappable_reads_htseqoriginal, dict_gt_genelength)
		if fpkm_htseqoriginal != "NA" and fpkm_peakrescue_final_counts != "NA":
			fh_final_out.write("%s\t%i\t%i\t%i\t%1.2f\t%i\t%1.2f\t%1.2f\t%1.2f\t%i\t%1.2f\n" % (ensid, htseq_orig_unique, unique_disambiguated, ambiguous_c, \
														amb_unique_prop, mm_c, mm_prop, \
														final_count, fpkm_peakrescue_final_counts, tag_peak, fpkm_htseqoriginal))
		elif  fpkm_htseqoriginal == "NA" and fpkm_peakrescue_final_counts == "NA":
			fh_final_out.write("%s\t%i\t%i\t%i\t%1.2f\t%i\t%1.2f\t%1.2f\t%s\t%i\t%s\n" % (ensid, htseq_orig_unique, unique_disambiguated, ambiguous_c, \
														amb_unique_prop, mm_c, mm_prop, \
														final_count, fpkm_peakrescue_final_counts, tag_peak, fpkm_htseqoriginal))
		elif fpkm_htseqoriginal == "NA":
			fh_final_out.write("%s\t%i\t%i\t%i\t%1.2f\t%i\t%1.2f\t%1.2f\t%1.2f\t%i\t%s\n" % (ensid, htseq_orig_unique, unique_disambiguated, ambiguous_c, \
														amb_unique_prop, mm_c, mm_prop, \
														final_count, fpkm_peakrescue_final_counts, tag_peak, fpkm_htseqoriginal))
		elif fpkm_peakrescue_final_counts == "NA":
			fh_final_out.write("%s\t%i\t%i\t%i\t%1.2f\t%i\t%1.2f\t%1.2f\t%s\t%i\t%1.2f\n" % (ensid, htseq_orig_unique, unique_disambiguated, ambiguous_c, \
														amb_unique_prop, mm_c, mm_prop, \
														final_count, fpkm_peakrescue_final_counts, tag_peak, fpkm_htseqoriginal))
	return dict_gene_all_values


#####################################################################################################################
#####################################################################################################################
## -- MAIN FUNCTION WHICH PERFORMS ALL THE STEPS REQUIRED FOR THE PROBABILISTIC ASSIGNMENT 
## -- CALLED FROM MAIN() 
#####################################################################################################################
#####################################################################################################################

def peakrescue_probabilistic_assignment(main_dir, peak_filename, mappings_reads2genes_filename, gene_length_filename, readtype):

	## os.chdir(main_dir) ## output dir
	list_tag_read_types = ["multimappers", "ambiguous_unique"]
	dict_readtype_proportions2gene_mapping = {}
	tag_all_genes_peak = "tag_all_genes_peak"

	## ----------------------------------------------------------------------------------------------------------------
	## -- MAPPING READ NAME -> GENE NAMES (ADDITIONAL OUTPUT FILE FROM MODIFIED HTSeq)
	## -- OUTPUT FILE: AMBIGUOUSLY MAPPED READS (EITHER: AMBIGUOUS UNIQUE OR MULTIMAPPED READS) PROPORTIONS ASSIGNED TO EACH GENE
	## ----------------------------------------------------------------------------------------------------------------
	passed_read_type = readtype
	contribution_output_filename = os.path.join(main_dir, "results_peakrescue_readtype_%s.tsv" % passed_read_type)
	
	if passed_read_type == "multimappers":
		mapping_multimapped_rn_gn = mappings_reads2genes_filename
		# -- Peak-based contribution calculated for mulitmapped reads
		#fh_read_names_gene_names = open(mapping_multimapped_rn_gn, 'r')
		fh_read_names_gene_names = open(mappings_reads2genes_filename, 'r')
		
		dict_readtype_proportions2gene_mapping[passed_read_type] = contribution_output_filename ## output_mappedreads_proportions
		fh_tmp_out = open(contribution_output_filename, 'w')
		fh_tmp_out.write("EnsemblID\tMultimapper proportions(based on UniquePeak contribution)\n")
		
		#contribution_output_filename = output_mappedreads_proportions
	elif passed_read_type == "ambiguous_unique":
		mapping_amb_unique_rn_gn = mappings_reads2genes_filename
		# -- Peak-based contribution calculated for the ambiguous uniquely mapped reads
		#fh_read_names_gene_names = open(mapping_amb_unique_rn_gn, 'r')
		fh_read_names_gene_names = open(mappings_reads2genes_filename, 'r')

		dict_readtype_proportions2gene_mapping[passed_read_type] =  contribution_output_filename ## output_ambiguous_unique_reads_proportions 
		fh_tmp_out = open(contribution_output_filename, 'w')
		fh_tmp_out.write("EnsemblID\tAmbiguous uniquely mapped reads' proportions(based on UniquePeak contribution)\n")
		#contribution_output_filename = output_ambiguous_unique_reads_proportions

	## ----------------------------------------------------------------------------------------------------------------
	## -- UNIQULEY MAPPED READS-BASED PEAK COVERAGE
	## ----------------------------------------------------------------------------------------------------------------
	fh_unique_peak = open(peak_filename, 'r')
	dict_gene_unique_peak = get_dict_gene_name_value(fh_unique_peak)
	list_dict_gene_parameter_values = [dict_gene_unique_peak]
	print "dict_gene_unique_peak = %s" % dict_gene_unique_peak

	## ----------------------------------------------------------------------------------------------------------------
	## -- GET AMBIGUOUSLY MAPPED READS' CONTRIBUTION PER GENE (USING PEAK COVERAGE) 
	## ----------------------------------------------------------------------------------------------------------------
	if passed_read_type in list_tag_read_types:
		list_dict_gene_proportions = weighting_ambiguously_mapped_reads(fh_read_names_gene_names, list_dict_gene_parameter_values)
		print "Number of genes with peak coverage contribution > 0: len(dict_proportions) = %i" % (len(list_dict_gene_proportions[0]))

		save_ambiguous_proportions_per_gene_as_tsv(fh_tmp_out, list_dict_gene_proportions)
		
		fh_read_names_gene_names.close()
		fh_unique_peak.close()
		fh_tmp_out.close()

	## ----------------------------------------------------------------------------------------------------------------
	## COMBINE HTSeq OUTPUT WITH BOTH: MULTIMAPPED READS AND AMBIGUOUS UNIQUELY MAPPED READS PROPORTIONS TO GENES
	## ----------------------------------------------------------------------------------------------------------------
	# -- Gene length based on global transcript
	#fh_genelength = open( os.path.join(main_dir, global_transcript_based_genelength) )
	fh_genelength = open( gene_length_filename )
	dict_gt_genelength  = get_dict_gene_name_value(fh_genelength)
	fh_genelength.close()

	# -- Account for all the genes in the contribution output
	if passed_read_type in list_tag_read_types:
		fh_input_contribution = open(contribution_output_filename, 'r')
		dict_genes_with_contribution = get_dict_gene_name_value(fh_input_contribution, "proportions") ## tag: "proportions" to remove header line.
		for g in dict_gt_genelength.keys():
			if g not in dict_genes_with_contribution.keys():
				dict_genes_with_contribution[g] = 0
		list_genes_sorted = [k for k in dict_genes_with_contribution.keys()]
		list_genes_sorted.sort()
		fh_out = open( "%s_all_genes.out" % os.path.splitext(contribution_output_filename)[0], 'w')
		for g in list_genes_sorted:
			fh_out.write( "%s\t%s\n" % (g, dict_genes_with_contribution[g]))
		fh_out.close()   
		fh_input_contribution.close()
		# -- 

	# -- Add genes with 0 peak in ensid-peak output
	if tag_all_genes_peak == "tag_all_genes_peak":
		fh_peak_new_output = open( "%s_all_genes.out" % os.path.splitext(peak_filename)[0], 'w')
		#list_genes_with_peak_not_0 = dict_gene_unique_peak.keys()
		for g in dict_gt_genelength.keys():
			if g not in dict_gene_unique_peak.keys():
				dict_gene_unique_peak[g] = 0
		list_genes_sorted = [k for k in dict_gene_unique_peak.keys()]
		list_genes_sorted.sort()
		for g in list_genes_sorted:
			fh_peak_new_output.write( "%s\t%s\n" % (g, dict_gene_unique_peak[g]))
		fh_peak_new_output.close()   

	# -- Deprecated - get final combined output in main pipeline script. 
	#fh_final_out = open(peak_rescue_workflow_output, 'w')
	#dict_gene_all_values = store_combined_outputs_htseq_and_rescued_assigned_mm_proportions(fh_final_out, fh_htseq, fh_htseq_original, fh_mm_proportions, fh_amb_unique_proportions, dict_gene_unique_peak, dict_gt_genelength)
	#fh_htseq.close()
	#fh_mm_proportions.close() 
	#fh_amb_unique_proportions.close()
	#fh_final_out.close()

#####################################################################################################################
## -- FUNCTION MAIN: 
## -- COLLECTS PASSED PARAMETERS & CALLS 'peakrescue_probabilistic_assignment()'
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
         "(choices: ambiguous_unique, multimapped; default: multimapped)" )
            
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
   		peakrescue_probabilistic_assignment(opts.outdir, opts.peak_filename, opts.mappings_reads2genes_filename, opts.gene_length_filename, opts.readtype)
   except:
      sys.stderr.write( "Error: %s\n" % str( sys.exc_info()[1] ) )
      sys.stderr.write( "[Exception type: %s, raised in %s:%d]\n" % 
         ( sys.exc_info()[1].__class__.__name__, 
           os.path.basename(traceback.extract_tb( sys.exc_info()[2] )[-1][0]), 
           traceback.extract_tb( sys.exc_info()[2] )[-1][1] ) )
      sys.exit( 1 )


if __name__ == '__main__':
	main()
	
