### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###
### some functions were copied from KIPEs: https://doi.org/10.3390/plants9091103 and MaMYB: https://doi.org/10.1371/journal.pone.0239275 ###

### WARNING: do not use underscores in the bait MYB IDs ###

__version__ = "v0.11"

__usage__ = """
					python3 MYB_annotator.py
					--baits <MYB_SEQ_FILE>
					--info <MYB_CLASSIFICATION_FILE>
					--out <OUTPUT_DIR>
					--subject <SUBJECT_FILE (peptide,transcript,genomic sequences)> | --subjectdir <SUBJECT_FOLDER_WITH_SEQ_FILES>
					
					optional:
					--mode <TREE_BUILDER>(fasttree|raxml)[fasttree]
					--refmybs <REF_MYB_FILE>
					--cpu <NUMBER_OF_THREADS>[4]
					--cdsinput <CHANGES_EXPECTED_INPUT_TO_CDS>
					
					--mafft <PATH_TO_MAFFT>[mafft]
					--blastp <PATH_TO_AND_INCLUDING_BINARY>[blastp]
					--makeblastdb <PATH_TO_AND_INCLUDING_BINARY>[makeblastdb]
					
					--fasttree <PATH_TO_FASTTREE>[fasttree]
					--raxml <PATH_TO_RAXML>[raxml]				
					
					--simcutp <BLASTP_SIMILARITY_CUTOFF>[0.8]
					--poscutp <BLASTP_HIT_NUMBER_PER_BAIT_CUTOFF>[100]
					--lencutp	<BLASTP_MIN_LENGTH_CUTOFF>[50]
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, glob, sys, re, subprocess, dendropy
from operator import itemgetter
try:
	import hashlib
except ImportError:
	pass

# --- end of imports --- #


def load_BLAST_results( blast_result_file, similarity_cutoff, possibility_cutoff, length_cutoff ):
	"""! @brief load BLAST results """
	
	valid_blast_hits = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if float( parts[2] ) > similarity_cutoff:	#similarity is sufficient
				if float( parts[3] ) > length_cutoff:	#substantial part of query is matched
					try:
						valid_blast_hits[ parts[1] ].append( { 'gene': parts[0], 'score': float( parts[-1] ) } )
					except KeyError:
						valid_blast_hits.update( { parts[1]: [ { 'gene': parts[0], 'score': float( parts[-1] ) } ] } )
			line = f.readline()
	
	# --- reduce BLAST hit number to given number of candidate possibilities ---- #
	final_valid_blast_hits = {}
	for key in list(valid_blast_hits.keys()):
		hits = sorted( valid_blast_hits[ key ], key=itemgetter( 'score' ) )[::-1]
		genes = []
		for hit in hits:
			if hit['gene'] not in genes:
				if len( genes ) < possibility_cutoff:
					genes.append( hit['gene'] )
		final_valid_blast_hits.update( { key: genes } )
	
	return final_valid_blast_hits


def load_alignment( aln_file, tmp_mapping ):
	"""! @brief load alignment and replace query IDs by real sequence names """
	
	sequences = {}
	
	with open( aln_file ) as f:
		header = f.readline()[1:].strip()
		try:
			header = tmp_mapping[ header ]
		except KeyError:
			pass
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					try:
						header = tmp_mapping[ header ]
					except KeyError:
						pass
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def alignment_trimming( aln_file, cln_aln_file, occupancy ):
	"""! @brief remove all alignment columns with insufficient occupancy """
	
	alignment = load_alignment( aln_file, {} )
	# --- if there is an alignment (expected case) 
	if len( list(alignment.keys()) ) > 0:
		# --- identify valid residues in aligned sequences (columns with sufficient occupancy) --- #
		valid_index = []
		for idx, aa in enumerate( list(alignment.values())[0] ):
			counter = 0
			for key in list(alignment.keys()):
				if alignment[ key ][ idx ] != "-":
					counter += 1
			if counter / float( len( list(alignment.keys()) ) ) > occupancy:
				valid_index.append( idx )
		
		# --- generate new sequences --- #
		with open( cln_aln_file, "w" ) as out:
			for key in list(alignment.keys()):
				seq = alignment[ key ]
				new_seq = []
				for idx in valid_index:
					new_seq.append( seq[ idx ] )
				out.write( ">" + key + '\n' + "".join( new_seq ) + '\n' )
	# --- just in case the alignment file is empyt (is this possible?) ---#
	else:
		with open( cln_aln_file, "w" ) as out:
			out.write( "" )


def split_into_ingroup_and_outgroup( tree_file, in_list, out_list, neighbour_cutoff, mean_factor_cutoff, min_neighbour_cutoff ):
	"""! @brief split subject sequences into intgroup and outgroup based on reference MYBs and MYB-like sequences """

	# --- preparation of data structure --- #
	groups_around_ref_gene = {}
	for gene in ( in_list+out_list ):
		groups_around_ref_gene.update( { gene: [] } )
	
	# --- find node objects of reference genes --- #
	tree = dendropy.Tree.get_from_path( tree_file, "newick" )
	pdm = dendropy.PhylogeneticDistanceMatrix.from_tree( tree )
	my_mean_nearest_taxon_distance = pdm.mean_nearest_taxon_distance()
	
	ref_node_objects = {}
	for node in tree.taxon_namespace:
		try:
			groups_around_ref_gene[ node.label ]
			ref_node_objects.update( { node.label: node } )
		except KeyError:
			pass
	
	ref_gene_nodes = []
	ref_gene_nodes_dict_to_check = {}
	for gene in ( in_list+out_list ):
		ref_gene_nodes.append( ref_node_objects[ gene ] )
		ref_gene_nodes_dict_to_check.update( { ref_node_objects[ gene ]: None } )
	
	results = {}
	for i, t1 in enumerate( tree.taxon_namespace ):
		try:
			ref_gene_nodes_dict_to_check[ t1 ]
		except KeyError:	#only run analysis for non-reference sequences
			path_distances = []
			patristic_distances = {}
			for t2 in ref_gene_nodes:	#calculate distance to all other sequences in tree
				path_distance = pdm.path_edge_count( t1, t2)
				patr_distance = pdm.patristic_distance( t1, t2 )
				path_distances.append( { 'key': t2.label, 'val': path_distance } )
				patristic_distances.update( { t2.label: patr_distance } )
			in_counter = 0
			out_counter = 0
			sorted_distances = sorted( path_distances, key=itemgetter('val') )
			for each in sorted_distances[ : min( [ len( path_distances ), neighbour_cutoff ] ) ]:
				patr = patristic_distances[ each['key'] ]
				if patr < mean_factor_cutoff*my_mean_nearest_taxon_distance:	#exclude outliers on extremely long branches
					if each['key'] in in_list:	#check if smalles path_distances are to in- or outgroup baits
						in_counter += 1
					else:
						out_counter += 1
			if in_counter+out_counter > min_neighbour_cutoff:
				results.update( { t1.label: { 'score': float( in_counter ) / ( in_counter + out_counter ), 'in': in_counter, 'out': out_counter } } )
			else:
				results.update( { t1.label: { 'score': 0.0, 'in': in_counter, 'out': out_counter } } )
			#score ranges from 0 (non-MYB) to 1 (MYB)
	return results


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def load_bait_MYB_anno( info_file ):
	"""! @brief load MYB IDs into two list (in and out) """
	
	in_list, out_list = [], []
	with open( info_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if parts[1] == "in":
				in_list.append( parts[0] )
			else:
				out_list.append( parts[0] )
			line = f.readline()
	return in_list, out_list


def translate( seqs ):
	"""! @brief translates the given nucleotide sequence into peptide and splits at each star (stop codon) """
	
	genetic_code = {	'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I',
								'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T',
								'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N',
								'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H',
								'ACG': 'T', 'CCG': 'P', 'AGT': 'S', 'CAG': 'Q', 'CAA': 'Q',
								'CCC': 'P', 'TAG': '*', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C',
								'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R',
								'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'TGA': '*', 'GGA': 'G',
								'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S',
								'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S',
								'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A',
								'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T',
								'TTG': 'L', 'CGT': 'R', 'TAA': '*', 'CGC': 'R'
							}
	
	final_peptide_seqs = {}
	for key in seqs.keys():
		seq = seqs[ key ].upper()
		peptide = []
		for i in range( int( len( seq ) / 3.0 ) ):
			codon = seq[i*3:i*3+3]
			try:
				peptide.append( genetic_code[ codon ] )
			except:
				peptide.append( "*" )
		final_peptide_seqs.update( { key: "".join( peptide ) } )
	return final_peptide_seqs


def clean_input_FASTA_file( raw_subject_file, subject_file, mapping_table, cds_input ):
	"""! @brief clean input FASTA file """
	
	forbidden_characters = [ ";", ":", "(", ")", "_" ]
	
	with open( mapping_table, "w" ) as out:
		out.write( "InitialID\tCleanID\n" )
		sequences = {}
		with open( raw_subject_file ) as f:
			header = f.readline()[1:].strip()
			out.write( header + "\t" )
			if " " in header:
				header = header.split(' ')[0]
			for each in forbidden_characters:
				header = header.replace( each, "-" )
			header = header.encode("ascii", "ignore").decode()	#removal of non-ASCII characters
			out.write( header + "\n" )
			seq = []
			line = f.readline()
			while line:
				if line[0] == '>':
						sequences.update( { header: "".join( seq ) } )
						header = line.strip()[1:]
						out.write( header + "\t" )
						if " " in header:
							header = header.split(' ')[0]
						for each in forbidden_characters:
							header = header.replace( each, "-" )
						header = header.encode("ascii", "ignore").decode()
						out.write( header + "\n" )
						seq = []
				else:
					seq.append( line.strip() )
				line = f.readline()
			sequences.update( { header: "".join( seq ) } )
	
	if cds_input:
		sequences = translate( sequences )
	
	with open( subject_file, "w" ) as out:
		for key in sequences.keys():
			out.write( '>' + key + "\n" + sequences[ key ] + "\n" )


def load_ref_mybs( ref_mybs_file ):
	"""! @brief load IDs from given file """
	
	refmybs = {}
	with open( ref_mybs_file, "r" ) as f:
		line = f.readline()
		while line:
			if "\t" in line:
				parts = line.strip().split()
			else:
				refmybs.update( { line.strip(): line.strip() } )
			line = f.readline()
	return refmybs


def myb_group_assignment( ref_mybs, tree_file, myb_candidates ):
	"""! @brief assign new MYBs to reference MYBs e.g. the A.thaliana MYBs """
	
	new2ref_mapping_table, new_per_ref_myb = {}, {}
	
	#new2ref_mapping_table = { candiate1: ref1, candidate2: ref1, candiate3: ref14, ... }
	#new_per_ref_myb = { ref1: [ candidate1, candidate2 ], ref2: [candidate15], ref3: [], ... }
	
	# --- preparation of data structure --- #
	my_ref_mybs = list( sorted( ref_mybs.keys() ) )
	for gene in my_ref_mybs:	#reference MYBs
		new_per_ref_myb.update( { gene: [] } )
	
	for gene in myb_candidates:	#candidate genes of new species
		new2ref_mapping_table.update( { gene: None } )
	
	# --- find node objects of reference genes --- #
	tree = dendropy.Tree.get_from_path( tree_file, "newick" )
	pdm = dendropy.PhylogeneticDistanceMatrix.from_tree( tree )
	ref_node_objects = {}
	new_node_objects = {}
	for node in tree.taxon_namespace:
		try:
			new_per_ref_myb[ node.label ]
			ref_node_objects.update( { node.label: node } )
		except KeyError:
			try:
				new2ref_mapping_table[ node.label ]
				new_node_objects.update( { node.label: node } )
			except KeyError:
				pass
	
	ref_gene_nodes = []
	ref_gene_nodes_dict_to_check = {}
	candidate_gene_nodes = []
	canidate_gene_nodes_dict_to_check = {}
	for gene in my_ref_mybs:
		ref_gene_nodes.append( ref_node_objects[ gene ] )
		ref_gene_nodes_dict_to_check.update( { ref_node_objects[ gene ]: None } )
	for gene in myb_candidates:
		candidate_gene_nodes.append( new_node_objects[ gene ] )
		canidate_gene_nodes_dict_to_check.update( { new_node_objects[ gene ]: None } )
	
	for i, t1 in enumerate( candidate_gene_nodes ):
		#print ( "ref matching - " + str( i+1 ) + "/" + str( len( candidate_gene_nodes ) ) )
		distances = []
		for t2 in ref_gene_nodes:	#calculate distance to all other sequences in tree
			distance = pdm.path_edge_count( t1, t2)
			distances.append( distance )
		ref_myb = my_ref_mybs[ distances.index( min( distances ) ) ]
		new2ref_mapping_table[ t1.label ] = ref_myb
		new_per_ref_myb[ ref_myb ].append( t1.label )

	return new2ref_mapping_table, new_per_ref_myb


def load_myb_classification_from_file( tmp_result_table ):
	"""! @brief load MYB classification from file """
	
	myb_classification = {}
	with open( tmp_result_table, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			myb_classification.update( { parts[0]: float( parts[1] ) } )
			line = f.readline()
	return myb_classification


def R2R3_comain_check( seqs ):
	"""! @brief screen sequences for R2R3-MYB domain """
	
	R2R3_domain = {}
	for key in sorted( seqs.keys() ):
		seq = seqs[ key ]
		try:
			match = re.findall( "\w{5}W\w{85,100}W\w{7}", seq )[0]
			R2R3_domain.update( { key: match } )
		except:
			R2R3_domain.update( { key: "" } )
	return R2R3_domain


def check_MYB_IDs_across_files( bait_seq_file, info_file, ref_mybs_file ):
	"""! @brief check MYB IDs across the different files """
	
	myb_status = True
	# --- check FASTA file for forbidden characters --- #
	forbidden_characters = [ ";", ":", "(", ")", "_" ]
	seqs = load_sequences( bait_seq_file )
	header_string = "".join( seqs.keys() )
	for each in forbidden_characters:
		if each in header_string:
			sys.stderr.write( "Forbidden character detected in MYB IDs (bait FASTA file): " + each + "(occurrences:" + str( header_string.count( each ) ) + ")\n" )
			sys.stderr.flush()
			myb_status = False
	
	# --- check structure of info file --- #
	info_myb_IDs = {}
	with open( info_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts ) < 2:
				sys.stderr.write( "Issue in MYB info file (number of columns not 2): " + line + "\n" )
				sys.stderr.flush()
			else:
				info_myb_IDs.update( { parts[0]: None } )
				if parts[1] not in [ "in", "out" ]:
					sys.stderr.write( "Issue in MYB info file (unexpected status; only 'in' and 'out' permitted): " + line + "\n" )
					sys.stderr.flush()
			line = f.readline()
	
	# --- compare IDs between bait and info file --- #
	missing_in_fasta = []
	for myb in info_myb_IDs.keys():
		try:
			seqs[ myb ]
		except KeyError:
			missing_in_fasta.append( myb )
	if len( missing_in_fasta ) > 0:
		sys.stderr.write( "Unmatched MYB IDs (missing in bait FASTA file) :" + ";".join( missing_in_fasta ) + "\n" )
		sys.stderr.flush()
		myb_status = False
	
	missing_in_info = []
	for myb in seqs.keys():
		try:
			info_myb_IDs[ myb ]
		except KeyError:
			missing_in_info.append( myb )
	if len( missing_in_info ) > 0:
		sys.stderr.write( "Unmatched MYB IDs (missing in info file) :" + ";".join( missing_in_info ) + "\n" )
		sys.stderr.flush()
		myb_status = False
	
	if len( ref_mybs_file ) > 0:	#only try to check if file actually exists
		missing_ref_mybs = []
		with open( ref_mybs_file, "r" ) as f:
			line = f.readline()
			while line:
				x = line.strip()
				try:
					seqs[ x ]
				except KeyError:
					missing_ref_mybs.append( x )
				line = f.readline()
		if len( missing_ref_mybs ) > 0:
			sys.stderr.write( "Reference MYB IDs missing in FASTA and/or info file) :" + ";".join( missing_ref_mybs ) + "\n" )
			sys.stderr.flush()
			myb_status = False
		
	return myb_status


def md5_calculator( input_file ):
	"""! @brief calculate md5sum of given file """
	
	with open( input_file, "rb" ) as f:
		content = f.read()
	try:
		return hashlib.md5( content ).hexdigest()
	except NameError:
		return "n/a"


def generate_documentation_file( 	doc_file, bait_seq_file, info_file, output_folder, raw_subject_file,
															mode, blastp, makeblastdb, mafft, cpu, raxml, fasttree, ref_mybs_file,
															similarity_cutoff_p, possibility_cutoff_p, length_cutoff_p, cds_input
														):
	"""! @brief write documentation file with specified inputs and parameters """
	
	with open( doc_file, "w" ) as out:
		out.write( "Please cite XXX when using MYB_annotator.py.\n\n" )
		out.write( "MYB_annotator.py version: " + __version__ + "\n" )
		bait_seq_file_md5 = md5_calculator( bait_seq_file )
		out.write( "MYB bait file: " + bait_seq_file + "\t" + bait_seq_file_md5 + "\n" )
		info_file_md5 = md5_calculator( info_file )
		out.write( "MYB info file: " + info_file + "\t" + info_file_md5 + "\n" )
		raw_subject_file_md5 = md5_calculator( raw_subject_file )
		out.write( "Subject FASTA file: " + raw_subject_file + "\t" + raw_subject_file_md5 + "\n" )
		out.write( "Output folder: " + output_folder + "\n" )
		
		#--- optional --- #
		out.write( "Tool for tree construction: " + mode + "\n" )
		out.write( "CPUs: " + str( cpu ) + "\n" )
		if len( ref_mybs_file ) > 0:
			ref_myb_file_md5 = md5_calculator( ref_mybs_file )
			out.write( "Reference MYB file: " + ref_mybs_file + "\t" + ref_myb_file_md5 + "\n" )
		else:
			out.write( "Reference MYB file: n/a\n" )
		if cds_input:
			out.write( "Type of input: CDS\n" )
		else:
			out.write( "Type of input: PEP\n" )
		
		# --- paths to tools --- #
		out.write( "blastp path: " + blastp + "\n" )
		out.write( "makeblastdb path: " + makeblastdb + "\n" )
		out.write( "mafft path: " + mafft + "\n" )
		out.write( "raxml path: " + raxml + "\n" )
		out.write( "fasttree path: " + fasttree + "\n" )
		
		# ---- BLAST filter criteria --- #
		out.write( "Minimal BLASTp hit similarity cutoff: " + str( similarity_cutoff_p ) + "\n" )
		out.write( "Maximal number of BLASTp hits per bait: " + str( possibility_cutoff_p ) + "\n" )
		out.write( "Minimal BLASTp hit alignment length: " + str( length_cutoff_p ) + "\n" )
		
		# --- add tool versions --- #
		try:
			mafft_version_raw = subprocess.Popen( args=mafft + " --version", stderr=subprocess.PIPE, shell=True )
			mafft_version = mafft_version_raw.stderr.read()
			out.write ( "MAFFT version: " + str( mafft_version )[2:-3] + "\n" )	#remove characters introduced through binary
		except:
			out.write ( "MAFFT version detection failed.\n" )	#if no MAFFT installation was detected
		out.write ( "FastTree version: PLEASE_ADD_MANUALLY\n"  )	#version not available via command
		try:
			raxml_version_raw = subprocess.Popen( args=raxml + " --version", stdout=subprocess.PIPE, shell=True )
			raxml_version = str( raxml_version_raw.stdout.read() ).strip()
			out.write ( "RAxML version: " + ( raxml_version[4:65]) + "...\n" )	#remove characters introduced through binary
		except:
			out.write ( "RAxML version detection failed.\n" )	#if no RAxML installation was detected


def main( arguments ):
	"""! @brief run everything """
	
	bait_seq_file = arguments[ arguments.index('--baits')+1 ]
	info_file = arguments[ arguments.index('--info')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	if '--subject' in arguments:
		raw_subject_files = [ arguments[ arguments.index('--subject')+1 ] ]
	else:
		subject_file_dir = arguments[ arguments.index('--subjectdir')+1 ]
		if not subject_file_dir[-1] == "/":
			subject_file_dir + "/"
		extensions = [ ".fasta", ".fa", ".fas", ".FASTA", ".FA", ".FAS" ]
		raw_subject_files = [ ]
		for each in extensions:
			raw_subject_files += glob.glob( subject_file_dir + "*" + each )
		raw_subject_files = list( sorted( raw_subject_files ) )
	
	if '--mode' in arguments:
		mode = arguments[ arguments.index('--mode')+1 ]
		if mode not in [ "fasttree", "raxml" ]:
			mode = "fasttree"
	else:
		mode = "fasttree"
	
	if '--blastp' in arguments:
		blastp = arguments[ arguments.index('--blastp')+1 ]
	else:
		blastp = "blastp"
	if '--makeblastdb' in arguments:
		makeblastdb = arguments[ arguments.index('--makeblastdb')+1 ]
	else:
		makeblastdb = "makeblastdb"
	if '--mafft' in arguments:
		mafft = arguments[ arguments.index('--mafft')+1 ]
	else:
		mafft = "mafft"
	if '--cpu' in arguments:
		cpu = int( arguments[ arguments.index('--cpu')+1 ] )
	else:
		cpu = 4
	if '--raxml' in arguments:
		raxml = arguments[ arguments.index('--raxml')+1 ]
	else:
		raxml = "raxml"
	if "--fasttree" in arguments:
		fasttree = arguments[ arguments.index('--fasttree')+1 ]
	else:
		fasttree = "fasttree"
	
	if '--refmybs' in arguments:
		ref_mybs_file = arguments[ arguments.index('--refmybs')+1 ]
	else:
		ref_mybs_file = ""
	
	# --- BLAST hit cutoffs --- #
	if "--simcutp" in arguments:
		similarity_cutoff_p = float( arguments[ arguments.index('--simcutp')+1 ] )
	else:
		similarity_cutoff_p=0.8
	if '--poscutp' in arguments:
		possibility_cutoff_p = int( arguments[ arguments.index('--poscutp')+1 ] )
	else:
		possibility_cutoff_p=100
	if '--lencutp' in arguments:
		length_cutoff_p = int( arguments[ arguments.index('--lencutp')+1 ] )
	else:
		length_cutoff_p=75
	
	if '--cdsinput' in arguments:
		cds_input = True
	else:
		cds_input = False
	
	neighbour_cutoff=10	#numbers of closest neightbour that is considered in ingroup/outgroup classification
	mean_factor_cutoff=3	#X*average nearest neighbor distance
	min_neighbour_cutoff = 0	#minimal number of valid bait sequences (ingroup+outgroup) in range - 1 
	
	if output_folder[-1] != "/":
		output_folder += "/"
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	for jidx, raw_subject_file in enumerate( raw_subject_files ):	#use jidx to generate unique IDs for all jobs
		# --- prepare output folder for each job if there are multiple --- #
		if len( raw_subject_files ) == 1:
			job_output_folder = output_folder
		else:
			job_ID = raw_subject_file.split('/')[-1].split('.')[0]
			job_output_folder = output_folder + str( jidx ).zfill(5) + "_" + job_ID + "/"
		
		if not os.path.exists( job_output_folder ):
			os.makedirs( job_output_folder )
	
		# --- validation of inputs --- #
		#check if all MYB baits are listed in the info file
		subject_file = job_output_folder + "clean_subject_sequences.fasta"
		mapping_table = job_output_folder + "raw_subject_to_clean_subject_mapping_table.txt"
		if not os.path.isfile( subject_file ):
			clean_input_FASTA_file( raw_subject_file, subject_file, mapping_table, cds_input )	#remove illegal characters from subject sequence headers
		MYB_check_status = check_MYB_IDs_across_files( bait_seq_file, info_file, ref_mybs_file )
		if not MYB_check_status:
			sys.exit( "ERROR: analysis is stopped due to inconstistency of MYB IDs between files" )
		
		result_folder = job_output_folder + "RESULTS/"
		if not os.path.exists( result_folder ):
			os.makedirs( result_folder )
		
		doc_file = result_folder + "0_documentation.txt"
		generate_documentation_file( 	doc_file, bait_seq_file, info_file, job_output_folder, raw_subject_file,
															mode, blastp, makeblastdb, mafft, cpu, raxml, fasttree, ref_mybs_file,
															similarity_cutoff_p, possibility_cutoff_p, length_cutoff_p, cds_input
														)
			
		
		# --- find initial candidates --- #
		blast_result_file = job_output_folder + "blast_results.fasta"
		if not os.path.isfile( blast_result_file ):
			blast_db = job_output_folder + "blastdb"
			p = subprocess.Popen( args= makeblastdb + " -in " + subject_file + " -out " + blast_db + " -dbtype prot", shell=True )
			p.communicate()
			
			p = subprocess.Popen( args= "blastp -query " + bait_seq_file + " -db " + blast_db + " -out " + blast_result_file + " -outfmt 6 -evalue 0.001 -num_threads " + str( cpu ), shell=True )
			p.communicate()
		
		blast_results = load_BLAST_results( blast_result_file, similarity_cutoff_p, possibility_cutoff_p, length_cutoff_p )	#load valid BLASTp results
		
		subject_sequences = load_sequences( subject_file )
		candidate_file = result_folder + "1_initial_candidates.fasta"
		with open( candidate_file, "w" ) as out:
			for each in blast_results.keys():
				out.write( '>' + each + "\n" + subject_sequences[ each ] + "\n" )
		
		# --- construct phylogenetic tree --- #
		aln_input_file = job_output_folder + "alignment_input.fasta"
		if not os.path.isfile( aln_input_file ):
			p = subprocess.Popen( args= "cat " + bait_seq_file + " " + candidate_file + " > " + aln_input_file, shell=True )
			p.communicate()
		
		aln_file = job_output_folder + "alignment_input.fasta.aln"
		if not os.path.isfile( aln_file ):
			p = subprocess.Popen( args= mafft + " --quiet " + aln_input_file + " > " + aln_file, shell=True )
			p.communicate()
		
		cln_aln_file = job_output_folder + "alignment_input.fasta.aln.cln"
		if not os.path.isfile( cln_aln_file ):
			alignment_trimming( aln_file, cln_aln_file, occupancy=0.1 )
		
		if mode == "raxml":	#RAxML
			prefix = job_output_folder + "RAxML_tree"
			tree_file = prefix + ".raxml.bestTree"
			if not os.path.isfile( tree_file ):
				p = subprocess.Popen( args= " ".join( [ raxml, "--all --threads " + str( cpu ) + " --model LG+G8+F --msa", cln_aln_file, "--prefix", prefix ] ), shell=True )
				p.communicate()
		else:	#FastTree2
			tree_file = job_output_folder + "FastTree_tree.tre"
			if not os.path.isfile( tree_file ):
				p = subprocess.Popen( args= " ".join( [ fasttree, "-wag  -nopr -nosupport <", cln_aln_file, ">", tree_file ] ), shell=True )
				p.communicate()
		
		# --- analyze tree file --- #
		clean_mybs_file = result_folder + "2a_clean_MYBs.fasta"
		tmp_result_table = result_folder + "2b_in_out_MYB_analysis_results.txt"
		if not os.path.isfile( tmp_result_table ):
			in_list, out_list = load_bait_MYB_anno( info_file )
			sys.stdout.write( "Number of ingroup MYB baits: " + str( len( in_list ) ) + "\n" )
			sys.stdout.write( "Number of outgroup MYB baits: " + str( len( out_list ) ) + "\n" )
			sys.stdout.flush()
			myb_classification = split_into_ingroup_and_outgroup( tree_file, in_list, out_list, neighbour_cutoff, mean_factor_cutoff, min_neighbour_cutoff )
			#dictionary with subject IDs: values are in the range of 0 (no MYB) to 1 (MYB)
			
			with open( clean_mybs_file, "w" ) as out:
				with open( tmp_result_table, "w" ) as out2:
					out2.write( "ID\tScore\n" )
					candidate_order = list( sorted( myb_classification.keys() ) )
					for candidate in candidate_order:
						if myb_classification[ candidate ]['score'] > 0.5:
							out.write( '>' + candidate + "\n" + subject_sequences[ candidate ] + "\n" )
						out2.write( "\t".join( list( map( str, [ 	candidate,
																						myb_classification[ candidate ]['score'],
																						myb_classification[ candidate ]['in'],
																						myb_classification[ candidate ]['out']
																					] ) ) ) + "\n" )
		else:
			myb_classification = load_myb_classification_from_file( tmp_result_table )
		clean_mybs = load_sequences( clean_mybs_file )
		
		# --- find closest reference MYB --- #
		if len( ref_mybs_file ) > 0:	#only performed if reference MYB file is provided
			group_around_ref_myb_file = result_folder + "3a_group_around_ref_MYBs.txt"	#produce table sorted by reference MYBs
			new_2_ref_myb_mapping_file = result_folder + "3b_new_2_ref_myb_mapping_file.txt"	#produce table sorted by subject sequences
			if not os.path.isfile( new_2_ref_myb_mapping_file ):
				ref_mybs = load_ref_mybs( ref_mybs_file )	#load IDs and trivial name from additional text file (dictionary)
				new2ref_mapping_table, new_per_ref_myb = myb_group_assignment( ref_mybs, tree_file, clean_mybs.keys() )
		
				with open( group_around_ref_myb_file, "w" ) as out:
					out.write( "RefMYB\tNewMYBs\n" )
					gene_order = list( sorted( new_per_ref_myb.keys() ) )
					for gene in gene_order:
						out.write( gene + "\t" + ";".join( new_per_ref_myb[ gene ] ) + "\n" )		
				
				with open( new_2_ref_myb_mapping_file, "w" ) as out:
					out.write( "NewMYB\tRefMYBs\n" )
					gene_order = list( sorted( new2ref_mapping_table.keys() ) )
					for gene in gene_order:
						out.write( gene + "\t" + new2ref_mapping_table[ gene ] + "\n" )
		
		# --- check for presence of R2R3 MYB domain --- #
		myb_domain_check_file = result_folder + "4_R2R3_MYB_domain_check.txt"	#produce table with R2R3-MYB domain status and sequence
		if not os.path.isfile( myb_domain_check_file ):
			clean_candidate_myb_sequences = load_sequences( clean_mybs_file )
			r2r3_domains = R2R3_comain_check( clean_candidate_myb_sequences )	#based on banana MYB paper: https://doi.org/10.1371/journal.pone.0239275
			proper_R2R3_MYB_counter = 0
			with open( myb_domain_check_file, "w" ) as out:
				out.write( "GeneID\tR2R3-MYB domain status\tR2R3-MYB domain\n" )
				candidates = list( sorted( clean_candidate_myb_sequences.keys() ) )
				for candidate in candidates:
					if len( r2r3_domains[ candidate ] ) > 0:
						out.write( "\t".join( [ candidate, "1", r2r3_domains[ candidate ] ] ) + "\n" )
						proper_R2R3_MYB_counter += 1
					else:
						out.write( "\t".join( [ candidate, "0", "" ] ) + "\n" )
			sys.stdout.write( "Number of proper R2R3-MYBs: " + str( proper_R2R3_MYB_counter ) + "\n" )
			sys.stdout.flush()
		
		# --- construct a final tree --- #
		
		


if '--baits' in sys.argv and '--info' in sys.argv and '--out' in sys.argv and '--subject' in sys.argv:
	main( sys.argv )
elif '--baits' in sys.argv and '--info' in sys.argv and '--out' in sys.argv and '--subjectdir' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
