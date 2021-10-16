# Automatic annotation of MYBs

### Please get in touch if you need help running the MYB annotator on your own dataset: [Boas Pucker (email)](mailto:b.pucker@tu-braunschweig.de?subject=[GitHub]MYB_annotator) ###

## Background ##
This tool allows an automatic identification and analysis of MYBs in a CDS or peptide sequence collection. Since the quality of structural annotation has substantially improved during the last years, the focus can now shift towards functional annotation. RNA-seq and long read sequencing of transcripts were fundamental for this development. There are numerous publications about R2R3-MYBs in various plant species. This tools allows the automatic investigation of the MYB gene family in any new species. The provided MYB sequences should be sufficient as baits in most plant species thus no manual preparation of input files is required. Results are a FASTA file containing the bona fide MYB sequences and a summary table with additional details. MYB candidates are screened for the presence of the conserved repeats and also for many other motifs. The function of all MYB candidates is predicted based on orthology to a previously characterized sequence (often in _Arabidopsis thaliana_).

## Installation ##

There is no installation of this tool required. Downloading and executing the script on a Linux system is sufficient. There is currently no support for other operating systems. Most required modules should be included in the initial Python installation. Requirements are (Linux installation instructions):

[Python3](https://www.python.org/) (sudo apt-get install python3.8)

[dendropy](https://dendropy.readthedocs.io/en/main/) (python3 -m pip install git+https://github.com/jeetsukumaran/DendroPy.git)



## Usage ##

```
Usage
python3 MYB_annotator.py --baits <FILE> --info <FILE> --out <DIR> [--subject <FILE> | --subjectdir <DIR>]

Mandatory:
--baits         STR     MYB bait sequence file
--info          STR     MYB info file
--out           STR     Output folder
--subject       STR     Subject file (OR: --subjectdir    STR     Subject dir)
					
Optional:
--mode          STR     Tool for tree construction (fasttree|raxml)[fasttree]
--refmybs       STR     Reference MYB file
--cpu           INT     Number of threads [4]
--cdsinput      NONE    Changes expected input to CDS
--keepnames     NONE    Prevents splitting of sequence names at first space

--ath           STR     FASTA file with A. thaliana MYBs
--name          STR     Prefix of output file names
--collapse      NONE    Reduces paralogs to one representative
--motifs        STR     File with motifs to check in candidate sequences
--cpub          INT     Number of threads for BLASTp[cpu]
--cpur          INT     Number of threads for RAxML[cpu]
					
--mafft         STR     Path to MAFFT [mafft]
--blastp        STR     Path to blastp [blastp]
--makeblastdb   STR     Path to makeblastdb [makeblastdb]			
--fasttree      STR     Path to FastTree [fasttree]
--raxml         STR     Path to RAxML [raxml]				
					
--simcutp       FLOAT   BLASTp similarity cutoff [0.8]
--poscutp       INT     Max number of BLASTp hits per bait [100]
--lencutp	INT     Min BLASTp alignment length [50]
```

`--baits` specifies a FASTA file that contains the MYB bait sequences. Some of these sequences are _bona fide_ MYBs and other are MYB-like sequences that are used as outgroup.

`--info` specifies a text file that contains a table with all MYB IDs (matching the baits FASTA file) in the first column. Columns are TAB-separted. The second column contains an assignment to the ingroup ('in') or outgroup ('out'). Ingroup and outgroup baits are used to classify MYB candidates as _bona fide_ MYBs or MYB-likes. For example, MYB-like/CDC5 sequences of _A. thaliana_ are marked as outgroups in the baits_MYBs_info.txt. They are considered the closest known relatives of MYB sequences.

`--out` specifies the output folder. This folder will be created if it does not exist already. All files generated during the analysis will be placed in this folder. There is a subfolder called 'RESULTS' that contains the important final files.

`--subject` specifies the FASTA file that should be analysed. The default expectation is that this file contains peptide sequences. However, it is also possible to supply a collection of coding sequences (CDS) if the `--cdsinput` flag is set.

`--subjectdir` specifies a folder with FASTA files. All files need to contain sequences of the same type (PEP or CDS). Supported file extensions are '.fasta', '.fas', '.fa', '.FASTA', '.FAS', '.FA', '.fna', '.FNA', '.cds', '.CDS', '.pep', and '.PEP'.

`--mode` specifies the tool for the construction of phylogenetic trees. [RAxML](https://doi.org/10.1093/bioinformatics/btz305) ('raxml') or [FastTree2](https://doi.org/10.1371/journal.pone.0009490) ('fasttree') can be specified. Default is 'fasttree'. Running the tree construction with RAxML is expected to result in slightly more accurate results, but would increase the run time from minutes to days in many cases.

`--refmybs` specifies a text file with one MYB ID per line. All IDs listed in this file need to be present in the baits FASTA file and the info file. These IDs can be used for the functional annotation of the newly identified MYBs. This assignment is performed in two ways leading to two different output tables. The functional annotation is lifted from previously characterized sequences to newly identified orthologs. 

`--cpu` specifies threads to use for the BLASTp search and also for the RAxML tree construction. Default: 4.

`--cpub` specifies specifically the number of threads to use for the BLASTp search. This option can be used if the number of cores to use for BLASTp and RAxML is different. Default is the value of `--cpu`.

`--cpur` specifies specifically the number of threads to use for the RAxML analysis. This option is only relevant if RAxML is used for the tree construction. This option can be used if the number of cores to use for BLASTp and RAxML is different. Default is the value of `--cpu`.

`--cdsinput` changes the input expectation from PEP to CDS. The file needs to be supplied through the `--subject`. This flag can be included at any place in the command and does not need a value.

`--keepnames` this flag deactivates the trimming of the names of input sequences at the first space. This flag can be included at any place in the command and does not need a value.

`--ath` specifies a FASTA file with Arabidopsis thaliana MYB sequences that are used to construct a final tree with the candidate sequences and these reference MYB sequences.

`--name` specifies a string that will be attached as prefix to all file names in the result folder. This option is useful if the result files of different analyses will be combined in one folder in a later step of the project. This name should NOT use any special characters. Just using the numbers 0-9 and the normal characters a-z or A-Z should be save.

`--collapse` this flag triggers and additional analysis step. Group of paralogs will be represented by a single (the longest) sequence. This is helpful if a transcriptome assembly is analyzed which might have many isoforms.

`--motifs` specifies a file with motif sequences and names. All candidate sequences will be screened for these motifs. Many subgroups are characterized by conserved motifs outside the R2R3 domain. This option allows an automatic analysis of all candidates and summarizes the results in a table.

`--mafft` specifies the full path to [MAFFT](https://doi.org/10.1093/molbev/mst010) including the binary.

`--blastp` specifies the full path to [BLASTp](https://dx.doi.org/10.1093%2Fnar%2F25.17.3389) including the binary.

`--fasttree` specifies the full path to [FastTree2](https://doi.org/10.1371/journal.pone.0009490). This is required if fasttree is not included in the PATH variable. If in doubt, please specify the full path.

`--raxml` specifies the full path to [RAxML](https://doi.org/10.1093/bioinformatics/btz305).

`--simcutp` specifies the minimal similarity of BLASTp hits to be considered in the initial identification of MYB candidates. Default: 0.8 (80%).

`--poscutp` specifies the maximal number of hits (possibilities) per bait sequence in the BLASTp analysis. Default: 100.

`--lencutp` specifies the minimal alignment length of BLASTp hits to be considered in the initial identification of MYB candidates. Default: 50.

WARNING: If errors occur or if parameters are changed, it is necessary to delete all output files to get clean files based on the new parameters. Existing files are usually not overwritten.


## Example ##

This example should help to understand how the script is used. The required MYB baits FASTA file, the MYB info, and the optional reference MYB file (_Arabidopsis thaliana_ MYBs) are included in this repository.

```
python3 /my_folder/MYB_annotator.py \
--baits /my_folder/data/bait_MYBs.fasta \
--info /my_folder/data/bait_MYBs.txt \
--out /my_folder/test_analysis/ \
--subject /my_folder/data/Croton_tiglium_peptides.fasta \
--mode raxml \
--refmybs /my_folder/data/AthRefMYBs.txt \
--raxml /my_folder/bin/RAxML9/raxml-ng \
```

## Result files ##

While running the script, several temporary files are generated, but there are not important. Please only look at files located in the 'RESULTS' folder. Here is a list of the possible result files. The generation of some files depends on setting specific options. Therefore, it is possible that you do not see all these files in your result folder. If you specified a prefix (--name), all file names with start with that prefix and thus differ in this respect from the file names listed here.

* **00_documentation.txt** this file stores the command used to run the analysis and also the versions of the tools used for the individual steps. Information contained in this file should be included in the method section when writing a publication about the results of this analysis.

* **01_initial_candidates.fasta** this FASTA file contains all sequences in your subject set that have some similarity to the bait MYB sequences. This selection depends on the setting that you used for the BLASTp search. Usually, this is a very large collection of sequences that have only small similarity to MYBs. Many of these sequences will be discarded in the next step, but this collection is important to achieve high sensitivity. Truncated MYBs could be of high interest and should not be discarded due to short alignments.

* **02a_clean_MYBs.fasta** this FASTA file contains all MYB candidates that fall into the clade of the bait MYB sequences i.e. all MYB-like sequences should be removed at this stage. The accuracy might be low at the basis of the MYB clade, but this approach should allow the inclusion of truncated MYBs that are sitting deep in one of the MYB subgroups.

* **02b_in_out_MYB_analysis_results.txt** this file provides detailed information about the evidence for inclusion or exclusion, respectively, of each initial MYB candidate. The decision is based on phylogenetic distances between the sequences in your subject set and the bait MYBs supplied. If a candidate is close to many bait sequences that are considered bona fide MYBs, the sequence will be included. If a candidate is located in a clade with MYB-like bait sequences, it will NOT be considered a MYB.
	* Score = Ratio of ingroup matches to outgroup matches.
	* IngroupMatches = Number of nearest neighbours (only baits considered) in the tree that are ingroup members.
	* OutgroupMatches = Number of nearest neighbours  (only baits considered) in the tree that are outgroup members.

* **03a_group_around_ref_MYBs.txt** this file is a table which assigns all MYB candidates in the subject species to the provides reference sequences. This file is only generated if such reference sequences are provided (--refmybs).
 
* **03b_new_2_ref_myb_mapping_file.txt** this file is table which assigns the best fitting reference MYB (hopefully the ortholog) to each of the candidate MYBs in your subject data set. This file is only generated if such reference sequences are provided (--refmybs).
	* EdgeDistance = Number of edges between two leaves (sequences) on the tree.
	* PatristicDistance = Cumulative length of all nodes between two leaves (sequences) on the tree.

* **04a_MYB_domain_check.txt** this table contains the results of a search for the MYB domains (1R, R23R, 3R, and others) in all of the clean MYB candidates. Regulator expressions are used to identify 3R, R2R3, or 1R motifs. If none of them is detected, the candidate is classified as other/pseudo-MYB. Here are the regular expressions used for the three domains:

	* R1 = '\w{3,4}W\w{17,21}W\w{17,21}W\w{5,8}'

	* R2 = '\w{5}[WF]{1}\w{18,21}W\w{15,27}[WY]{1}\w{4}'	#F at pos1 and Y at pos3 are very rare

	* R3 = "\w{5}[WLIMF]{1}\w{14,21}W\w{17,21}[WYF]{1}\w{4}"	#diversity at pos1 is high, but po2 and pos3 are conserved

	These MYB domain patterns are based on [Feng et al., 2017](https://doi.org/10.1093/gbe/evx056), [Du et al., 2015](https://doi.org/10.1038/srep11037), and [Pucker et al., 2020](https://doi.org/10.1371/journal.pone.0239275).


* **04b_motif_check.txt** this table contains the results of a search for a set of provided motifs in all of the clean MYB candidates. This file is only generated if a file with motifs was supplied (--motif).

* **04c_MYB_domain_check.fasta** this file contains all MYB domains identified in the clean MYB candidates.

* **04d_MYB_domain_check.doc.txt** this file quanifies the number of clean MYB candidates containing a certain domain. Statistics collected in this file can be combined into a summary table if more than one subject file is analyzed.

* **05** this tree is generated based on all MYB bait sequences and the identified clean MYB candidates of your subject file. The tool used for the tree construction depends on your settings.

* **06** this tree is generated based on all _Arabidopsis thaliana_ and the identified clean MYB candidates of your subject file. The tool used for the tree construction depends on your settings.

* **07a_repr_MYBs.txt** this file is generated if you choose to represent groups of paralogs by only one representative sequence (--collapse). 

* **07b_repr_MYBs.fasta** this file is generated if you choose to represent groups of paralogs by only one representative sequence (--collapse). The sequences contained in this FASTA file represent one group of paralogos. All singletons are included as well.

* **08a_repr_ath_MYBs.fasta** this file contains all representative sequences (see #7) and also the _Arabidopsis thaliana_ MYB sequences. Therefore, this file is only generated if groups of paralogs are collapsed to one representative sequence (--collapse) and the necessary _Arabidopsis thaliana_ MYB sequences are provided (--ath).

* **08b** this tree is generated based on all _Arabidopsis thaliana_ and the representative MYB candidate sequences. The tool used for the tree construction depends on your settings. This tree file is only generated if the file '08a' was generated.

* **08d_MYB_domain_check.txt** please seee the full description at #04a for details. This file is the result of the same analysis, but restricted to the representative sequences described above.

* **08e_MYB_domain_check.fasta** please seee the full description at #04c for details. This file is the result of the same analysis, but restricted to the representative sequences described above.

* **08f_MYB_domain_check.doc.txt** please seee the full description at #04d for details. This file is the result of the same analysis, but restricted to the representative sequences described above.

* **4_domain_detection_summary.txt** this file is a summary of the MYB domain check results described in #04b across all analyzed species. This file is only generated if more than one species are investigated.

* **8_domain_detection_summary.txt** this file is a summary of the MYB domain check results described in #08f across all analyzed species. It is similar to the one described above, but restricted to the representative sequences of paralog groups. This file is only generated if more than one species are investigated.



## Data sets ##

The MYB bait sequences are intended to cover a broad phylogenetic range. The bait sequence set comprises the MYB doamins of:

_Arabidopsis thaliana_ ([Stracke et al., 2001](https://doi.org/10.1016/S1369-5266(00)00199-0))

_Beta vulgaris_ ([Stracke et al., 2014](https://doi.org/10.1186/s12870-014-0249-8))

_Musa acuminata_ ([Pucker et al., 2020](https://doi.org/10.1371/journal.pone.0239275))

and a collection by [Du et al., 2015](https://doi.org/10.1038/srep11037): _Medicago truncatula_, _Populus trichocarpa_, _Citrus sinensis_, _Vitis vinifera_, _Solanum lycopersicum_, _Solanum tuberosum_, _Aquilegia coerulea_, _Oryza sativa_, _Zea mays_, _Amborella trichopoda_, _Picea abies_, _Selaginella moellendorfii_, _Physcomitrella patens_, _Chlamydomonas reinhardtii_, _Volvox carteri_, _Micromonas pusilla_, _Ostreococcus lucimarinus_, and _Cyanidioschyzon merolae_


## References ##

Pucker, B., Reiher, F. and Schilbert, H. M. Automatic identification of players in the flavonoid biosynthesis with application on the biomedicinal plant _Croton tiglium_. Plants 2020, 9, 1103. doi:[10.3390/plants9091103](https://www.mdpi.com/2223-7747/9/9/1103/htm).

Pucker B., Pandey A., Weisshaar B. and Stracke R. (2020). The R2R3-MYB gene family in banana (_Musa acuminata_): genome-wide identification, classification and expression patterns. PLOS ONE 15(10): e0239275. doi:[10.1371/journal.pone.0239275](https://doi.org/10.1371/journal.pone.0239275).

Stracke R., Werber M., Weisshaar B. (2001). The R2R3-MYB gene family in _Arabidopsis thaliana_. Current Opinion in Plant Biology. 2001;4(5):447-456. 

Stracke, R., Holtgräwe, D., Schneider, J., Pucker, B., Sörensen, T.R., and Weisshaar, B. (2014). Genome-wide identification and characterisation of R2R3-MYB genes in sugar beet (Beta vulgaris). BMC Plant Biol. 14: 249. doi:[10.1186/s12870-014-0249-8](https://doi.org/10.1186/s12870-014-0249-8).

Sukumaran J, Holder MT. DendroPy: a Python library for phylogenetic computing. Bioinformatics. 2010 Jun 15;26(12):1569-71. doi: [10.1093/bioinformatics/btq228](https://doi.org/10.1093/bioinformatics/btq228). Epub 2010 Apr 25. PMID: 20421198.

Du, H., Liang, Z., Zhao, S. et al. The Evolutionary History of R2R3-MYB Proteins Across 50 Eukaryotes: New Insights Into Subfamily Classification and Expansion. Sci Rep 5, 11037 (2015). doi:[10.1038/srep11037](https://doi.org/10.1038/srep11037)

