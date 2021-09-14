# Automatic annotation of MYBs

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
--cdsinput      STR     Changes expected input to CDS
					
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

`--subjectdir` specifies a folder with FASTA files. All files need to contain sequences of the same type (PEP or CDS). Supported file extensions are '.fasta', '.fas', '.fa', '.FASTA', '.FAS', and '.FA'.

`--mode` specifies the tool for the construction of phylogenetic trees. [RAxML](https://doi.org/10.1093/bioinformatics/btz305) ('raxml') or [FastTree2](https://doi.org/10.1371/journal.pone.0009490) ('fasttree') can be specified. Default is 'fasttree'. Running the tree construction with RAxML is expected to result in slightly more accurate results, but would increase the run time from minutes to days in many cases.

`--refmybs` specifies a text file with one MYB ID per line. All IDs listed in this file need to be present in the baits FASTA file and the info file. These IDs can be used for the functional annotation of the newly identified MYBs. This assignment is performed in two ways leading to two different output tables. The functional annotation is lifted from previously characterized sequences to newly identified orthologs. 

`--cpu` specifies threads to use for the BLASTp search and also for the RAxML tree construction. Default: 4.

`--cdsinput` changes the input expectation from PEP to CDS. The file needs to be supplied through the `--subject`. This flag can be included at any place in the command and does not need a value.

`--mafft` specifies the full path to [MAFFT](https://doi.org/10.1093/molbev/mst010) including the binary.

`--blastp` specifies the full path to [BLASTp](https://dx.doi.org/10.1093%2Fnar%2F25.17.3389) including the binary.

`--fasttree` specifies the full path to [FastTree2](https://doi.org/10.1371/journal.pone.0009490). This is required if fasttree is not included in the PATH variable. If in doubt, please specify the full path.

`--raxml` specifies the full path to [RAxML](https://doi.org/10.1093/bioinformatics/btz305).

`--simcutp` specifies the minimal similarity of BLASTp hits to be considered in the initial identification of MYB candidates. Default: 0.8 (80%).

`--poscutp` specifies the maximal number of hits (possibilities) per bait sequence in the BLASTp analysis. Default: 100.

`--lencutp` specifies the minimal alignment length of BLASTp hits to be considered in the initial identification of MYB candidates. Default: 50.



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

