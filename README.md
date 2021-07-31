# Automatic annotation of MYBs

## Background ##
This tool allows an automatic identification and analysis of MYBs in a CDS or peptide sequence collection. Since the quality of structural annotation has substantially improved during the last years, the focus can now shift towards functional annotation. RNA-seq and long read sequencing of transcripts were fundamental for this development. There are numerous publications about R2R3-MYBs in various plant species. This tools allows the automatic investigation of the MYB gene family in any new species. Results are a FASTA file containing the bona fide MYB sequences and a summary table with additional details.

## Installation ##

There is no installation of this tool required. Downloading and executing the script is sufficient. Most required modules should be included in the initial Python installation. Requirements are (Linux installation instructions):

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
--subject       STR     Subject file
--subjectdir    STR     Subject dir
					
Optional:
--mode          STR     Tool for tree construction (fasttree|raxml)[raxml]
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

`--info` specifies a text file that contains a table with all MYB IDs (matching the baits FASTA file) in the first column. Columns are TAB-separted. The second column contains an assignment to the ingroup ('in') or outgroup ('out'). Ingroup and outgroup baits are used to classify MYB candidates as _bona fide_ MYBs or MYB-likes.

`--out` specifies the output folder. This folder will be created if it does not exist already. All files generated during the analysis will be placed in this folder. There is a subfolder called 'RESULTS' that contains the important final files.

`--subject` specifies the FASTA file that should be analysed. The default expectation is that this file contains peptide sequences. However, it is also possible to supply a collection of coding sequences (CDS) if the `--cdsinput` flag is set.

`--subjectdir` specifies a folder with FASTA files. All files need to contain sequences of the same type (PEP or CDS).

`--mode` specifies the tool for the construction of phylogenetic trees. [RAxML](https://doi.org/10.1093/bioinformatics/btz305) ('raxml') or [FastTree2](https://doi.org/10.1371/journal.pone.0009490) ('fasttree') can be specified. Default is 'raxml'.

`--refmybs` specifies a text file with one MYB ID per line. All IDs listed in this file need to be present in the baits FASTA file and the info file. These IDs can be used for the functional annotation of the newly identified MYBs. This assignment is performed in two ways leading to two different output tables.

`--cpu` specifies threads to use for the BLASTp search and also for the RAxML tree construction. Default: 4.

`--cdsinput` changes the input expectation from PEP to CDS. The file needs to be supplied through the `--subject`. This flag can be included at any place in the command.

`--mafft` specifies the full path to [MAFFT](https://doi.org/10.1093/molbev/mst010) including the binary.

`--blastp` specifies the full path to [BLASTp](https://dx.doi.org/10.1093%2Fnar%2F25.17.3389) including the binary.

`--fasttree` specifies the full path to [FastTree2](https://doi.org/10.1371/journal.pone.0009490).

`--raxml` specifies the full path to [RAxML](https://doi.org/10.1093/bioinformatics/btz305).

`--simcutp` specifies the minimal similarity of BLASTp hits to be considered in the initial identification of MYB candidates. Default: 0.8 (80%).

`--poscutp` specifies the maximal number of hits per bait sequence in the BLASTp analysis. Default: 100.

`--lencutp` specifies the minimal alignment length of BLASTp hits to be considered in the initial identification of MYB candidates. Default: 50.



## Example ##



## References ##

Pucker, B., Reiher, F. and Schilbert, H. M. Automatic identification of players in the flavonoid biosynthesis with application on the biomedicinal plant _Croton tiglium_. Plants 2020, 9, 1103. doi:[10.3390/plants9091103](https://www.mdpi.com/2223-7747/9/9/1103/htm)

Pucker B., Pandey A., Weisshaar B. and Stracke R. (2020). The R2R3-MYB gene family in banana (_Musa acuminata_): genome-wide identification, classification and expression patterns. PLOS ONE 15(10): e0239275. doi:[10.1371/journal.pone.0239275](https://doi.org/10.1371/journal.pone.0239275)


