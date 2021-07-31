# MYB_annotator

This tool allows an automatic identification and analysis of MYBs in a CDS or peptide sequence collection. Since the quality of structural annotation has substantially improved during the last years, the focus can now shift towards functional annotation. There are numerous publications about R2R3-MYBs in various plant species. This tools allows the automatic investigation of the MYB gene family in any new species. Results are a FASTA file containing the bona fide MYB sequences and a summary table with additional details.

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
					
--mafft <PATH_TO_MAFFT>[mafft]
--blastp <PATH_TO_AND_INCLUDING_BINARY>[blastp]
--makeblastdb <PATH_TO_AND_INCLUDING_BINARY>[makeblastdb]			
--fasttree <PATH_TO_FASTTREE>[fasttree]
--raxml <PATH_TO_RAXML>[raxml]				
					
--simcutp <BLASTP_SIMILARITY_CUTOFF>[0.8]
--poscutp <BLASTP_HIT_NUMBER_PER_BAIT_CUTOFF>[100]
--lencutp	<BLASTP_MIN_LENGTH_CUTOFF>[50]
```

`--gff` specifies a GFF3 file. It is important that the IDs of mRNAs in this file are matching the IDs of the coding sequences in the FASTA file.



## References ##

Pucker, B., Reiher, F. and Schilbert, H. M. Automatic identification of players in the flavonoid biosynthesis with application on the biomedicinal plant _Croton tiglium_. Plants 2020, 9, 1103. doi:[10.3390/plants9091103](https://www.mdpi.com/2223-7747/9/9/1103/htm)

Pucker B., Pandey A., Weisshaar B. and Stracke R. (2020). The R2R3-MYB gene family in banana (_Musa acuminata_): genome-wide identification, classification and expression patterns. PLOS ONE 15(10): e0239275. doi:[10.1371/journal.pone.0239275](https://doi.org/10.1371/journal.pone.0239275)


