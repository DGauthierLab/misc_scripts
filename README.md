# misc_scripts
Basic script for identifying regions of deletion in bacterial genomes.  
Currently set for >=750 bp regions with coverage <5X

Requires:
samtools http://www.htslib.org/
SHRiMP http://compbio.cs.toronto.edu/shrimp/
bedtools https://bedtools.readthedocs.io/en/latest/

Requires paired-end illumina reads, a .fasta reference genome, and a genomefile with format: 

\<chromosome name\>  \<length\>
