# SPLITREADER

SPLITREADER is a bioinformatic pipeline dedicated to the discovery of non-reference TE insertions with Target Site Duplications (TSDs) through the use of Illumina short genome sequence reads and a reference genome.

This pipeline requieres as an input a .bam file of mapped reads in the reference genome as well as a file of TE annotations. Better performance is obtained if the expected size for the TSDs are also provided

If you use this software, please cite:
The Arabidopsis thaliana mobilome and its impact at the species level. Quadrana L, Bortolini Silveira A, Mayhew GF, LeBlanc C, Martienssen RA, Jeddeloh JA, Colot V. 
Elife. 2016 Jun 3;5. pii: e15716. doi: 10.7554/eLife.15716. PubMed PMID: 27258693.


Dependencies:
SPLITREADER requieres the following dependencies:
SAMtools (http://samtools.sourceforge.net/)
Bowtie2 (https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/)
Picard tools (> Java 1.8) (https://github.com/broadinstitute/picard/releases/tag/2.4.1)
bedtools (https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz)


SPLITREADER pipeline consistes in four steps: 
1- Extraction of reads not mapping to the reference genome (containing the SAM flag 4)\n 
2- Forced mapping to a collection collection of 5’ and 3’ TE sequence extremities (300bp) obtained from the reference genome (using the module TEXTRACTION). SPLITREADER then identify all reads with one end (≥20nt) mapping to a TE extremity (by locating reads where the CIGAR string contained the ‘S’ character with a value equal or greater to 20)
3- SPLITREADER are recursively soft clipped from one end by 1nt and mapped to the provided reference genome using Bowtie2 until the soft-clipped read length reached 20nt.
4- Read clusters composed of predefined number of reads clipped from the same extremity and overlapping with read clusters composed of reads clipped from the other extremity were taken to indicate the presence of a bona fide TE insertion only if the size of the overlap was equal or less than 2-fold longer than that reported for TSDs for the corresponding TE family
