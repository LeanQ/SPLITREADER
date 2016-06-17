# SPLITREADER

SPLITREADER is a bioinformatic pipeline dedicated to the discovery of non-reference TE insertions with Target Site Duplications (TSDs) through the use of short-reads sequencies and a reference genome.

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
