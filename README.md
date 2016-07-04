# SPLITREADER

SPLITREADER is a bioinformatic pipeline dedicated to the discovery of non-reference TE insertions with Target Site Duplications (TSDs) through the use of Illumina short genome sequence reads and a reference genome.

##Table of contents
- [**Installation**](#installation)
 - [Dependencies](#Dependencies)
- [**SPLITREADER main steps**](#steps)
- [**Input data**](#inputs)

- [**Usage**](#usage)
 - [Prepare TEs](#prepare-tes)
 - [Detect TE insertions](#detect-te-insertions)
- [**Output**](#Output)
- [**Cite**](#cite)



##Installation
You can either download the code by clicking the ZIP link on this webpage or clone the project using:

	git clone https://github.com/LeanQ/SPLITREADER.git

###Dependencies

SPLITREADER requieres the following dependencies:

* SAMtools (v1.2 or higher) (http://samtools.sourceforge.net/)
* Bowtie2 (v2.2.9 or higher) https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/)
* Picard tools (> Java 1.8) (https://github.com/broadinstitute/picard/releases/tag/2.4.1)
* bedtools (v2.20.1 or higher) (https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz)

##SPLITREADER main steps

SPLITREADER pipeline consists in four steps: 

1- Extraction of reads not mapping to the reference genome (using SAM flag 4) 

2- Forced mapping to a collection of reference TE sequences (constructed using prepare_TEs.sh as indicated in (#prepare-tes) section). SPLITREADER then identify all reads with one end (≥20nt) mapping to a TE extremity (by locating reads where the CIGAR string starting or ending with ‘S’ characterwith a value equal or greater to 20)

3- Reads are then recursively soft clipped by 1nt from one end and mapped to the reference genome using Bowtie2 until the soft-clipped read length reached 20nt.

4- Read clusters composed by reads clipped from the same extremity and overlapping with read clusters composed of reads clipped from the other extremity were taken to indicate the presence of a bona fide TE insertion only if the size of the overlap was equal or less than 2-fold longer than that reported for the corresponding TE family.

##Input data

This pipeline requieres as an input:

1. A BAM file of reads mapped to the reference genome. 
2. Tab-delimited file containing TE coordiantes across the genome and the expected size for the TSDs. The file should have the following format: 

	- Chromosome name
	- TE start
	- TE end
	- TE_ID
	- TSD

*NOTE: Better performance is obtained if the expected size for the TSDs are provided. If no TSD is provided, it will be set to 3.*

##Usage

Use the test data in order to quickly test the scripts. To do so, one can edit the configuration file with the correct path to the dependecies and run the following command:

	 bash SPLITREADER-beta1.2.sh -c SPLITREADER_configuration_file_example.txt

##Prepare TEs

Before runing TE detection, one must create a Bowtie2 index for each input TE. 

    bash prepare_TEs.sh -i <Input TE coordiantes described above> 
						-g <Genome.fa> 
						-d <directory used to store created TE indexes. If it contains a collection of indexed TE, the script will check if the index already exists> 
						-p </path/to/bowtie2-build> 

##Detect TE insertions
Once the TEs are indexed, one can edit the configuration file and procceed with the discovery of non-reference TE insertions using the main script: 

	bash SPLITREADER-beta1.2.sh -c SPLITREADER_configuration_file.txt

##Output

SPLITREADER outputs the location of the detected non-reference TE insertions and a BAM file for visualisation purposes. The output file have the following format: 

	- Chromosome name
	- Insertion start
	- Insertion end
	- Reconstructed TSD size
	- # Reads supporting the 5'end of the TE insertion 
	- # Reads supporting the 3'end of the TE insertion 

##Cite

If you use this software, please cite:

The Arabidopsis thaliana mobilome and its impact at the species level. Quadrana L, Bortolini Silveira A, Mayhew GF, LeBlanc C, Martienssen RA, Jeddeloh JA, Colot V. 
Elife. 2016 Jun 3;5. pii: e15716. doi: 10.7554/eLife.15716. PubMed PMID: 27258693.





