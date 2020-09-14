# SPLITREADER

SPLITREADER is a bioinformatic pipeline dedicated to the discovery of non-reference TE insertions with Target Site Duplications (TSDs) through the use of Illumina short genome sequence reads and a reference genome.

## Table of contents
- [SPLITREADER](#splitreader)
  * [Installation](#installation)
    + [Dependencies](#dependencies)
  * [SPLITREADER main steps](#splitreader-main-steps)
  * [Input data](#input-data)
  * [Usage](#usage)
	+ [Prepare index files for reference genome](#prepare-index-files-for-reference-genome)
	+ [Prepare TEs](#prepare-tes)
	+ [Detect TE insertions](#detect-te-insertions)
	+ [Output](#output)
  * [Cite](#cite)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>


## Installation
You can either download the code by clicking the ZIP link on this webpage or clone the project using:

```
git clone https://github.com/LeanQ/SPLITREADER.git
```

### Dependencies

SPLITREADER requires the following softwares:

* SAMtools (v1.2 or higher) (http://samtools.sourceforge.net/)
* Bowtie2 (v2.2.9 or higher) https://sourceforge.net/projects/bowtie-bio/files/bowtie2/)
* Picard tools (> Java 1.8) (https://broadinstitute.github.io/picard/)
* bedtools (v2.20.1 or higher) (https://bedtools.readthedocs.io/en/latest/content/installation.html)

## SPLITREADER main steps

SPLITREADER pipeline consists in four steps: 

1- Extraction of reads not mapping to the reference genome (using SAM flag 4) 

2- Forced mapping to a collection of reference TE sequences (constructed using `prepare_TEs.sh` as indicated [below](#prepare-tes)). SPLITREADER then identify all reads with one end (>=20nt) mapping to a TE extremity (by locating reads where the CIGAR string starting or ending with 'S' character with a value equal or greater to 20)

3- Reads are then recursively soft clipped by 1nt from one end and mapped to the reference genome using Bowtie2 until the soft-clipped read length reached 20nt.

4- Read clusters composed by reads clipped from the same extremity and overlapping with read clusters composed of reads clipped from the other extremity were taken to indicate the presence of a bona fide TE insertion only if the size of the overlap was equal or less than 2-fold longer than that reported for the corresponding TE family.

## Input data

This pipeline requieres as an input:

1. A BAM file of reads mapped to the reference genome
2. A tab-delimited file containing TE coordinates across the genome and the expected size for the TSDs. The file should have the following format: 

	- Chromosome name (should be the same nomenclature as in the reference genome fasta file)
	- TE start
	- TE end
	- TE_ID
	- TSD

For instance as indicated in the `./Test_data/TE_list/test_TE_coordinates.bed` file:

```
cat test_TE_coordinates.bed
#Chr start end  TE_ID TSD
Chr1 21524995 21525295 TE1  5
Chr1 21529551 21529851 TE1  5
Chr3 13369174 13369474 TE1  5
Chr3 13373808 13374108 TE1  5
Chr3 22059535 22059835 TE1  5
Chr3 22064029 22064329 TE1  5
Chr3 22695566 22695866 TE1  5
Chr3 22700222 22700522 TE1  5
Chr5 4208083 4208383 TE1 5
Chr5 4212784 4213084 TE1 5
```

*NOTE: Better performance is obtained if the expected size for the TSDs are provided. If no TSD is provided, it will be set to 3.*

## Usage

### Prepare index files for reference genome
Use the test data in order to quickly test the scripts. To do so, one needs to download the TAIR10 *Arabidopsis thaliana* fasta file and generate the bowtie2 indexes.

```
# Download the fasta file
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas

# Add Chr as prefix for each chromosome
sed -i 's/>\([1-5]\)/>Chr\1/' TAIR10_chr_all.fas

# Index the file with bowtie2 (use prefix TAIR10 for the index files)
bowtie2-build TAIR10_chr_all.fas TAIR10
```

The generated index files should be:

```
TAIR10.1.bt2
TAIR10.2.bt2
TAIR10.3.bt2
TAIR10.4.bt2
TAIR10.rev.1.bt2
TAIR10.rev.2.bt2
```

Edit the configuration file with the correct path to the dependencies `SPLITREADER_configuration_file_example.txt`.

### Prepare TEs

Before runing TE detection, one must create a Bowtie2 index for each input TE. 

```
bash prepare_TEs.sh -i <Input TE coordiantes described above>
	-g <Genome.fa>
	-d <directory used to store created TE indexes. If it contains a collection of indexed TE, the script will check if the index already exists>
	-p </path/to/bowtie2-build> 
```

In our example, it will be:

```
bash prepare_TEs.sh -i ./Test_data/TE_list/test_TE_coordinates.bed \
	-g ./TAIR10.fa \
	-d ./Test_data/TE_indexes  \
	-p /path/to/bowtie2-build
```

### Detect TE insertions
Once the TEs are indexed, one can edit the configuration file and procceed with the discovery of non-reference TE insertions using the main script: 

```
bash SPLITREADER-beta1.2.sh -c SPLITREADER_configuration_file.txt
```

### Output

SPLITREADER outputs the location of the detected non-reference TE insertions and a BAM file for visualisation purposes. The output file have the following format: 

	- Chromosome name
	- Insertion start
	- Insertion end
	- Reconstructed TSD size
	- # Reads supporting the 5'end of the TE insertion 
	- # Reads supporting the 3'end of the TE insertion 

For our example output, the expected output files are:

```
test-TE1-insertion-sites.bed
test-TE1-split.bam
test-TE1-split.bam.bai
```

And the content of `test-TE1-insertion-sites.bed` should be:

```
cat test-TE1-insertion-sites.bed
Chr1 23655713 23655723 TE1  9  8 3
Chr1 23655713 23655723 TE1  9  8 4
```

## Cite

If you use this software, please cite:

The Arabidopsis thaliana mobilome and its impact at the species level. Quadrana L, Bortolini Silveira A, Mayhew GF, LeBlanc C, Martienssen RA, Jeddeloh JA, Colot V. 
Elife. 2016 Jun 3;5. pii: e15716. doi: 10.7554/eLife.15716. PubMed PMID: 27258693. (https://elifesciences.org/content/5/e15716)





