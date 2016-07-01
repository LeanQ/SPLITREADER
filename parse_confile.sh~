#!/bin/bash


:<<'comment'
display_usage() { 
        echo -e "\nparse_confile.sh: parse SPLITREADER configuration file.\n"
        echo -e "Usage: bash parse_confile.sh -c SP_config.conf\n" 
        echo "Options:"
        echo 
} 

# if less than two arguments supplied, display usage 
if [  $# -eq 0 ] ; then  display_usage; exit 1; fi 

#if [[ ! -s $CONFIG_FILE ]]; then >&2 echo -e "[ERROR]: Empty/non-existing parameter file or parameter file not correctly provided (use -c option)."; exit 1; fi
while getopts "c:" opt
do

case $opt in

c) CONFIG_FILE=$OPTARG ;;
?) echo "Unknown option -$OPTARG ";;

esac
done
comment
#CONFIG_FILE=$1
CONFIG_FILE=SPLITREADER_configuration_file_example.txt


########################### Parsing conig file
grep -v "#" $CONFIG_FILE | awk '/./' > tmp.conf

typeset -A config # init array
config=( # set default values in config array
[OutputDir]=""
[TmpDir]=""
[in]=""
[listTE]=""
[LENGTH]=""
[READS]=""
[maxcov]=""
[CORES]=""
[SequencesDir]=""
[GenomeIndexFile]=""
[Bowtie2Dir]=""
[samtoolsDir]=""
[bedtoolsdir]=""
[picardDir]=""
)

while read line
do
    if echo $line | grep -F = &>/dev/null
    then
        varname=$(echo "$line" | cut -d '=' -f 1)
        config[$varname]=$(echo "$line" | cut -d '=' -f 2-)
    fi
done < tmp.conf


OutputDir=${config[OutputDir]}
TmpDir=${config[TmpDir]}
in=${config[in]}
listTE=${config[listTE]}
LENGTH=${config[LENGTH]}
READS=${config[READS]}
maxcov=${config[maxcov]}
CORES=${config[CORES]}
SequencesDir=${config[SequencesDir]}
GenomeIndexFile=${config[GenomeIndexFile]}
Bowtie2Dir=${config[Bowtie2Dir]}
samtoolsDir=${config[samtoolsDir]}
bedtoolsdir=${config[bedtoolsdir]}
picardDir=${config[picardDir]}
############################ END parsing config file


cat > tmp.conf <<EOL

###########################
SPLITREADER RUN PARAMETERS
###########################

Working directories:
===================
Output dir: $OutputDir
Temporary file: $TmpDir

Inputs and main parameters:
==========================
Input BAM file: $in
List of TEs: $listTE

Read length (bp): $LENGTH
Maximum coverage: $maxcov
Number of threads: $CORES

Indexes:
========
Bowtie2 index (reference genome): $GenomeIndexFile
Directory of Bowtie2 indexes for each TE: $SequencesDir

Tools:
======
Path to Bowtie2: $Bowtie2Dir
Path to Samtools: $samtoolsDir
Path to Bedtools: $bedtoolsdir
Path to Picard tools: $picardDir
EOL

cat tmp.conf
rm -f tmp.conf

