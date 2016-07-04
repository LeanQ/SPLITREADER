#!/bin/bash

############### USAGE/HELP

#TODO: 
#in test data put only TE_coordinates.bed
#create folder TAIR10_TE_indexes


display_usage() { 
        echo -e "\nprepare_TEs.sh: Preparation of TEs sequences and indexes.\n"
        echo -e "Usage: bash prepare_TEs.sh -i <TE_coordiantes.bed> -g <Genome.fa> -d <directory of existing TE indexed provided by SPLITREADER> -p </path/to/bowtie2-build>\n" 
        echo 
		echo -e "Example command: bash prepare_TEs.sh -i Test_data/TE_list/test_TE_coordinates.bed -g Genome.fa -d Test_data/TE_indexes -p /usr/local/bin/bowtie2-build"
		echo "Contact:quadrana(ar)biologie.ens.fr"
		echo
} 

# if no arguments supplied, display usage 
if [  $# -eq 0 ] ; then  display_usage; exit 1; fi 
# check whether user had supplied -h or --help . If yes display usage 
if [[ ( $1 == "--help") ||  $1 == "-h" ]]; then display_usage; exit 0; fi 


while getopts "i:g:d:p:" opt
do

case $opt in

i) TEs=$OPTARG ;;
g) GENOME=$OPTARG ;;
d) INDEXES=$OPTARG ;;
p) BOWTIE2=$OPTARG ;;
?) echo "Unknown option -$OPTARG ";;

esac
done


#### Some tests: 
#if TE file not correctly formated (we need 5 cols : chr-start-end-ID-TSD)
a=`grep -v "#" $TEs| awk '/./' |wc -l`
b=`grep -v "#" $TEs | awk '/./ && NF==5 {print}'|wc -l`
if [ $a != $b ]; then >&2 echo "[ERROR]: $TEs file is not formated correctly. Number of columns must be 5." ; exit 1; fi
# if dir indexes dos not exists
if [ ! -d $INDEXES ]; then  >&2 echo "[ERROR]:`file $INDEXES`" ; exit 1; fi
# remove tailing \ if present
INDEXES=`echo "$INDEXES" | sed 's|/| |g'`
##### END tests


#create tmp file to store .fa files
if [ ! -d TMP_FA ]; then mkdir TMP_FA; fi

####################################################################
#STEP1: 
#Extract sub file of coordiantes for each TE based on the ID column
####################################################################

#get list of TE names
grep -v "#" $TEs | awk '/./' | cut -f 4 | sort -u > TE_names.txt

# get the coordiantes for each TE
cat TE_names.txt | while read ID
do

echo "----> Processing $ID "

# Check if an index for ID already exists
if [ -f  $INDEXES/$ID.1.bt2 ]; then 

echo "Bowtie2 index for $ID already exists in TE_indexes folder. Skip $ID ! "

else
#Extract fastq file
echo "Extracting $ID sequences.."
awk -v ID=$ID '$4 == ID' $TEs | bedtools getfasta  -fi $GENOME -bed - -fo TMP_FA/$ID.fa
echo "Creating Bowtie2 indexes.."
$BOWTIE2 TMP_FA/$ID.fa $INDEXES/$ID > /dev/null 2>&1	
fi
done

rm -rf TE_names.txt TMP_FA


