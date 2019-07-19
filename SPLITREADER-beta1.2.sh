#!/bin/bash 

#############################################
#                                           #  
#                                           #
#            SPLITREADER beta1.2            #
#                                           #
#                                           # 
#############################################

##Questions or comments to quadrana(ar)biologie.ens.fr

############### USAGE/HELP
display_usage() { 
        echo -e "\nSPLITREADER: Split-read pipeline for the identification of non-reference TE insertions with TSDs — .\n"
        echo -e "Usage: bash SPLITREADER-beta1.2.sh -c SPLITREADER_configuration_file_example.txt\n" 
        echo 
		echo "Contact:quadrana(ar)biologie.ens.fr"
		echo
} 

# if less than two arguments supplied, display usage 
if [  $# -eq 0 ] ; then  display_usage; exit 1; fi 
# check whether user had supplied -h or --help . If yes display usage 
if [[ ( $1 == "--help") ||  $1 == "-h" ]]; then display_usage; exit 0; fi 


while getopts "c:" opt
do

case $opt in

c) CONFIG_FILE=$OPTARG ;;
?) echo "Unknown option -$OPTARG ";;

esac
done
#CONFIG_FILE=$1
#CONFIG_FILE=SPLITREADER_configuration_file_example.txt
if [[ ! -s $CONFIG_FILE ]]; then >&2 echo -e "[ERROR]: Empty/non-existing parameter file or parameter file not correctly provided (use -c option)."; exit 1; fi




########################### Parsing conig file
grep -v "#" $CONFIG_FILE | awk '/./' > tmp.conf

typeset -A config # init array
config=( # set default values in config array
[OutputDir]=""
[Tmp]=""
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
Tmp=${config[Tmp]}
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
Temporary directory: $Tmp

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

###############################
END SPLITREADER RUN PARAMETERS
###############################

EOL

cat tmp.conf
rm -f tmp.conf

##################### MAIN ###############################
#Get bam name
tmp=`basename $in`
BAMname=${tmp%.bam}

# PID of this batch
IDPID=$$
# Temporary directory
TmpDir=$Tmp/QD-$IDPID
mkdir -p $TmpDir
  
#Extracting unmapped reads
echo "Extracting unmapped reads from $in"

# Check if bam file contain SE or PE data (0 if SE, 1 if PE) 
pe=`$samtoolsDir/samtools view -c -f 1 $in | awk '{print $1}'` 

# Get unmapped reads from bam file
$samtoolsDir/samtools view -f 4 -u $in > $TmpDir/$BAMname.bam 2>> $TmpDir/log.txt

# Convert the bam files of unmapped reads into fastq files
# If single end data
if [ -z "$pe" ]; then
	java -jar $picardDir/picard.jar SamToFastq I=$TmpDir/$BAMname.bam \
	   	FASTQ=$TmpDir/$BAMname.fastq 2>> $TmpDir/log.txt

# If paired end data
else
	java -jar $picardDir/picard.jar SamToFastq I=$TmpDir/$BAMname.bam \
	   	FASTQ=$TmpDir/$BAMname.1.fastq SECOND_END_FASTQ=$TmpDir/$BAMname.2.fastq 2>> $TmpDir/log.txt
	
	# Concatenate the 2 mates of each unmapped read
	cat $TmpDir/$BAMname.1.fastq $TmpDir/$BAMname.2.fastq > $TmpDir/$BAMname.fastq

	# rm -f $TmpDir/$BAMname.1.fastq
	# rm -f $TmpDir/$BAMname.2.fastq
fi


##############
#Starting the SPLITREADER pipeline for each TE in the TE-indormation.txt file

# Get unique names for TEs
grep -v "#" $listTE | awk '/./' | cut -f 4,5 | sort -u > TE_names.txt


# Loop over each TE name
cat TE_names.txt | while read line ; do
	STARTTIME=$(date +%s)
	
	# Attribute unique name for temporary folder (same PID as parent directory)
	TmpResultsDir=$TmpDir/results-$IDPID
	mkdir -p $TmpResultsDir
	
	# Get name of TE and size of TSD
	TE=`echo $line | awk '{print $1}'`
	TSD=`echo $line | awk '{print $2}'`

	#set TSD to 3 if no TSD is provided		
	if [ -z "$TSD" ]; then 
		$TSD=3
	fi

	echo -e "\n##### RUNNING SPLIT-READ ANALYSIS ON $TE ######\n"    

	# Selecting split-reads by mapping the unmapped reads over TE extremities
	echo "Selecting split-reads"

	$Bowtie2Dir/bowtie2 -x $SequencesDir/$TE -U $TmpDir/$BAMname.fastq \
   -S $TmpResultsDir/$BAMname-$TE.sam --local --very-sensitive \
   --threads $CORES 2>> $TmpDir/log.txt 
        
       
    #############################################################
    ###filter soft-clipped reads with at least 20nt softclipped at 5' or 3' read's end 
    
	
		echo -e "\n#######filter soft-clipped reads with at least 20nt softclipped at 5' or 3' read's end######\n"    
	
	$samtoolsDir/samtools view -F 4 -S $TmpResultsDir/$BAMname-$TE.sam | \
		awk '$6~/^[2-8][0-9]S/ || $6~/[2-8][0-9]S$/ {print $1}'| \
		sed 's/\/2$//' | sed 's/\/1$//' | \
		awk '{print $1"/1""\n"$1"/2"}' | \
		sort -u > $TmpResultsDir/reads.name 2>> $TmpDir/log.txt

	java -jar $picardDir/picard.jar FilterSamReads I=$TmpResultsDir/$BAMname-$TE.sam \
		FILTER=includeReadList READ_LIST_FILE=$TmpResultsDir/reads.name \
		OUTPUT=$TmpResultsDir/$BAMname-$TE-selected.sam \
		2>> $TmpResultsDir/log.txt 2>> $TmpDir/log.txt

	java -jar $picardDir/picard.jar SamToFastq I=$TmpResultsDir/$BAMname-$TE-selected.sam \
		FASTQ=$TmpResultsDir/$BAMname-$TE-split.fastq \
		2>> $TmpResultsDir/log.txt 2>> $TmpDir/log.txt

    ################################
   
    # rm -f $TmpResultsDir/$BAMname-$TE-split.sam


	echo -e "\n##### Estimating max read size ######\n" 

    ### Estimating max read size (If necessary)
	if [ -z "$LENGTH" ]; then
		LENGTH=`awk 'NR%4 == 2 {print length($0)}' $TmpResultsDir/$BAMname-$TE-split.fastq | sort | tail -1 `  
		length=$((LENGTH-20))
		echo "Maximum Read length: $LENGTH [Estimated] "
	else
		length=$((LENGTH-20))
		echo "Maximum Read length: $LENGTH [User defined] "
     fi
    
	
	
   ##Recursive split-reads mapping
   # step 1 for 3' read extremity: begining the loop.
	
	echo -e "\n##### Recursive split-reads mapping ######\n" 
	
	echo -e "\nMapping 5' split-reads on reference genome\n"
	echo -n "Progresssion:["
  
	echo $TmpResultsDir/$BAMname-$TE-split.fastq

	# -U for unpaired reads
	$Bowtie2Dir/bowtie2 -x $GenomeIndexFile -U $TmpResultsDir/$BAMname-$TE-split.fastq \
	  -S $TmpResultsDir/$BAMname-$TE-local.sam --local \
	  --very-sensitive --threads $CORES --quiet 2>> $TmpDir/log.txt

	$samtoolsDir/samtools view -H -S $TmpResultsDir/$BAMname-$TE-local.sam > \
	  $TmpResultsDir/$BAMname-$TE-split-local-up.sam 2>> $TmpDir/log.txt

	cat $TmpResultsDir/$BAMname-$TE-split-local-up.sam > $TmpResultsDir/$BAMname-$TE-split-local-down.sam 

	$samtoolsDir/samtools view -F 4 -S $TmpResultsDir/$BAMname-$TE-local.sam | \
	  awk '$6~/^[0-9][0-9]S/ {print $0}' >> $TmpResultsDir/$BAMname-$TE-split-local-down.sam 2>> $TmpDir/log.txt

	$samtoolsDir/samtools view -F 4 -S $TmpResultsDir/$BAMname-$TE-local.sam | \
	  awk '$6~/[0-9][0-9]S$/ {print $0}' >> $TmpResultsDir/$BAMname-$TE-split-local-up.sam 2>> $TmpDir/log.txt

	$samtoolsDir/samtools view -Sbu -q 5 $TmpResultsDir/$BAMname-$TE-split-local-down.sam | \
	  $samtoolsDir/samtools sort - -o $TmpResultsDir/$BAMname-$TE-split-local-down.bam 2>> $TmpDir/log.txt

	$samtoolsDir/samtools view -Sbu -q 5 $TmpResultsDir/$BAMname-$TE-split-local-up.sam | \
	  $samtoolsDir/samtools sort - -o $TmpResultsDir/$BAMname-$TE-split-local-up.bam 2>> $TmpDir/log.txt

	# Uncommented lines	
	# rm -f $TmpResultsDir/$BAMname-$TE-split-local-up.sam
	# rm -f $TmpResultsDir/$BAMname-$TE-split-local-down.sam
	# rm -f $TmpResultsDir/$BAMname-$TE-local.sam
	  
	echo -e "\n##### Check point 1 ######\n"   
   
   ############################################
   
   length=$(($((LENGTH/2))-1))
   
	$Bowtie2Dir/bowtie2 -x $GenomeIndexFile -U $TmpResultsDir/$BAMname-$TE-split.fastq -S $TmpResultsDir/$BAMname-$TE-splitjunction-5-$length.sam --un $TmpResultsDir/$BAMname-$TE-split-5-$length -5 $length --mp 13 --rdg 8,5 --rfg 8,5 --very-sensitive --quiet 2>> $TmpDir/log.txt

	# This line was commented out before
	#$samtoolsDir/samtools view -Sbu -F 4  $TmpResultsDir/$BAMname-$TE-splitjunction-5-$length.sam | $samtoolsDir/samtools sort - $TmpResultsDir/$BAMname-$TE-splitjunction-5-$length
	# Add output argument -o
	$samtoolsDir/samtools view -Sbu -F 4  $TmpResultsDir/$BAMname-$TE-splitjunction-5-$length.sam | \
		$samtoolsDir/samtools sort - -o $TmpResultsDir/$BAMname-$TE-splitjunction-5-$length.bam
	# rm -f $TmpResultsDir/$BAMname-$TE-splitjunction-5-$length.sam
     
    for ((i=$((length+1)); $i<$((LENGTH-18)); i=$i+1)); do
    echo -n "|" 
      previous=$(($i-1));

      $Bowtie2Dir/bowtie2 -x $GenomeIndexFile -U $TmpResultsDir/$BAMname-$TE-split-5-$previous -S $TmpResultsDir/$BAMname-$TE-splitjunction-5-$i.sam --un $TmpResultsDir/$BAMname-$TE-split-5-$i -5 $i --mp 13 --rdg 8,5 --rfg 8,5 --very-sensitive --quiet 2>> $TmpDir/log.txt
      $samtoolsDir/samtools view -Sbu -F 4 $TmpResultsDir/$BAMname-$TE-splitjunction-5-$i.sam | $samtoolsDir/samtools sort - -o $TmpResultsDir/$BAMname-$TE-splitjunction-5-$i.bam 2>> $TmpDir/log.txt
      # rm -f $TmpResultsDir/$BAMname-$TE-splitjunction-5-$i.sam
      # rm -f $TmpResultsDir/$BAMname-$TE-split-5-$previous
    done
	echo -e "\n##### Check point 2 ######\n" 
    
  $Bowtie2Dir/bowtie2 -x $GenomeIndexFile -U $TmpResultsDir/$BAMname-$TE-split.fastq -S $TmpResultsDir/$BAMname-$TE-splitjunction-3-$length.sam --un $TmpResultsDir/$BAMname-$TE-split-3-$length -3 $length --mp 13 --rdg 8,5 --rfg 8,5 --very-sensitive --quiet 2>> $TmpDir/log.txt
  # Line commented out before
  $samtoolsDir/samtools view -Sbu -F 4  $TmpResultsDir/$BAMname-$TE-splitjunction-3-$length.sam | $samtoolsDir/samtools sort - -o $TmpResultsDir/$BAMname-$TE-splitjunction-3-$length.bam
 #   rm -f $TmpResultsDir/$BAMname-$TE-splitjunction-3-$length.sam
     
    for ((i=$((length+1)); $i<$((LENGTH-18)); i=$i+1)); do
		echo -n "|" 
		previous=$(($i-1));

		$Bowtie2Dir/bowtie2 -x $GenomeIndexFile -U $TmpResultsDir/$BAMname-$TE-split-3-$previous -S $TmpResultsDir/$BAMname-$TE-splitjunction-3-$i.sam --un $TmpResultsDir/$BAMname-$TE-split-3-$i -3 $i --mp 13 --rdg 8,5 --rfg 8,5 --very-sensitive --quiet 2>> $TmpDir/log.txt
		$samtoolsDir/samtools view -Sbu -F 4 $TmpResultsDir/$BAMname-$TE-splitjunction-3-$i.sam | $samtoolsDir/samtools sort - -o $TmpResultsDir/$BAMname-$TE-splitjunction-3-$i.bam 2>> $TmpDir/log.txt
		# rm -f $TmpResultsDir/$BAMname-$TE-splitjunction-3-$i.sam
		# rm -f $TmpResultsDir/$BAMname-$TE-split-3-$previous
	done
	
	echo -n "]"
	echo -e "\n"

	echo -e "\n##### Check point 3 ######\n" 

    # Post-treatment:
    
    # merge all the recursive mappings from either the 3' and 5'
        
    #rm -f $TmpResultsDir/$BAMname-$TE-split.fastq

    # Merging and sorting bam files
    $samtoolsDir/samtools merge -f -u $TmpResultsDir/$BAMname-$TE-split-5.bam $TmpResultsDir/$BAMname-$TE-splitjunction-5-*.bam 2>> $TmpDir/log.txt
    $samtoolsDir/samtools merge -f -u $TmpResultsDir/$BAMname-$TE-split-3.bam $TmpResultsDir/$BAMname-$TE-splitjunction-3-*.bam 2>> $TmpDir/log.txt
    $samtoolsDir/samtools sort $TmpResultsDir/$BAMname-$TE-split-5.bam -o $TmpResultsDir/$BAMname-$TE-split-5.bam 2>> $TmpDir/log.txt
    $samtoolsDir/samtools sort $TmpResultsDir/$BAMname-$TE-split-3.bam -o $TmpResultsDir/$BAMname-$TE-split-3.bam 2>> $TmpDir/log.txt
    
   # rm -f $TmpResultsDir/$BAMname-$TE-splitjunction-[35]-*.bam
    
    # merge reads that were cliped at the 3' and mapped in the + strand with those clipped at the 5' and mapped on the - strand
    # merge reads that were cliped at the 3' and mapped in the - strand with those clipped at the 5' and mapped on the + strand
    
      echo "Searching for reads clusters..."
    $samtoolsDir/samtools view -F 16 -bh $TmpResultsDir/$BAMname-$TE-split-5.bam > $TmpResultsDir/$BAMname-$TE-split-5+.bam 2>> $TmpDir/log.txt 2>> $TmpDir/log.txt
    $samtoolsDir/samtools view -f 16 -bh $TmpResultsDir/$BAMname-$TE-split-5.bam > $TmpResultsDir/$BAMname-$TE-split-5-.bam 2>> $TmpDir/log.txt 2>> $TmpDir/log.txt
    $samtoolsDir/samtools view -F 16 -bh $TmpResultsDir/$BAMname-$TE-split-3.bam > $TmpResultsDir/$BAMname-$TE-split-3+.bam 2>> $TmpDir/log.txt 2>> $TmpDir/log.txt
    $samtoolsDir/samtools view -f 16 -bh $TmpResultsDir/$BAMname-$TE-split-3.bam > $TmpResultsDir/$BAMname-$TE-split-3-.bam 2>> $TmpDir/log.txt 2>> $TmpDir/log.txt

	echo -e "\n##### Check point 4 ######\n" 
    #Merge the 5' and 3' clusters to create the downstream and upstream cluster

   $samtoolsDir/samtools merge -f -u $TmpResultsDir/$BAMname-$TE-up.bam $TmpResultsDir/$BAMname-$TE-split-5-.bam $TmpResultsDir/$BAMname-$TE-split-3+.bam > $TmpResultsDir/$BAMname-$TE-split-local-up.bam 2>> $TmpDir/log.txt
   $samtoolsDir/samtools merge -f -u $TmpResultsDir/$BAMname-$TE-down.bam $TmpResultsDir/$BAMname-$TE-split-5+.bam $TmpResultsDir/$BAMname-$TE-split-3-.bam > $TmpResultsDir/$BAMname-$TE-split-local-down.bam 2>> $TmpDir/log.txt
   samtools sort $TmpResultsDir/$BAMname-$TE-down.bam -o $TmpResultsDir/$BAMname-$TE-down.bam 2>> $TmpDir/log.txt
   samtools sort $TmpResultsDir/$BAMname-$TE-up.bam -o $TmpResultsDir/$BAMname-$TE-up.bam 2>> $TmpDir/log.txt
   
	echo -e "\n##### Check point 5 ######\n" 

    #Calculate the coverage over mapped regions - filter regions according to minimum and maximum read-depth
    $samtoolsDir/samtools depth $TmpResultsDir/$BAMname-$TE-up.bam | awk -v M=$maxcov '$3<(M) {print $1 "\t" $2 "\t"$2"\t"$3}' > $TmpResultsDir/$BAMname-$TE-up.bed 2>> $TmpDir/log.txt
    $samtoolsDir/samtools depth $TmpResultsDir/$BAMname-$TE-down.bam | awk -v M=$maxcov ' $3<(M) {print $1 "\t" $2 "\t"$2"\t"$3}' > $TmpResultsDir/$BAMname-$TE-down.bed 2>> $TmpDir/log.txt

    #merge cluster of covered regions - filter out clusters that are longer than read-length
    sort -k 1,1 -k2,2n $TmpResultsDir/$BAMname-$TE-up.bed | $bedtoolsdir/mergeBed -i stdin -c 4 -o max > $TmpResultsDir/$BAMname-$TE-up-merge.bed 2>> $TmpDir/log.txt
    sort -k 1,1 -k2,2n $TmpResultsDir/$BAMname-$TE-down.bed | $bedtoolsdir/mergeBed -i stdin -c 4 -o max  > $TmpResultsDir/$BAMname-$TE-down-merge.bed 2>> $TmpDir/log.txt
    
     
    # rm -f $TmpResultsDir/$BAMname-$TE-down.bed
    # rm -f $TmpResultsDir/$BAMname-$TE-up.bed
    
    # rm -f $TmpResultsDir/$BAMname-$TE-split-5+.bam
    # rm -f $TmpResultsDir/$BAMname-$TE-split-5-.bam
    # rm -f $TmpResultsDir/$BAMname-$TE-split-3+.bam
    # rm -f $TmpResultsDir/$BAMname-$TE-split-3-.bam
  

    #searching for overlapping clusters meeting the expected TSD size 
      echo "Searching overlaps and defining insertions..."
      
      $bedtoolsdir/intersectBed -a $TmpResultsDir/$BAMname-$TE-up-merge.bed -b $TmpResultsDir/$BAMname-$TE-down-merge.bed -wo | awk -v tsd=$TSD -v te=$TE -v rea=$READS '($6-$2)>10 && ($7-$3)>10 && $9>=tsd && ($4+$8)>=rea {print $1 "\t" $6 "\t" $3 "\t"te"\t" ($9-1)"\t"$4"\t"$8}' >> $OutputDir/$BAMname-$TE-insertion-sites.bed
      
      ###only if excluding inner pericentromeres and donnors TEs:
      #awk '$1!="" {print $0}' $SequencesDir/$TE.txt > $TmpResultsDir/$TE-donnor.txt
      #$bedtoolsdir/intersectBed -a $TmpResultsDir/$BAMname-$TE-up-merge.bed -b $TmpResultsDir/$BAMname-$TE-down-merge.bed -wo | awk -v tsd=$TSD -v te=$TE '($6-$2)>10 && ($7-$3)>10 && $9>=tsd {print $1 "\t" $6 "\t" $3 "\t"te"\t" ($9-1)"\t"$4"\t"$8}' | intersectBed -a stdin -b $TmpResultsDir/$TE-donnor.txt -v | intersectBed -a stdin -b $Cent -v >> $OutputDir/$BAMname-insertion-sites.bed

    INSERTIONS=$(wc -l $OutputDir/$BAMname-$TE-insertion-sites.bed | awk '{print $1}')
      
      echo "Split-read analyis done: $INSERTIONS putative insertions identified..."

    ###merging bam files and moving them to the output folder 
    $samtoolsDir/samtools merge -f $TmpResultsDir/$BAMname-$TE-split.bam $TmpResultsDir/$BAMname-$TE-up.bam $TmpResultsDir/$BAMname-$TE-down.bam 2>> $TmpDir/log.txt
    $samtoolsDir/samtools sort $TmpResultsDir/$BAMname-$TE-split.bam -o $OutputDir/$BAMname-$TE-split.bam 2>> $TmpDir/log.txt
    $samtoolsDir/samtools index $OutputDir/$BAMname-$TE-split.bam 2>> $TmpDir/log.txt

    #rm -r -f $TmpResultsDir

    ENDTIME=$(date +%s)
      echo "It takes $((ENDTIME-STARTTIME)) seconds to analyse $TE."

done

mv $TmpDir/log.txt $OutputDir
#rm -r -f $TmpDir
