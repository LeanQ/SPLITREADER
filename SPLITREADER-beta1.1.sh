#!/bin/bash 
###SPLITREADER beta1.1 
###Leandro Quadrana 

#$1 -> basename of the bam file (i.e. if the file is "sequencing.bam", you should enter only "sequencing")

#TE-> should be indicated in the first column of the TE-information.txt file located in SequencesDir
#TSD -> should be indicated in the second column of the TE-information.txt file located in SequencesDir


#############################################################
#IN THE FOLLOWING VARIABLE YOU CAN EXPLICITE THE MINIMUM READ LEGTH. By default this is 100nt
LENGTH=100

#### If not specified, the program calculate teh longer read


#############################################################


#############################################################
#IN THE FOLLOWING VARIABLE YOU CAN EXPLICITE THE MINIMUM NUMBER OF SPLIT-READS IN EACH EXTRIMITY OF THE TE. By default is 5 reads
READS=5
#### If not specified, the program calculate it. To this end, the program calculates the minimum number of reads as 3 standard deviation under the mean whole genome coverage
####This value should be at least 3, if not, it is forced to be 3

#############################################################

#############################################################
#IN THE FOLLOWING VARIABLE YOU CAN EXPLICITE THE NUMBER OF THREADS YOU WANT TO USE FOR MAPPING. By default is 2 threads
CORES=2
#############################################################



# Path to configure (by default it is the current directory)

InputDir=./
OutputDir=./

#list of TEs to analize
listDir=/.

#Bowtie2 index files for each TE sequence
SequencesDir=/.

#Bowtie2 index file for reference genome
GenomeIndexFile=./TAIR10index

#Coordinates of inner pericentromeres 
Cent=./centromeres.bed


# PID of this batch
IDPID=$$
# Temporary directory
TmpDir=./QD-$IDPID
mkdir -p $TmpDir


###If required, stimating whole genome coverage and TE insertion tresholds
  if [ -z "$READS" ]
     then
     samtools depth $InputDir/$1.bam > $TmpDir/coverage.temp
      co=`awk '{sum+=$3} END { print int(sum/NR)}' $TmpDir/coverage.temp`
      sd=`awk '{sum+=$3; sumsq+=$3*$3} END { print int(sqrt(sumsq/NR - (sum/NR)**2))}' $TmpDir/coverage.temp`
      cov=$((co-3*sd))
      maxcov=$((co+3*sd))
      rm -f $TmpDir/coverage.temp
      
       if [ $cov -gt 6 ]
	then
	  READS=$cov
	  "Minimum number of Reads in split-reads cluster: $READS [Estimated] "
	else
	  READS=6
	  "Minimum number of Reads in split-reads cluster:  $READS [Default]"
	fi
      else
      echo "Minimum number of Reads in split-reads cluster:  $READS [User defined]"
  fi

  
#Extracting unmapped reads
echo "Extracting unmapped reads from $1"

samtools view -f 4 -u $InputDir/$1.bam | bamToFastq -bam stdin -fq1 $TmpDir/$1.1.fastq -fq2 $TmpDir/$1.2.fastq

end=`wc -l $listDir/TE-information.txt | awk '{print $1}'`


#Starting the SPLITREADER pipeline for each TE in the TE-indormation.txt file
for ((l=1; $l<=$end; l=$l+1)); do
    STARTTIME=$(date +%s)

    IDPID2=$PPID
    TmpResultsDir=$TmpDir/results-$IDPID2
    mkdir -p $TmpResultsDir

    TE=`sed -n "${l}p" $listDir/TE-information.txt | awk '{print $1}'`
    TSD=`sed -n "${l}p" $listDir/TE-information.txt | awk '{print $2}'`
   
   echo "##### RUNNING SPLIT-READ ANALYSIS ON $TE ######"    
   echo ""
   
   
   # Selecting split-reads by mapping the unmapped reads over TE extremities
   
   echo "Selecting split-reads"
   
    bowtie2 -x $SequencesDir/$TE -U $TmpDir/$1.1.fastq,$TmpDir/$1.2.fastq -S $TmpResultsDir/$1-$TE.sam --local --very-sensitive --threads $CORES 
   
    T=$(awk '$1~/@/ {x++} END {print x}' $TmpResultsDir/$1-$TE.sam)

    samtools view -h -S $TmpResultsDir/$1-$TE.sam | head -$T > $TmpResultsDir/$1-$TE-split.sam 
    samtools view -F 4 -S $TmpResultsDir/$1-$TE.sam | awk '$6~S {print $0}' >> $TmpResultsDir/$1-$TE-split.sam
    rm -f $TmpResultsDir/$1-$TE.sam

    java -jar /groups/a2e/Leandro/split-analysis/SamToFastq.jar INPUT=$TmpResultsDir/$1-$TE-split.sam FASTQ=$TmpResultsDir/$1-$TE-split.fastq
    rm -f $TmpResultsDir/$1-$TE-split.sam

    ###Stimating max read size (If necessary)
    
    if [ -z "$LENGTH" ]
     then
      LENGTH=`awk 'NR%4 == 2 {print length($0)}' $TmpResultsDir/$1-$TE-split.fastq | sort | tail -1 `  
      length=$(($LENGTH-20))
      echo "Maximum Read legth: $LENGTH [Stimated] "
      else
      length=$(($LENGTH-20))
      echo "Maximum Read legth: $LENGTH [User defined] "
     fi
    
         
    ###Recursive split-reads mapping
    # step 1 for 3' read extremity: begining the loop.
    
      
   echo "Mapping 5' split-reads on reference genome"
   echo -ne "Progresssion: ["
    bowtie2 -x $GenomeIndexFile -U $TmpResultsDir/$1-$TE-split.fastq -S $TmpResultsDir/$1-$TE-splitjunction-3-20.sam --un $TmpResultsDir/$1-$TE-split-3-20 --quiet -3 20 --mp 13 --rdg 8,5 --rfg 8,5 --very-sensitive --threads $CORES
    samtools view -q 5 -Sbu $TmpResultsDir/$1-$TE-splitjunction-3-20.sam | samtools sort - $TmpResultsDir/$1-$TE-splitjunction-3-20
    rm -f $TmpResultsDir/$1-$TE-splitjunction-3-20.sam
    
    
    #  step 2 for 3' read extremity: loop for the mapping from read-length LENGTH-20 up to 20
    
    for ((i=21; $i<=$length; i=$i+1)); do
    
      previous=$(($i-1));

      bowtie2 -x $GenomeIndexFile -U $TmpResultsDir/$1-$TE-split-3-$previous -S $TmpResultsDir/$1-$TE-splitjunction-3-$i.sam --un $TmpResultsDir/$1-$TE-split-3-$i --quiet -3 $i --mp 13 --rdg 8,5 --rfg 8,5 --very-sensitive
      samtools view -q 5 -Sbu $TmpResultsDir/$1-$TE-splitjunction-3-$i.sam | samtools sort - $TmpResultsDir/$1-$TE-splitjunction-3-$i 
      rm -f $TmpResultsDir/$1-$TE-splitjunction-3-$i.sam
      rm -f $TmpResultsDir/$1-$TE-split-3-$previous
    done
    
    
    echo "Mapping 3' split-reads on reference genome"
    echo -ne "Progresssion: ["
    
    
    # step 1 for 5' read extremity: begining the loop.
 
    bowtie2 -x $GenomeIndexFile -U $TmpResultsDir/$1-$TE-split.fastq -S $TmpResultsDir/$1-$TE-splitjunction-5-20.sam --un $TmpResultsDir/$1-$TE-split-5-20 --quiet -5 20 --mp 13 --rdg 8,5 --rfg 8,5 --very-sensitive 
    samtools view -q 5 -Sbu $TmpResultsDir/$1-$TE-splitjunction-5-20.sam | samtools sort - $TmpResultsDir/$1-$TE-splitjunction-5-20
    rm -f $TmpResultsDir/$1-$TE-splitjunction-5-20.sam
    
    
    #  step 2 for 5' read extremity: loop for the mapping from read-length LENGTH-20 up to 20
  
  
    for ((i=21; $i<=$length; i=$i+1)); do
    
      previous=$(($i-1));

      bowtie2 -x $GenomeIndexFile -U $TmpResultsDir/$1-$TE-split-5-$previous -S $TmpResultsDir/$1-$TE-splitjunction-5-$i.sam --un $TmpResultsDir/$1-$TE-split-5-$i --quiet -5 $i --mp 13 --rdg 8,5 --rfg 8,5 --very-sensitive --threads $CORES
      samtools view -q 5 -Sbu $TmpResultsDir/$1-$TE-splitjunction-5-$i.sam | samtools sort - $TmpResultsDir/$1-$TE-splitjunction-5-$i
      rm -f $TmpResultsDir/$1-$TE-splitjunction-5-$i.sam
      rm -f $TmpResultsDir/$1-$TE-split-5-$previous
    done
    
    
    
    # Post-treatment:
    
    # merge all the recursive mappings from either the 3' and 5'
        
    rm -f $TmpResultsDir/$1-$TE-split.fastq

    # Merging and sorting bam files
    samtools merge -u $TmpResultsDir/$1-$TE-split-5.bam $TmpResultsDir/$1-$TE-splitjunction-5-*.bam
    samtools merge -u $TmpResultsDir/$1-$TE-split-3.bam $TmpResultsDir/$1-$TE-splitjunction-3-*.bam
    samtools sort $TmpResultsDir/$1-$TE-split-5.bam $TmpResultsDir/$1-$TE-split-5
    samtools sort $TmpResultsDir/$1-$TE-split-3.bam $TmpResultsDir/$1-$TE-split-3
    samtools index $TmpResultsDir/$1-$TE-split-5.bam
    samtools index $TmpResultsDir/$1-$TE-split-3.bam 
    rm -f $TmpResultsDir/$1-$TE-splitjunction-[35]-*.bam
    
    # merge reads that were cliped at the 3' and mapped in the + strand with those clipped at the 5' and mapped on the - strand
    # merge reads that were cliped at the 3' and mapped in the - strand with those clipped at the 5' and mapped on the + strand
    
      echo "Searching for reads clusters..."
    samtools view -F 16 -bh $TmpResultsDir/$1-$TE-split-5.bam > $TmpResultsDir/$1-$TE-split-5+.bam
    samtools view -f 16 -bh $TmpResultsDir/$1-$TE-split-5.bam > $TmpResultsDir/$1-$TE-split-5-.bam
    samtools view -F 16 -bh $TmpResultsDir/$1-$TE-split-3.bam > $TmpResultsDir/$1-$TE-split-3+.bam
    samtools view -f 16 -bh $TmpResultsDir/$1-$TE-split-3.bam > $TmpResultsDir/$1-$TE-split-3-.bam

    
    #Merge the 5' and 3' clusters to create the downstream and upstream cluster

    samtools merge $TmpResultsDir/$1-$TE-up.bam $TmpResultsDir/$1-$TE-split-5-.bam $TmpResultsDir/$1-$TE-split-3+.bam
    samtools merge $TmpResultsDir/$1-$TE-down.bam $TmpResultsDir/$1-$TE-split-5+.bam $TmpResultsDir/$1-$TE-split-3-.bam

    #Calculate de coverage over mapped regions - filter regions according to minimum and maximum read-depth
    samtools depth $TmpResultsDir/$1-$TE-up.bam | awk -v m=$READS M=$maxcov '$3>=(m/2) $$ $3<(M/2) {print $1 "\t" $2 "\t"$2"\t"$3}' > $TmpResultsDir/$1-$TE-up.bed
    samtools depth $TmpResultsDir/$1-$TE-down.bam | awk -v m=$READS M=$maxcov '$3>=(m/2) $$ $3<(M/2) {print $1 "\t" $2 "\t"$2"\t"$3}' > $TmpResultsDir/$1-$TE-down.bed

    #merge cluster of covered regions - filter out clusters that are longer than read-length
    sort -k 1,1 -k2,2n $TmpResultsDir/$1-$TE-up.bed | mergeBed -i stdin | awk -v L=$length '$3-$2<L {print $0}' > $TmpResultsDir/$1-$TE-up-merge-temp.bed
    sort -k 1,1 -k2,2n $TmpResultsDir/$1-$TE-down.bed | mergeBed -i stdin | awk -v L=$length '$3-$2<L {print $0}' > $TmpResultsDir/$1-$TE-down-merge-temp.bed
    
    #recover max coverage information over clusters 
    intersectBed -a $TmpResultsDir/$1-$TE-up-merge-temp.bed -b $TmpResultsDir/$1-$TE-up.bed -wo | sort -k 1,1 -k2,2n | bedtools groupBy -c 7 -o max > $TmpResultsDir/$1-$TE-up-merge.bed
    intersectBed -a $TmpResultsDir/$1-$TE-down-merge-temp.bed -b $TmpResultsDir/$1-$TE-down.bed -wo | sort -k 1,1 -k2,2n | bedtools groupBy -c 7 -o max > $TmpResultsDir/$1-$TE-down-merge.bed
    

    rm -f $TmpResultsDir/$1-$TE-split-5+.bam
    rm -f $TmpResultsDir/$1-$TE-split-5-.bam
    rm -f $TmpResultsDir/$1-$TE-split-3+.bam
    rm -f $TmpResultsDir/$1-$TE-split-3-.bam
    rm -f $TmpResultsDir/$1-$TE-up.bam
    rm -f $TmpResultsDir/$1-$TE-down.bam


    #searching for overlapping clusters meeting the expected TSD size 
      echo "Searching overlaps and defining insertions..."
      awk '$1!="" {print $0}' $SequencesDir/$TE.txt > $TmpResultsDir/$TE-donnor.txt
      intersectBed -a $TmpResultsDir/$1-$TE-up-merge.bed -b $TmpResultsDir/$1-$TE-down-merge.bed -wo | awk -v tsd=$TSD -v te=$TE '($6-$2)>10 && ($7-$3)>10 && $9>=tsd {print $1 "\t" $6 "\t" $3 "\t"te"\t" ($9-1)"\t"$4"\t"$8}' | intersectBed -a stdin -b $TmpResultsDir/$TE-donnor.txt -v | intersectBed -a stdin -b $Cent >> $OutputDir/$1-insertion-sites.bed

    INSERTIONS=$(wc -l $OutputDir/$1-$TE-insertion-sites.bed | awk '{print $1}')
      
      echo "Split-read analyis done: $INSERTIONS putative insertions identified..."

    cp $TmpResultsDir/$1-$TE-split-[35].bam* $OutputDir/ ###moving splitread mapped files to the output folder 

    rm -f $TmpResultsDir/$1-$TE-up-merge.bed
    rm -f $TmpResultsDir/$1-$TE-down-merge.bed
    rm -f $TmpResultsDir/$TE-donnor.txt
    rm -r -f $TmpResultsDir

    ENDTIME=$(date +%s)
      echo "It takes $((ENDTIME-STARTTIME)) seconds to analyse $TE."

done

rm -r -f $TmpDir
