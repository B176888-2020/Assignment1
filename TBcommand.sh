#!/bin/bash

################ Arguments ################
# Input the arguments and make the shellscript a "command"
while getopts f:r:b:o: option
do
case "${option}"
in
f) FQFILE=${OPTARG};;
r) REFGENOME=${OPTARG};;
b) BEDFILE=${OPTARG};;
o) OUTPUT=${OPTARG};;
esac
done

################ Preparation ################
# Make directories to store intermediate documents
mkdir -p ${OUTPUT}/qcResult/; mkdir -p ${OUTPUT}/refData; mkdir -p ${OUTPUT}/interVar;

# Variables stored the sample names
samples=$(cat ${FQFILE}/fqfiles | awk '{print $3}' | cut -c -6 | sort | uniq)
slender=$(cat ${FQFILE}/fqfiles | awk '$2=="Slender" {print $3}' | cut -c -6 | sort | uniq)
stumpy=$(cat ${FQFILE}/fqfiles | awk '$2=="Stumpy" {print $3}' | cut -c -6 | sort | uniq)

# Preload Functions
# Function used to calcualte the mean of two types of dataset
meanCal () {
    # Initialization
    sum=0; mean=0;

    # Calculate the mean
    for sl in $1;
    do
    slCount=$(cat ${OUTPUT}/interVar/$sl.sam.bam.sorted.txt | awk -v gene="$2" '$4==gene {print $7}')
    sum=$((sum + $slCount))
    done
    mean=$(awk -v suma="$sum" 'BEGIN{print suma/3}');
    
    # Assign the mean value to proper variable
    if [[ "$1" == "$slender" ]]
    then 
    slenderMean=$mean
    elif [[ "$1" == "$stumpy" ]]
    then 
    stumpyMean=$mean
    else
    echo "Warning: Unknown Type"
    fi
}

# Function used to generate the summary of count data
geneMean (){
    # Initialisation
    slenderMean=0;
    stumpyMean=0;
    
    # Mean for slender samples and stumpy samples
    echo "Processing Gene $gene ...."
    meanCal "$slender" "$gene" # Actually this function could get two values at the same time with additional tag column to distinct `slender` and `stumpy`
    meanCal "$stumpy" "$gene" # However, it would require intermediate documents that may not be required in the assignment. To keep it simple, just keep it. 

    # Generate the output
    echo -e "${gene}\t${slenderMean}\t${stumpyMean}" >> countStat.txt
}


################ Main Pipeline Process ################
# Quality Control with `FastQC`
## The quality for all the data
zcat ${FQFILE}/*fq.gz | fastqc --outdir=${OUTPUT}/qcResult/ stdin:allSamples
## The quality for each fastq.gz data
fastqc -t 6 -noextract --outdir=${OUTPUT}/qcResult/ ${FQFILE}/*fq.gz 

# Uncompress the reference genome
gunzip -c ${REFGENOME}/Tb927_genome.fasta.gz > ${OUTPUT}/refData/Tb927_genome.fasta

# Build the index for the reference genome
bowtie2-build ${OUTPUT}/refData/Tb927_genome.fasta ${OUTPUT}/refData/Tb927_genome.fasta.index

# Align the data to the reference genome
for sample in $samples;
do 
    # Use bowtie2 to 
    echo -e "Processing $sample..."
    bowtie2 -p12 -x ${OUTPUT}/refData/Tb927_genome.fasta.index -1 "${FQFILE}/${sample}_1.fq.gz" -2 "${FQFILE}/${sample}_2.fq.gz" -S "${OUTPUT}/interVar/${sample}.sam"
done

# Use samtools to convert .sam to sorted.bam and .bai files
find ${OUTPUT}/interVar/*.sam | parallel "samtools view -bS {} -o {}.bam"
find ${OUTPUT}/interVar/*.bam | parallel "samtools sort {} -o {}.sorted"
find ${OUTPUT}/interVar/*.sorted | parallel "samtools index {}"

# BEDtools to count the gene-aligned sequences
find ${OUTPUT}/interVar/*.sorted | parallel "bedtools multicov -bams {} -bed ${BEDFILE}/Tbbgenes.bed > {}.txt"

# Gene count data Summary
geneName=$(cat ${OUTPUT}/interVar/*.sam.bam.sorted.txt | awk '$5=="gene" {print $4}' | sort | uniq)
> countStat.txt
for gene in $geneName;
do
    ((i=i%12)); ((i++==0)) && wait
    geneMean &
done

echo "The analysis process has been DONE and the final result should be produced in the countStat.txt document in this directory."
