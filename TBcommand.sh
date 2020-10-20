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
# Function used to generate the summary of count data
geneMean (){
    # Initialisation
    slenderMean=0;
    stumpyMean=0;
    
    # Mean for slender samples and stumpy samples
    echo "Processing Gene $gene ...."
    slenderMean=$(cat ${OUTPUT}/interVar/slender.txt | awk -v gene="$1" '$4==gene {sum=0; for( i = 7; i <= NF; i++ ){sum+=$i};} END{print sum/(NF-6)}')
    stumpyMean=$(cat ${OUTPUT}/interVar/stumpy.txt | awk -v gene="$1" '$4==gene {sum=0; for( i = 7; i <= NF; i++ ){sum+=$i};} END{print sum/(NF-6)}')

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
bedtools multicov -bams $(eval echo ${OUTPUT}/interVar/{$(echo ${slender} | tr " " ",")}.sam.bam.sorted) -bed ${BEDFILE}/Tbbgenes.bed > ${OUTPUT}/interVar/slender.txt
bedtools multicov -bams $(eval echo ${OUTPUT}/interVar/{$(echo ${stumpy} | tr " " ",")}.sam.bam.sorted) -bed ${BEDFILE}/Tbbgenes.bed > ${OUTPUT}/interVar/stumpy.txt

# Gene count data Summary
geneName=$(cat slender.txt stumpy.txt | awk '$5=="gene" {print $4}' | sort | uniq)
> countStat.txt
for gene in $geneName;
do
    ((i=i%12)); ((i++==0)) && wait
    geneMean &
done

echo "The analysis process has been DONE and the final result should be produced in the countStat.txt document in this directory."
