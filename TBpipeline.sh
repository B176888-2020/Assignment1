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
slenderN=$(cat ${FQFILE}/fqfiles | awk '$2=="Slender" {slN+=1}; END{print slN}')
stumpy=$(cat ${FQFILE}/fqfiles | awk '$2=="Stumpy" {print $3}' | cut -c -6 | sort | uniq)
stumpyN=$(cat ${FQFILE}/fqfiles | awk '$2=="Stumpy" {stN+=1}; END{print stN}')

################ Main Pipeline Process ################
# Quality Control with `FastQC`
zcat ${FQFILE}/*fq.gz | fastqc --outdir=${OUTPUT}/qcResult/ stdin:allSamples # The quality check for all the data
fastqc -t 6 -noextract --outdir=${OUTPUT}/qcResult/ ${FQFILE}/*fq.gz # The quality check for each fastq.gz data

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
bedtools multicov -bams $(eval echo ${OUTPUT}/interVar/{$(echo ${slender} ${stumpy} | tr " " ",")}.sam.bam.sorted) -bed ${BEDFILE}/Tbbgenes.bed > ${OUTPUT}/interVar/counts.txt
> ${OUTPUT}/countStat.txt
cat $OUTPUT/interVar/counts.txt | awk -v slNawk="$slenderN" -v stNawk="$stumpyN" '
    BEGIN{       
        FS="\t"; OFS="\t";
    }
    {
        x=0; w=0
        for (i=7;i<=6+slNawk;i++){
            x=x+$i;
        }
        for (j=7+slNawk;j<=6+slNawk+stNawk;j++){
            w=w+$j;
        }
        print $4,(x/slNawk),(w/stNawk)
    }
    ' >> ${OUTPUT}/countStat.txt

echo "The analysis process has been DONE and the final result should be produced in the countStat.txt document in this directory."
