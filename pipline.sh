#!/bin/bash

# TODO: make it a program like bowtie2 that could be used in a command way

# Quick Check on the data by `fastqc`
zcat ./data/*fq.gz | fastqc stdin:fastqc_output

# Uncompress the reference genome
gunzip -c ./refData/Tb927_genome.fasta.gz > ./refData/Tb927_genome.fasta


# Build the index for the reference genome
bowtie2-build ./refData/Tb927_genome.fasta ./refData/Tb927_genome.fasta.index

# Loop Method
samples=$(find ./data/ -type f -iname "*_1.fq.gz" -execdir sh -c 'printf "%s\n" "${0%_1.fq*}"' {} ';' | cut -d"/" -f2 | sort | uniq)

for sample in $samples;
do 
# Use bowtie2 to 
echo -e "Processing $sample..."
bowtie2 -p12 -x ./refData/Tb927_genome.fasta.index -1 "./data/${sample}_1.fq.gz" -2 "./data/${sample}_2.fq.gz" -S "./InterVar/${sample}.sam"
done


# Use samtools to convert .sam to .bam
find ./InterVar/*.sam | parallel "samtools view -bS {} -o {}.bam"
find ./InterVar/*.bam | parallel "samtools sort {} -o {}.sorted"
find ./InterVar/*.sorted | parallel "samtools index {}"

# TODO: all .bams should be calculate and to do the following method
# BEDtools to count the gene-aligned sequences
find ./InterVar/*.sorted | parallel "bedtools multicov -bams {} -bed ./refData/Tbbgenes.bed > {}.txt"

# TODO: Final Summary




# Alternatvie Trials
###### loop or Expansion; Possible var: counter for the prefix and ######

# Parallel Method: `parallel` > two input problem ; 'bowtie2' > only one output for multisample, strange
## pair1=$(ls ./data/*1.fq.gz | sort | awk '{printf "%s,",$0;}');
## pair2=$(ls ./data/*2.fq.gz | sort | awk '{printf "%s,",$0;}');
## names=$(find ./data/ -type f -iname "*_1.fq.gz" -execdir sh -c 'printf "%s\n" "${0%_1.fq*}"' {} ';' | cut -d"/" -f2 | sort | awk '{printf "%s,",$0;}');
## bowtie2 -p12 -x ./refData/Tb927_genome.fasta.index -1 $pair1 -2 $pair2 -S $names.sam

# GNU Parallel should do a much beter job but I dont know how to add specific two inputs at the same time with particular arguments for the next pipe
# find ./data/*1.fq.gz ./data/*2.fq.gz | sort | parallel "bowtie2 -x ./refData/Tb927_genome.fasta.index -1 {$1} -2 {$2} -S {$1}.sam;"



