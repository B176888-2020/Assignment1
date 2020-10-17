#!/bin/bash

# TODO: make it a program like bowtie2 that could be used in a command way
# TODO: make the directory more tidy during the processing

# Quick Check on the data by `fastqc`
zcat ./data/*fq.gz | fastqc stdin:fastqc_output;

# Uncompress the reference genome
gunzip -k Tb927_genome.fasta.gz;


# Use bowtiew2 to align 
bowtie2-build Tb927_genome.fasta Tb927_genome.fasta.index;

# TODO: auto-get the number of the pairs and assign to `num_pairs` 

# TODO: use `num_pairs` to loop 
###### loop or Expansion; Possible var: counter for the prefix and ######
bowtie2 -x Tb927_genome.fasta.index -1 ./data/216_L8_1.fq.gz -2 ./data/216_L8_2.fq.gz -S 216_L8.sam;
# Use samtools to convert .sam to .bam
samtools view -bS 216_L8.sam > 216_L8.bam

# TODO: Generate Counts data
bedtools 
###### loop or Expansion ######





