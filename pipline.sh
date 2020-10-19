#!/bin/bash

# TODO: make it a program like bowtie2 that could be used in a command way
mkdir -p ./data
mkdir -p ./qcResult/
mkdir -p ./refData
mkdir -p ./interVar


# Variables for different datasets
# one-line method: samples=$(find ./data/ -type f -iname "*_1.fq.gz" -execdir sh -c 'printf "%s\n" "${0%_1.fq*}"' {} ';' | cut -d"/" -f2 | sort | uniq)
samples=$(cat ../refData/fqfiles | awk '{print $3}' | cut -c -6 | sort | uniq)
slender=$(cat ../refData/fqfiles | awk '$2=="Slender" {print $3}' | cut -c -6 | sort | uniq)
stumpy=$(cat ../refData/fqfiles | awk '$2=="Stumpy" {print $3}' | cut -c -6 | sort | uniq)

# Preload Functions
meanCal () {
    # Initialization
    sum=0;
    mean=0;

    # Calculate the sum
    for sl in $1;
    do
    slCount=$(cat ./interVar/$sl.sam.bam.sorted.txt | awk -v gene="$2" '$4==gene {print $7}')
    sum=$((sum + $slCount))
    done
    
    mean=$(awk -v suma="$sum" 'BEGIN{print suma/3}');
    
    # Assign the mean value
    if [[ "$1" == "$slender" ]]
    then 
    slenderMean=$mean
    elif [[ "$1" == "$stumpy" ]]
    then 
    stumpyMean=$mean
    else
    echo "Unknown Type"
    fi
}

geneMean (){
        # Initialisation
    slenderMean=0;
    stumpyMean=0;
    
    # Mean for slender samples and stumpy samples
    echo "Processing Gene $gene ...."
    meanCal "$slender" "$gene"
    meanCal "$stumpy" "$gene"

    # Generate the output
    echo -e "${gene}\t${slenderMean}\t${stumpyMean}" >> countStat.txt
}

# Main 
# Quick Check on the data by `fastqc`
zcat ../data/*fq.gz | fastqc --outdir=./qcResult/ stdin:fastqc_output

# Uncompress the reference genome
gunzip -c ../refData/Tb927_genome.fasta.gz > ./refData/Tb927_genome.fasta

# Build the index for the reference genome
bowtie2-build ../refData/Tb927_genome.fasta ./refData/Tb927_genome.fasta.index

for sample in $samples;
do 
    # Use bowtie2 to 
    echo -e "Processing $sample..."
    bowtie2 -p12 -x ./refData/Tb927_genome.fasta.index -1 "../data/${sample}_1.fq.gz" -2 "../data/${sample}_2.fq.gz" -S "./interVar/${sample}.sam"
done

# Use samtools to convert .sam to sorted.bam and .bai files
find ./interVar/*.sam | parallel "samtools view -bS {} -o {}.bam"
find ./interVar/*.bam | parallel "samtools sort {} -o {}.sorted"
find ./interVar/*.sorted | parallel "samtools index {}"

# BEDtools to count the gene-aligned sequences
find ./interVar/*.sorted | parallel "bedtools multicov -bams {} -bed ../refData/Tbbgenes.bed > {}.txt"

# Gene Summary
geneName=$(cat ./interVar/*.sam.bam.sorted.txt | awk '$5=="gene" {print $4}' | sort | uniq)

> countStat.txt

for gene in $geneName;
do
    ((i=i%12)); ((i++==0)) && wait
    geneMean &
done

echo "The analysis process has been DONE and the final result should be produced in the countStat.txt document in this directory."
