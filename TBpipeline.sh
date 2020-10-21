#!/bin/bash

################ Arguments ################
# Input the arguments and make the shell script a "command"
while getopts f:r:b:o:y option
do
case "${option}"
in
f) FQDIR=${OPTARG};;
r) REFGENOME=${OPTARG};;
b) BEDFILE=${OPTARG};;
o) OUTPUT=${OPTARG};;
y) QCPASS='pass';;
esac
done

################ Preparation ################
# Make directories to store intermediate documents
mkdir -p ${OUTPUT}qcResult/; mkdir -p ${OUTPUT}refData; mkdir -p ${OUTPUT}interVar;

# Variables stored the names and numbers of samples
samples=$(cat ${FQDIR}fqfiles | awk '{print $3}' | cut -c -6 | sort | uniq)
slender=$(cat ${FQDIR}fqfiles | awk '$2=="Slender" {print $3}' | cut -c -6 | sort | uniq)
slenderN=$(cat ${FQDIR}fqfiles | awk '$2=="Slender" {slN+=1}; END{print slN}')
stumpy=$(cat ${FQDIR}fqfiles | awk '$2=="Stumpy" {print $3}' | cut -c -6 | sort | uniq)
stumpyN=$(cat ${FQDIR}fqfiles | awk '$2=="Stumpy" {stN+=1}; END{print stN}')

################ Main Process ################
# Quality Control with `FastQC`
zcat ${FQDIR}*fq.gz | fastqc --extract --outdir=${OUTPUT}qcResult/ stdin:allSamples # The quality check for all the data
fastqc -t 6 --extract --outdir=${OUTPUT}qcResult/ ${FQDIR}*fq.gz # The quality check for each fastq.gz data

# Assess the number and quality of the raw sequence data
> ${OUTPUT}qcResult/qcResultSummary.txt
echo -e "\nFastQC Summary: "
for qcdata in $(eval echo ${OUTPUT}qcResult/*fastqc/summary.txt);
do
passN=$(cat $qcdata | grep -c "PASS")
warnN=$(cat $qcdata | grep -c "WARN")
failN=$(cat $qcdata | grep -c "FAIL")
echo -e "Data: $qcdata  PASS: $passN WARN: $warnN FAIL: $failN" | tee -a ${OUTPUT}qcResult/qcResultSummary.txt
done

# After getting the FastQC summary, decide whether to go through downstream processes
if [ "$QCPASS" = 'pass' ]; then
key=''
else
read -n1 -rsp $'\nIf everything is OK, please press SPACE to continue...or press ANY KEY to exit if you need to check input datasets or detailed outputs from FastQC.\n' key
fi

if [ "$key" = '' ]; then

# Uncompress the reference genome and build the index for it
gunzip -c ${REFGENOME}Tb927_genome.fasta.gz > ${OUTPUT}refData/Tb927_genome.fasta
bowtie2-build --threads 12 ${OUTPUT}refData/Tb927_genome.fasta ${OUTPUT}refData/Tb927_genome.fasta.index

# Align the sequence data to the reference genome
echo "Aligning the sequences..."
for sample in $samples;
do 
    echo -e "Processing $sample..."
    bowtie2 --sensitive -p12 -x ${OUTPUT}refData/Tb927_genome.fasta.index -1 "${FQDIR}${sample}_1.fq.gz" -2 "${FQDIR}${sample}_2.fq.gz" -S "${OUTPUT}interVar/${sample}.sam"
done

# Use samtools to convert .sam to sorted.bam and produce .bai files
echo "Transforming: SAM > BAM BAI..."
find ${OUTPUT}interVar/*.sam | parallel "samtools view -bS {} -o {}.bam"
find ${OUTPUT}interVar/*.bam | parallel "samtools sort {} -o {}.sorted"
find ${OUTPUT}interVar/*.sorted | parallel "samtools index {}"

# Use bedtools to count the gene-aligned sequences
echo "Generating counts data and statistical mean summary..."
bedtools multicov -bams $(eval echo ${OUTPUT}interVar/{$(echo ${slender} ${stumpy} | tr " " ",")}.sam.bam.sorted) -bed ${BEDFILE}Tbbgenes.bed > ${OUTPUT}interVar/counts.txt
# Generate the statistic mean of the count data
> ${OUTPUT}countStat.txt
cat ${OUTPUT}interVar/counts.txt | awk -v slNawk="$slenderN" -v stNawk="$stumpyN" '
    BEGIN{       
        FS="\t"; OFS="\t";
    }
    {
        slSum=0; stSum=0
        for (i=7;i<=6+slNawk;i++){slSum=slSum+$i;}
        for (j=7+slNawk;j<=6+slNawk+stNawk;j++){stSum=stSum+$j;}
        print $4,(slSum/slNawk),(stSum/stNawk)
    }
    ' >> ${OUTPUT}countStat.txt

echo "The analysis process has been DONE and the statistic mean summary of gene count data would be produced in the countStat.txt document in ${OUTPUT} directory."

else
    echo -e "\nReminder: The information above (qcResultSummary.txt) and other FastQC outputs are stored in the ${OUTPUT}qcResults/ directory."
    exit
fi
