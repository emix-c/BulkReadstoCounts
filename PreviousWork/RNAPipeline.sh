
#!/bin/bash
#fastq files pulled from SRA Toolkit -> pulled SRR15852443 
cd /Users/Emix/Desktop/RNASeq/SRAToolkit/sratoolkit.3.0.10-mac-x86_64/bin
prefetch SRR15852443
fasterq-dump SRR15852443 


#set working directory 
cd /Users/Emix/Desktop/RNASeq/


# STEP 1 : Quality Control via FASTQ + Trimming via Trimmetric 
fastqc SRAdata/SRR15852443_1.fastq -o SRAdata/ 
fastqc SRAdata/SRR15852443_2.fastq -o SRAdata/ 

echo "1st QC done!"  

java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -phred33 SRAdata/SRR15852443_1.fastq SRAdata/SRR15852443_2.fastq \
                SRAdata/SRR15852443_1_trimmed.fastq SRAdata/SRR15852443_1_untrim_fastq \
                SRAdata/SRR15852443_2_trimmed.fastq SRAdata/SRR15852443_2_untrim_fastq \
                ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10  TRAILING:10 

#want to trim 10 bases towards end of reads if below certain quality 
#convert quality scores to phred33 
#paired end reads 
#threads is # of processers i want trimmomatic to use 

echo "Trimmomatic finished running!"  

fastqc SRAdata/SRR15852443_1_trimmed.fastq -o SRAdata/ 
fastqc SRAdata/SRR15852443_2_trimmed.fastq -o SRAdata/ 

echo "2nd QC done!"  

# STEP 2 : Align via HISAT2 // need the reference genome indices 

hisat2 -x HISAT2/grch38/genome -1 SRAdata/SRR15852443_1_trimmed.fastq -2 SRAdata/SRR15852443_2_trimmed.fastq -S HISAT2/SRR15852443_trimmed.sam
#-q b/c reads are in fastq format 
#data is reverse stranded 
# provide indices after -x 
#genome indices are genome. (base name)

#redirect output to get a bam file instead of sam 
samtools view -b HISAT2/SRR15852443_trimmed.sam > HISAT2/SRR15852443_trimmed.bam
samtools sort HISAT2/SRR15852443_trimmed.bam -o HISAT2/SRR15852443_sorted.bam
echo "HISAT2 finished running!"  

mkdir quants_SRR15852443

# STEP 3 : Quantify aligned reads via FeatureCounts // need genome annotation files (gtf) via Ensembl 
#can tell which reads map to which part of the genome and can give the nmumber of reads that map to same region 
featureCounts -p -a /Users/Emix/Desktop/RNASeq/Homo_sapiens.GRCh38.111.gtf -o /Users/Emix/Desktop/RNASeq/quants_SRR15852443/featurecounts.txt /Users/Emix/Desktop/RNASeq/HISAT2/SRR15852443_sorted.bam

echo "featureCounts finished running!"  

