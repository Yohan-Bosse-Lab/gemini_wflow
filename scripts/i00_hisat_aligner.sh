#!/bin/bash

#dir
genomedir='/mnt/sde/renseb01/Documents/rnaseq/data/reference_genome/GRCh37'
genome_index='/mnt/sde/renseb01/Documents/rnaseq/data/reference_genome/genome_index/genome_index'

fastqdir='/home/renseb01/Documents/lord/raw_data/ctDNA_Gemini/GenomeQC/'
fastq_trim='/mnt/sde/renseb01/Documents/gemini_wflow/fastq_trim/'

annotation_gtf='/mnt/sde/renseb01/Documents/rnaseq/data/reference_genome/gencode.v19.annotation.gtf'
bamdir='/mnt/sde/renseb01/Documents/gemini_wflow/bams_hisat/'


#for loop
for R1 in $fastqdir/*_R1.fastq.gz
  do
    nameR1=${R1//'/home/renseb01/Documents/lord/raw_data/ctDNA_Gemini/GenomeQC/'}
    nameR2=${nameR1//'_R1'/'_R2'}

    trimgal_R1=${nameR1//'_R1.fastq.gz'/'_R1_val_1.fq.gz'}
    trimgal_R2=${nameR2//'_R2.fastq.gz'/'_R2_val_2.fq.gz'}
    cutadapt_R1=${nameR1//'_R1.fastq.gz'/'_cutadapt_R1.fq.gz'}
    cutadapt_R2=${nameR2//'_R2.fastq.gz'/'_cutadapt_R2.fq.gz'}

    output_sam=${nameR1//'_R1.fastq.gz'/'.sam'}
    output_bam=${output_sam//'sam'/'bam'}
    sorted_output_bam=${output_bam//'bam'/'sorted.bam'}

    echo "Running hisat: $file"
    echo $R1
    echo $nameR1
    echo $trimgal_R1
    echo $cutadapt_R1
    echo $output_sam
    echo $output_bam
    echo $sorted_output_bam

    #trim_galore
    trim_galore --cores 8 \
                -q 20 \
                -o $fastq_trim \
                --phred33 \
                --paired $fastqdir$nameR1 $fastqdir$nameR2  \
                --path_to_cutadapt cutadapt

    #remove the specific polyG sequences
    cutadapt -a GGGGGGGGGGGGGGG \
             -A GGGGGGGGGGGGGGG \
             --cores 8 \
             -o $fastq_trim$cutadapt_R1 \
             -p $fastq_trim$cutadapt_R2 \
             --discard-trimmed $fastq_trim$trimgal_R1 $fastq_trim$trimgal_R2

    #ref genome
    #hisat2-build /mnt/sde/renseb01/Documents/rnaseq/data/reference_genome/GRCh37.p13.genome.fa genome_index

    #align
    hisat2 -p 12 \
           -x $genome_index \
           -1 $fastq_trim$cutadapt_R1 \
           -2 $fastq_trim$cutadapt_R2 \
           -S $bamdir$output_sam

    #bam file
    samtools view -bS $bamdir$output_sam >$bamdir$output_bam

    #sorted bam
    samtools sort $bamdir$output_bam -o $bamdir$sorted_output_bam

    #index
    samtools index $bamdir$sorted_output_bam

    #a bit of clean-up
    rm  $bamdir$output_bam
    rm  $bamdir$output_sam
  done
