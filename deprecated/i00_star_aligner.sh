#!/bin/bash

#dir
genomedir='/mnt/sde/renseb01/Documents/rnaseq/data/reference_genome/GRCh37'
fastqdir='/mnt/sde/renseb01/Documents/gemini_wflow/fastq/'
annotation_gtf='/mnt/sde/renseb01/Documents/rnaseq/data/reference_genome/gencode.v19.annotation.gtf'
bamdir='/mnt/sde/renseb01/Documents/gemini_wflow/bams/'


#for loop
for R1 in $fastqdir/*R1.fastq.gz
  do
    echo $R1
    R2=${R1//'_R1'/'_R2'}
    echo $R2
    temp=${R1//'/mnt/sde/renseb01/Documents/gemini_wflow/fastq//'}
    file=${temp//'R1.fastq.gz'/}
    echo $temp
    echo $file
    echo "Running STAR: $file"

    #STAR cmd
    STAR --runMode alignReads \
    --genomeDir $genomedir \
    --readFilesCommand zcat \
    --limitBAMsortRAM 100000000000 \
    --readFilesIn $R1 $R2 \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbOverhang 99 \
    --outFilterMultimapNmax 100 \
    --outReadsUnmapped None \
    --quantMode GeneCounts \
   --sjdbGTFfile $annotation_gtf \
    --runThreadN 12 \
    --outFileNamePrefix $bamdir$file 1>$bamdir'star.log' 2>$bamdir'star.err'
  done



for R1 in $fastqdir/*R1.fastq.gz
  do
    temp=${R1//'/mnt/sde/renseb01/Documents/gemini_wflow/fastq//'}
    file=${temp//'_R1.fastq.gz'/}


    #cleanup
    mv $bamdir$file'_Aligned.sortedByCoord.out.bam' $bamdir$file'.bam'

    rm -r $bamdir*STARgenome*
    rm $bamdir*og*
    rm $bamdir*tab
  done
