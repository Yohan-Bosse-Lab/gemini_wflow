genomedir='/mnt/sde/renseb01/Documents/rnaseq/data/reference_genome/GRCh37'
fastqdir='/mnt/sde/renseb01/Documents/gemini_wflow/fastq/'
annotation_gtf='/mnt/sde/renseb01/Documents/rnaseq/data/reference_genome/gencode.v19.annotation.gtf'

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
    --outFileNamePrefix /mnt/sde/renseb01/Documents/gemini_wflow/bams/$file 1>/mnt/sde/renseb01/Documents/gemini_wflow/bams/star.log 2>/mnt/sde/renseb01/Documents/gemini_wflow/bams/star.err
  done
