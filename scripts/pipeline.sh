#0.1 Get some small fastq WGS files  as example datasets (say at ~0.1 or ~0.01X coverage).

#0.2 Get the  GRCh37 genome & annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/

#0.3 STAR index the genome.
STAR --genomeSAindexNbases 12 \
     --runThreadN 12 \
     --runMode genomeGenerate \
     --genomeDir GRCh37 \
     --genomeFastaFiles GRCh37.p13.genome.fa \
     --sjdbGTFfile gencode.v19.annotation.gtf \
     --sjdbOverhang 99 \


#0.4 Align reads against genome using STAR . Also index bam files
nohup ./star_aligner.sh  >star_aligner.log&
samtools index *bam




#2.Download gnomad (for a list of germline variants to exclude)

temp=../temp/
outDir=../outDir/gnomad
#-------------------------

Rscript ./pipeline/genno01-get-gnomad.R $tmpDir $outDir
