#--------
# Input
#-------------------------
tmpDir=../temp
outDir=../outDir/gnomad
#-------------------------

mkdir -p $tmpDir $outDir

# They now appear to be hosted here (the md5sum values match those from the files above):
#wget http://hgdownload-euro.soe.ucsc.edu/gbdb/hg38/gnomAD/vcf/gnomad.genomes.r3.0.sites.vcf.gz -O $tmpDir/gnomad.genomes.r3.0.sites.vcf.gz
#wget http://hgdownload-euro.soe.ucsc.edu/gbdb/hg38/gnomAD/vcf/gnomad.genomes.r3.0.sites.vcf.gz.tbi -O $tmpDir/gnomad.genomes.r3.0.sites.vcf.gz.tbi

Rscript ./i01-get-gnomad.R $tmpDir $outDir
