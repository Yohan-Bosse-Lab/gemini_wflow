
#a mini script to collect read statistics (idxstats) on files.


#paths
raw_files2 = list.files('/home/renseb01/Documents/lord/raw_data/ctDNA_Gemini/GenomeQC/',pattern = '*R1.fastq.gz',full.names = T)
raw_files1 = list.files('/home/renseb01/Documents/lord/raw_data/ctDNA_Gemini/IUCPQ/',pattern = '*R1.fastq.gz',full.names = T)
raw_files = c(raw_files1,raw_files2)

processed_files2 = list.files('/home/renseb01/Documents/lord/processed_data/ctDNA_Gemini/GenomeQC/',pattern = '*R1.fq.gz',full.names = T)
processed_files1 = list.files('/home/renseb01/Documents/lord/processed_data/ctDNA_Gemini/IUCPQ/',pattern = '*R1.fq.gz',full.names = T)
processed_files = c(processed_files1,processed_files2)

aligned_files =  list.files('/mnt/sde/renseb01/Documents/gemini_wflow/bams_hisat/',pattern = 'sorted.bam$',full.names = T)

#results_df
results = data.frame(raw_files = raw_files, processed_files = processed_files, aligned_files = aligned_files, raw = 0, processed = 0, aligned = 0, raw_p = 0, processed_p = 0, aligned_p = 0)

#for loop
for(i in 1:31){
  #raw
  cmd1 = paste0('samtools idxstats ',results$raw_files[i],' >idxstats.raw')
  system(cmd1)
  idxstats.raw = read.table('idxstats.raw')
  results$raw[i] = round(idxstats.raw[1,4]/1000000,2)

  #processed
  cmd2 = paste0('samtools idxstats ',results$processed_files[i],' >idxstats.processed')
  system(cmd2)
  idxstats.processed = read.table('idxstats.processed')
  results$processed[i] = round(idxstats.processed[1,4]/1000000,2)

  #aligned
  cmd3 = paste0('samtools idxstats ',results$aligned_files[i],' >idxstats.aligned')
  system(cmd3)
  idxstats.aligned = read.table('idxstats.aligned')
  results$aligned[i] = round(sum(idxstats.aligned[,3])/1000000,2)

  print(paste0('Done ',i,', The time is: ',Sys.time()))
}

system('rm idxstats.raw idxstats.processed idxstats.aligned')

#
results$aligned_p = round(results$aligned / results$processed * 50,2)
results$raw_p = round(results$processed / results$raw * 100,2)

write.table(results,'../outDir/alignment_stats.csv')
