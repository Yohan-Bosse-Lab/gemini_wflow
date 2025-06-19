
#a mini script to collect read statistics (idxstats) on files.


#paths
raw_files2 = list.files('/home/renseb01/Documents/lord/raw_data/ctDNA_Gemini/GenomeQC_data/',pattern = '*R1.fastq.gz',full.names = T)
raw_files1 = list.files('/home/renseb01/Documents/lord/raw_data/ctDNA_Gemini/IUCPQ_data/',pattern = '*R1.fastq.gz',full.names = T)
raw_files1 = raw_files1[-grep('filtered',raw_files1)]
raw_files = c(raw_files1,raw_files2)

processed_files2 = list.files('/home/renseb01/Documents/lord/processed_data/ctDNA_Gemini/GenomeQC_data/',pattern = '*R1.fq.gz',full.names = T)
processed_files1 = list.files('/home/renseb01/Documents/lord/processed_data/ctDNA_Gemini/IUCPQ_data/',pattern = '*R1.fq.gz',full.names = T)
processed_files = c(processed_files1,processed_files2)

aligned_files =  list.files('/mnt/sde/renseb01/Documents/gemini_wflow/bams_hisat/',pattern = 'sorted.bam$',full.names = T)

#results_df
results = data.frame(raw_files = raw_files, processed_files = processed_files, aligned_files = aligned_files, raw = 0, processed = 0, aligned = 0, raw_p = 0,
processed_p = 0, aligned_p = 0,insert_size = 0,coverage = 0, X1 = 0)

print(paste0('Start loop, Time is: ',Sys.time()))

#for loop
for(i in 1:nrow(results)){
  #coverage
    temp_coverage = NULL
    temp_X1 = NULL
  for (c in 1:3){
    chr = c("chr1","chr10","chr20")
    system(paste0("samtools coverage -r ",chr[c]," ",aligned_files[i]," >cov_temp"))
    cov_temp = read.table('cov_temp')
    temp_coverage = c(temp_coverage,cov_temp$V6)
    temp_X1 = c(temp_X1,cov_temp$V7)
  }

  print(paste0('Done covtemp, Time is: ',Sys.time()))
  
  #coverage
  results$coverage[i] = mean(temp_coverage)
  results$X1[i] = mean(temp_X1)


  #median insert size
  system(paste0("samtools view ",aligned_files[i] ,"| head -10000000 | awk '{print $9}' >temp_insertsize"))
  temp_insertsize = read.table('temp_insertsize')
  results$insert_size[i] = median(temp_insertsize[temp_insertsize[,1]>0,1])

  #raw
  cmd1 = paste0('samtools idxstats ',results$raw_files[i],' >idxstats.raw')
  system(cmd1)
  idxstats.raw = read.table('idxstats.raw')
  #results$raw[i] = round(idxstats.raw[1,4]/1000000,2)
  results$raw[i] = idxstats.raw[1,4]

  #processed
  cmd2 = paste0('samtools idxstats ',results$processed_files[i],' >idxstats.processed')
  system(cmd2)
  idxstats.processed = read.table('idxstats.processed')
  #results$processed[i] = round(idxstats.processed[1,4]/1000000,2)
  results$processed[i] = idxstats.processed[1,4]

  #aligned
  cmd3 = paste0('samtools idxstats ',results$aligned_files[i],' >idxstats.aligned')
  system(cmd3)
  idxstats.aligned = read.table('idxstats.aligned')
  #results$aligned[i] = round(sum(idxstats.aligned[,3])/1000000,2)
  results$aligned[i] = sum(idxstats.aligned[,3])

  print(paste0('Done ',i,', The time is: ',Sys.time()))
}

system('rm idxstats.raw idxstats.processed idxstats.aligned temp_insertsize cov_temp')

#
results$aligned_p = round(results$aligned / results$processed * 50,2)
results$raw_p = round(results$processed / results$raw * 100,2)

write.table(results,'../outDir/alignment_stats.csv')
