#!/bin/bash


#-------
# Input
#---------------------------------------------------------------------------------------------------------------
bamdir=/mnt/sde/renseb01/Documents/gemini_wflow/bams
outDir=../outDir
tmpDir=../temp
nProcesses=10
#---------------------------------------------------------------------------------------------------------------


#bamFile=$(ls -1v $bamDir/*bam | head -n $SGE_TASK_ID | tail -n 1)
#bamFile=$(ls -1v $bamDir/*bam |head -1)
#echo $bamFile 
#mkdir -p $outDir

#for loop
for bamfile in $bamdir/*bam
  do
    date
    echo $bamfile
#    samtools index $bamfile
    

    #Rscript ./genno05-count-sbs.R $bamFile $outDir $tmpDir $nProcesses
    Rscript ./i05-count-sbs.R $bamfile $outDir $tmpDir $nProcesses

  done
