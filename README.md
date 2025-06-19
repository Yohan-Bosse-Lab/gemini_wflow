### This repository contains updated code from Bruhm et al. (2023).
* Single-molecule genome-wide mutation profiles of cell-free DNA for non-invasive detection of cancer
* [paper](https://www.nature.com/articles/s41588-023-01446-3)

### This an updated version of the original repo
*  [repo](https://github.com/cancer-genomics/gemini_wflow)

### Modifications
* add scripts for STAR alignment and HISAT2 (better) from `.fastq` files prior to the gemini pipeline.
* add script for alignment statistics.
* add `print()` to facilitate debugging and time/CPU consumption.
* metadata re-structure
* minor refactoring 
* the fixed bins were calculated directly from the Figure in the paper Bruhm et al. 2023. 

### Notes (instructions)
* The code in `../scripts` generates GEMINI scores for a Control/Case cohort.
* Shell scripts are run sequentially (from `i00` to `i08`). When both a `bash` and `R` script available, run the bash. If only R, run the R code in an R session.
* Note that file paths in scripts are relative to the directory containing the script.
* Software versions, including for R, python, and their packages are included in ./software_versions (but this is likely outdated).
* The code assumes you have `..fastq` files, which are then aligned (against hg19) and (indexed) `..bams` are in the `.bam` directly.
