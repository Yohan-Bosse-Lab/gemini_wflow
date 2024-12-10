### This repository contains code from Bruhm et al., Nat Genet (In press).

### This an updated version of the original repo
*  [repo](https://github.com/cancer-genomics/gemini_wflow)
*  [paper](https://www.nature.com/articles/s41588-023-01446-3)


### Modifications
* add scripts for STAR alignment from `.fastq` files prior to the gemini pipeline.
* add `print()` to facilitate debugging and time/CPU consumption.
* metadata re-structure
* minor refactoring 


### Instructions (Notes)
* The code in `../scripts` generates GEMINI scores for a Control/Case cohort.
* Shell scripts were run when available, or directly from the R console if not.
* Note that file paths in scripts are relative to the directory containing the script.
* Software versions, including for R, python, and their packages are included in ./software_versions (but this is likely outdated).
* The code assumes you have `..fastq` files, which are then aligned (against hg19) and (indexed) `..bams` are in the `.bam` directly.
