library(tidyr)
suppressMessages(library(SummarizedExperiment))

# Reading in the 100kb bins and creating 2.5Mb bins
bins <- readRDS("../outDir/bins/bins_100kb.rds")
merge.grps <- readRDS("../outDir/bins/merge_grps_100kb_to_2.5mb_bins.rds")
bins <- split(bins, merge.grps) %>% reduce %>% unlist
names(bins) <- paste0("bin", 1:length(bins))

#-------------------------------------------------------------------------------
# Creating a RangedSumamrizedExperiment to hold bins, bin counts, and metadata
#-------------------------------------------------------------------------------
ctDir <- "../outDir/cts_hisat"
#ctDir = '../outDir/cts_trimseq_star'
#ctDir <- "../outDir/cts"
ctFiles <- list.files(ctDir,pattern = 'rds')
ctFiles = ctFiles[grep('NS',ctFiles,invert = F)]
#ctFiles  = ctFiles[c(21,26,30)]
print(ctFiles)
ct.list <- lapply(ctFiles, function(t) readRDS(file.path(ctDir, t)))

for (i in 1:length(ct.list)) {
  split.meta <- split(data.frame(mcols(ct.list[[i]])), merge.grps)
  meta <- lapply(split.meta, function(t) colSums(t)) %>% do.call("rbind", .)
  tmp.bins <- bins
  mcols(tmp.bins) <- meta
  ct.list[[i]] <- tmp.bins
}

ct_grps <- colnames(mcols(ct.list[[1]])) %>% .[-grep("known", .)]
assay.list <- vector("list", length(ct_grps)) %>% setNames(., ct_grps)
for (i in 1:length(ct_grps)) {
  x <- lapply(ct.list, function(t) mcols(t)[, ct_grps[i]]) %>% do.call("cbind", .)
  dimnames(x) <- list(names(bins),  gsub("_cts.rds", "", ctFiles))
  assay.list[[ct_grps[i]]] <- x
}

# Reading in sample metadata
#meta.df = read.csv('../data/metadata_new.csv')[1:31,]
meta.df = data.frame(id = colnames(assay.list[[1]]), patient_type = colnames(assay.list[[1]]), x = 'Case')


#meta.df = read.csv('../outDir/cts/metadata.csv')

# Creating a RangedSummarizedExperiment to hold the data
se <- SummarizedExperiment(assays = assay.list,
                           rowRanges = bins,
                           colData = meta.df)
#--------------------------------------------------------------------

#fixed_model bins
fixed_bins = readxl::read_xlsx('../data/fixed_bins.xlsx')

total_bins = table(seqnames(bins))

pos.bins = NULL
neg.bins = NULL 

for(i in 1:22){
  
  temp_cases = paste0('bin',as.numeric(strsplit(fixed_bins$orange[i],',')[[1]]) + (sum(total_bins[1:i]) - total_bins[i]))
  temp_controls = paste0('bin',as.numeric(strsplit(fixed_bins$green[i],',')[[1]]) + (sum(total_bins[1:i]) - total_bins[i]))
  
  pos.bins = c(pos.bins,temp_cases)  #orange is enriched in cancer
  neg.bins = c(neg.bins,temp_controls) #green is enriched in normal tissue
}

#
#pos.bins = paste0('bin',sample(1:sum(total_bins),114))  #orange is enriched in cancer
#neg.bins = paste0('bin',sample(1:sum(total_bins),114)) #green is enriched in normal tissue

#--------------------------------------------------------------------
# Generating regional differences in mutation frequencies using LOO
#--------------------------------------------------------------------
# Defining the types of mutations to analyze
types <- list(cg2at = c("cg", "cg2at"), cg2gc = c("cg", "cg2gc"), cg2ta = c("cg", "cg2ta"),
              ta2at = c("ta", "ta2at"), ta2cg = c("ta", "ta2cg"), ta2gc = c("ta", "ta2gc"))

# Creating empty columns in 'se' for each mutation type
feature.names <- names(types)
m <- matrix(nrow = ncol(se), ncol = length(feature.names), dimnames = list(NULL, paste0("multimf_", feature.names)))
colData(se) <- cbind(colData(se), m)

# Creating an empty list to save bins sets
metadata(se) <- vector("list", ncol(se)) %>% setNames(., colnames(se))
bin.list <- list(pos.bins = NA, neg.bins = NA)
for (i in 1:length(metadata(se))) {
  metadata(se)[[i]] <- list(cg2at = bin.list, cg2gc = bin.list, cg2ta = bin.list, ta2at = bin.list, ta2cg = bin.list, ta2gc = bin.list)
}

prints = data.frame()


for (i in 1:ncol(se)) {
  p1 = noquote(paste0(i, "/", ncol(se), ', Time is: ',Sys.time()))
  se.loo <- se[,-i]
  se.lo <- se[,i]

  for (j in 1:length(types)) {
    e.pyR1 <- paste0(types[[j]][1], ".pyR1")
    e.puR1 <- paste0(types[[j]][1], ".puR1")
    v.pyR1 <- paste0(types[[j]][2], ".pyR1") #variants?
    v.puR1 <- paste0(types[[j]][2], ".puR1")

    # Computing the regional mutation frequency in the held out sample using the bins with the largest differential mutation density
    if (names(types[j]) == "cg2at") {
      pos.dens.lo <- sum(assays(se.lo)[[v.pyR1]][pos.bins,]) / sum(assays(se.lo)[[e.pyR1]][pos.bins,])
      neg.dens.lo <- sum(assays(se.lo)[[v.pyR1]][neg.bins,]) / sum(assays(se.lo)[[e.pyR1]][neg.bins,])
      
      p2 = sum(assays(se.lo)[[v.pyR1]][pos.bins,])
      p3 = sum(assays(se.lo)[[v.pyR1]][neg.bins,])
      p4 = sum(assays(se.lo)[[e.pyR1]][pos.bins,])
      p5 = sum(assays(se.lo)[[e.pyR1]][neg.bins,])
      p6 = signif(pos.dens.lo*1e6,6)
      p7 = signif(neg.dens.lo*1e6,6)
      
      prints = rbind(prints,c(p1,p2,p3,p4,p5,p6,p7))
      colnames(prints) = c('Time','v_pos','v_neg','e_pos','e_neg','rate_pos','rate_neg')
     # print(mean(assays(se.lo)[[e.pyR1]][pos.bins,]))
     # print(mean(assays(se.lo)[[e.pyR1]][neg.bins,]))
      
    } else {
      pos.dens.lo <- (sum(assays(se.lo)[[v.pyR1]][pos.bins,]) + sum(assays(se.lo)[[v.puR1]][pos.bins,])) / (sum(assays(se.lo)[[e.pyR1]][pos.bins,]) + sum(assays(se.lo)[[e.puR1]][pos.bins,]))
      neg.dens.lo <- (sum(assays(se.lo)[[v.pyR1]][neg.bins,]) + sum(assays(se.lo)[[v.puR1]][neg.bins,])) / (sum(assays(se.lo)[[e.pyR1]][neg.bins,]) + sum(assays(se.lo)[[e.puR1]][neg.bins,]))
    }

    multimf.lo <- (pos.dens.lo - neg.dens.lo) * 1e6

    feature.name <- paste0("multimf_", names(types[j]))

    # Saving the regional mutation frequency for the held out sample
    colData(se)[i,feature.name] <- multimf.lo

    # Saving the bin sets
    metadata(se)[[i]][[names(types[j])]]$pos.bins <- pos.bins
    metadata(se)[[i]][[names(types[j])]]$neg.bins <- neg.bins

  } # End of j loop

} # End of i loop
#-------------------------------------------------------------------------------

colData(se)$multimf_cg2at



#-------------------------------------------------------------------------------
# Computing the GEMINI score [C>A] for each patient
#-------------------------------------------------------------------------------
#gemini_fit_cg2at <- glm(factor(patient_type, levels = c("Control", "Case")) ~ multimf_cg2at, data = colData(se), family = "binomial")
#se$gemini_score <- gemini_fit_cg2at$fitted.values
#print(se$gemini_score)

#gemini_df = data.frame(id = names(se$gemini_score), gemini_score = se$gemini_score, multimf_cg2at = se$multimf_cg2at)

write.csv(colData(se),'../outDir/gemini_score_hisat.csv')
#-------------------------------------------------------------------------------


