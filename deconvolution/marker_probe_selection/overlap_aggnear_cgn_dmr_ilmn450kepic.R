#==============================================================================
#
# Overlap DMR CGNs with 450K sites
#
#==============================================================================

#==============================================================================
# Directories
#==============================================================================

cgn_dmr_dir 	<- commandArgs(trailingOnly = TRUE)[1]
ilmn450k_dir 	<- commandArgs(trailingOnly = TRUE)[2]
bp_thresh 		<- commandArgs(trailingOnly = TRUE)[3]
coverage 		<- commandArgs(trailingOnly = TRUE)[4]
filen			<- commandArgs(trailingOnly = TRUE)[5]
out_dir			<- commandArgs(trailingOnly = TRUE)[6]

#==============================================================================
# Libraries
#==============================================================================

library(data.table)
library(GenomicRanges)
library(tidyverse)

#==============================================================================
# Read in data
#==============================================================================

print(paste("cgn_dmr_dir: ", cgn_dmr_dir, sep = " "))
print(paste("ilmn450k_dir: ", ilmn450k_dir, sep = " "))
print(paste("bp_thresh: ", bp_thresh, sep = " "))
print(paste("coverage thresh: ", coverage, sep = " "))
print(paste("filen: ", filen, sep = " "))
print(paste("out_dir: ", out_dir, sep = " "))

print("Read in data")

# CGN DMR file
cgn_dmr.gr <- readRDS(cgn_dmr_dir)

# manifest
ilmn450k.gr <- readRDS(ilmn450k_dir)
ilmn450k_rmmask.gr <- ilmn450k.gr[which(ilmn450k.gr$MASK_general == FALSE),c("gene")]
ilmn450k_rmmask.gr$Probe <- names(ilmn450k_rmmask.gr)
names(ilmn450k_rmmask.gr) <- NULL

bp_thresh <- as.numeric(bp_thresh)
coverage <- as.numeric(coverage)

#==============================================================================
# Sliding window approach
#==============================================================================

print(paste("Sliding window, based on threshold:", bp_thresh, "bp", sep = " "))
ilmn450k_rmmask_resize.gr <- resize(ilmn450k_rmmask.gr, bp_thresh/2, fix = "center") 
tmp.overlap <- findOverlaps(ilmn450k_rmmask_resize.gr, cgn_dmr.gr)
length(unique(queryHits(tmp.overlap))) # number of unique 450K probes covered
length(tmp.overlap) # number of sequenced sites included

ilmn450k_rmmask_resize_sum.tmp <- data.frame(Probe = ilmn450k_rmmask_resize.gr[queryHits(tmp.overlap)]$Probe,
          gene = ilmn450k_rmmask_resize.gr[queryHits(tmp.overlap)]$gene,
          methyl_C = cgn_dmr.gr[subjectHits(tmp.overlap)]$methyl_C, 
          methyl_unmethyl_C = cgn_dmr.gr[subjectHits(tmp.overlap)]$methyl_unmethyl_C)

ilmn450k_rmmask_resize_sum <- ilmn450k_rmmask_resize_sum.tmp %>% 
  group_by(Probe) %>% dplyr::summarise(methyl_C = sum(methyl_C), methyl_unmethyl_C = sum(methyl_unmethyl_C))

#==============================================================================
# Apply coverage threshold
#==============================================================================

ilmn450k_rmmask_resize_sum <- ilmn450k_rmmask_resize_sum %>% filter(methyl_unmethyl_C >= coverage)

#==============================================================================
# Add beta column
#==============================================================================

print("Add methylation mean column")
# Equation: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
# and from meffil: https://rdrr.io/github/perishky/meffil/src/R/get-beta.r
# - note alpha of 1 is default to make methylation estimate more stable
ilmn450k_rmmask_resize_sum$M_value <- log2((ilmn450k_rmmask_resize_sum$methyl_C + 1) / (ilmn450k_rmmask_resize_sum$methyl_unmethyl_C - ilmn450k_rmmask_resize_sum$methyl_C + 1))

print("Add methylation beta column")
# Equation: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
# and from meffil: https://rdrr.io/github/perishky/meffil/src/R/get-beta.r
# - note pseudocount of 100 to make methylation estimate more stable
ilmn450k_rmmask_resize_sum$beta <- ilmn450k_rmmask_resize_sum$methyl_C / (ilmn450k_rmmask_resize_sum$methyl_unmethyl_C + 100)

#==============================================================================
# Write output
#==============================================================================

print("Save output to .rds")
saveRDS(ilmn450k_rmmask_resize_sum, paste(out_dir, "/", filen, sep = ""))