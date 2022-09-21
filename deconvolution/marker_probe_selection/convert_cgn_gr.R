#==============================================================================
#
# Convert CGN file to GenomicRanges
#
#==============================================================================

#==============================================================================
# Directories
#==============================================================================

cgn_dir 		<- commandArgs(trailingOnly = TRUE)[1]
filen			<- commandArgs(trailingOnly = TRUE)[2]
out_dir			<- commandArgs(trailingOnly = TRUE)[3]

#==============================================================================
# Libraries
#==============================================================================

library(data.table)
library(GenomicRanges)

#==============================================================================
# Read in data
#==============================================================================

print("Read in data")
cgn <- fread(cgn_dir)

#==============================================================================
# GRanges format
#==============================================================================

print("GRanges format")
cgn.gr <- GRanges(
	seqnames = Rle(c(cgn$V1)),
	ranges = IRanges(start = cgn$V2, end = cgn$V2),
	strand = Rle(cgn$V3),
	CGN = cgn$V4,
	methyl_C = cgn$V5,
	methyl_unmethyl_C = cgn$V6)

#==============================================================================
# Autosomes only
#==============================================================================

print("Autosomes only to reduce size")
cgn.gr <- cgn.gr[which(seqnames(cgn.gr) %in% paste("chr", 1:22, sep = ""))]

#==============================================================================
# Add methylation beta column
#==============================================================================

print("Add methylation mean column")
# Equation: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
# and from meffil: https://rdrr.io/github/perishky/meffil/src/R/get-beta.r
# - note alpha of 1 is default to make methylation estimate more stable
cgn.gr$M_value <- log2((max(cgn.gr$methyl_C, 0) + 1) / (max(cgn.gr$methyl_unmethyl_C, 0) + 1))

print("Add methylation beta column")
# Equation: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
# and from meffil: https://rdrr.io/github/perishky/meffil/src/R/get-beta.r
# - note pseudocount of 100 to make methylation estimate more stable
cgn.gr$beta <- cgn.gr$methyl_C / (cgn.gr$methyl_C + cgn.gr$methyl_unmethyl_C + 100)

#==============================================================================
# Write output
#==============================================================================

print("Save output to .rds")
saveRDS(cgn.gr, paste(out_dir, "/", filen, sep = ""))
