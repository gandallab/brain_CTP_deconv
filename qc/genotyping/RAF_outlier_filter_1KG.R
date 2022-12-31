frq.file = commandArgs(trailingOnly = TRUE)[1]
ref.file = commandArgs(trailingOnly = TRUE)[2]
raf.out = commandArgs(trailingOnly = TRUE)[3]

library(data.table)

print("Reading in files")
raf <- fread(frq.file, header=T)
KG.raf <- fread(ref.file, header=T, showProgress=F)

# https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html
# - column 1 (id) - variant ID
# - column 2 (position) - base pair position
# - column 3 (a0) - allele labeled '0' in .hap file. This is the REF (or reference) allele.
# - column 4 (a1) - allele labeled '1' in .hap file. This is the ALT (or alternate) allele.
# - column 5 (TYPE) - SNP/INDEL/SV denotes type of biallelic variant
# - column 6-10 (AFR, AMR, EAS, EUR, SAS) - ALT allele frequency in continental groups. The mapping of populations to groups is given below.
# - column 11 (ALL) - ALT allele frequency across all 2,504 samples
# KG.raf <- fread("/u/home/c/chloeyap/shared-gandalm/brain_CTP/Data/genotyping/ref/1000GP_Phase3_combined_rsid.legend")
# colnames(KG.raf) <- c("id", "chr", "position", "a0", "a1", "TYPE", "AFR", "AMR", "EAS", "EUR", "SAS", "ALL")
# raf <- fread("/u/home/c/chloeyap/shared-gandalm/brain_CTP/Data/genotyping/Jaffe2018/Illumina_1M/LIBD_szControl_Genotype_QCd_Illumina_1M_Jaffe2018_genoQC-updated_MERGE_RAF.afreq")
# raf.out <- "/u/home/c/chloeyap/shared-gandalm/brain_CTP/Data/genotyping/Jaffe2018/Illumina_1M/LIBD_szControl_Genotype_QCd_Illumina_1M_Jaffe2018_genoQC-updated_MERGE_RAF_1KG"

# KG.raf$id <- unlist(lapply(strsplit(KG.raf$id, "\\:"), function(x) x[1]))

# Collate into one dataframe
raf$KGPanelAF = KG.raf[match(raf$ID, KG.raf$id),"ALL"]
raf$KGref = KG.raf[match(raf$ID, KG.raf$id),"a0"]

# Find the alleles where REF/ALT are the wrong way around, and recalculate ALT_FREQS
print("Finding mis-matched reference alleles")
foo <- which(raf$ALT == raf$KGref)
raf[foo,"ALT_FREQS"] <- 1-raf[foo,"ALT_FREQS"]
raf[foo,"ALT"] <- raf[foo,"REF"]
raf[foo,"REF"] <- raf[foo,"KGref"] # essentially, swapping columns around as ALT == KGref

# Output file to force plink to recode
swap <- raf[foo,c("ID", "REF")]

out <- raf[which(((raf$ALT_FREQS-raf$KGPanelAF)^2 > 0.2^2)  |  (raf$REF!=raf$KGref) ), ]$ID

print("Plotting allele frequencies")
png(paste(raf.out, "_raf.png", sep = ""))
plot(raf$KGPanelAF,raf$ALT_FREQS)
abline(0.2,1)
abline(-0.2,1)
dev.off()

print("Writing outputs")
write.table(out, file=paste(raf.out, ".afreqoutlier", sep = ""), 
	quote=F, sep="\t", row.names=F, col.names=F)
write.table(swap, file=paste(raf.out, ".afreq2recode", sep = ""), 
	quote=F, sep="\t", row.names=F, col.names=F)