#==============================================================================
# Script to generate heterozygosity and inbreeding coefficient
# plots in R for genotyping QC
#==============================================================================

# Directories
het_dir 	<- commandArgs(trailingOnly = TRUE)[1]
miss_dir 	<- commandArgs(trailingOnly = TRUE)[2]
anc_dir 	<- commandArgs(trailingOnly = TRUE)[3]
out_dir 	<- commandArgs(trailingOnly = TRUE)[4]
filen		<- commandArgs(trailingOnly = TRUE)[5]

# Read in data
het <- read.table(het_dir, as.is = T, header = T)
miss <- read.table(miss_dir, as.is = T, header = T)
colnames(miss) <- c("FID", "IID", "MISSING_CT", "OBS_CT", "F_MISS")

head(het)
head(miss)

# Libraries
library(plyr)
library(ggplot2)

# Calculate mean heterozygosity
# mean het = N-O/N
het$Het_mean <- (het$N.NM. - het$O.HOM.)/het$N.NM.

# Join the files together
d <- join(het, miss, by = c("IID", "FID"), type = "inner")

print("View joined files")
head(d)

# Add ancestry
if (anc_dir != "NA") {
	anc <- read.table(anc_dir, as.is = T, header = T)
	d <- join(d, anc, by = c("IID", "FID"))
	head(d)
}

# Plots 
# 1. mean heterozygosity vs missingness
# 2. inbreeding coefficient (F) vs missingness
if (anc_dir == "NA") {

	print("Not including ancestry information")

	print("Heterozygosity plot")
	het_plot <- ggplot(d, aes(x=MISSING_CT, y=Het_mean)) + 
	    geom_point()
	het_plot
	ggsave(paste(out_dir, "/hetmean_vs_missingness", filen, ".png", sep = ""))

	print("Inbreeding coefficient plot")
	inb_plot <- ggplot(d, aes(x=MISSING_CT, y=F)) + 
    geom_point()
	inb_plot
	ggsave(paste(out_dir, "/inbreeding_vs_missingness", filen, ".png", sep = ""))

} else if (anc_dir != "NA") {

	print("Including ancestry information")
	
	print("Heterozygosity plot")	
	het_plot <- ggplot(d, aes(x=MISSING_CT, y=Het_mean, color=population, alpha=0.8)) + 
	    geom_point()
	het_plot
	ggsave(paste(out_dir, "/hetmean_vs_missingness", filen, ".png", sep = ""))	

	print("Inbreeding coefficient plot")
	inb_plot <- ggplot(d, aes(x=MISSING_CT, y=F, color=population, alpha=0.8)) + 
    geom_point()
	inb_plot
	ggsave(paste(out_dir, "/inbreeding_vs_missingness", filen, ".png", sep = ""))

}



