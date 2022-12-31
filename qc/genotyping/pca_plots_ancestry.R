#==============================================================================
#
# PCA plots
# Scripts from Tian Lin
#
#==============================================================================

str.file    = commandArgs(trailingOnly = TRUE)[1]
ref.file    = commandArgs(trailingOnly = TRUE)[2]
pc.file     = commandArgs(trailingOnly = TRUE)[3]
cohort      = commandArgs(trailingOnly = TRUE)[4]
grm_prefix  = commandArgs(trailingOnly = TRUE)[5]
out         = commandArgs(trailingOnly = TRUE)[6]

# str.file <- "/QRISdata/Q0286/Chloe/ASD/Data/2_GWAS/1000GP_Phase3.sample"
# ref.file <- "/QRISdata/Q0286/Chloe/ASD/Data/2_GWAS/1000G_phase3_20130502_combined_snpsonly.05.common_pca2.eigenvec"
# pc.file <- "/90days/uqcyap3/ASD/Data/AUC_GSA_v1v2_recluster_merge/GSAv1v2_cleaned_aut_x_flipped_fixed_RAF_HRC.05.common_pca2_ACRC.proj.eigenvec"
# cohort <- "ACRC"
# grm_prefix <- "/90days/uqcyap3/ASD/Data/AUC_GSA_v1v2_recluster_merge/GSAv1v2_cleaned_aut_x_flipped_fixed_RAF_HRC.05.common_ACRC"
# out <- "/90days/uqcyap3/ASD/Data/AUC_GSA_v1v2_recluster_merge"

#==============================================================================
# Load files and libraries
#==============================================================================

library(ggplot2)
library(reshape2)
library(colorspace)

str <- read.table(str.file, header = T, as.is = T)
ref <- read.table(ref.file, as.is = T)
pc <- read.table(pc.file, as.is = T)

#==============================================================================
# Make plot
#==============================================================================

ref$str = str[match(ref$V2, str$ID),"GROUP"]
pc$str <- cohort

colnames(ref) = c("FID", "IID", "PC1", "PC2", "str")
colnames(pc) = c("FID", "IID", "PC1", "PC2", "str")

data <- rbind(data.frame(ref), data.frame(pc))

# Create plot and save
fig <- ggplot(data=data, aes(x=PC1, y=PC2, color=str)) + 
  geom_point(size=0.7) +
  scale_color_manual(values=rainbow_hcl(length(unique(data$str)), start = 0, end = 280, c = 100, l = 65))
ggsave(paste(out, "pc_plot.png", sep="/"), fig, height = 9, width = 9)

#------------------------------------------------------------------------------
# Cluster samples in 1000G populations
#------------------------------------------------------------------------------
# Define the threshold of PC1 and PC2 based on the mean value and 4 standard deviation of the PCs in that population. 
# In the figure below, blue lines are the thresholds of EURs, and red lines are for EAS.

## ancestry assignment
pc1_SD_eur <-   sd(data[data$str=="EUR",3])
pc2_SD_eur <-   sd(data[data$str=="EUR",4])   
pc1_mn_eur <- mean(data[data$str=="EUR",3])  
pc2_mn_eur <- mean(data[data$str=="EUR",4]) 

pc1_SD_eas <-   sd(data[data$str=="EAS",3])
pc2_SD_eas <-   sd(data[data$str=="EAS",4])   
pc1_mn_eas <- mean(data[data$str=="EAS",3])  
pc2_mn_eas <- mean(data[data$str=="EAS",4])     

pc1_SD_afr <-   sd(data[data$str=="AFR",3])
pc2_SD_afr <-   sd(data[data$str=="AFR",4])   
pc1_mn_afr <- mean(data[data$str=="AFR",3])  
pc2_mn_afr <- mean(data[data$str=="AFR",4])     

pc1_SD_sas <-   sd(data[data$str=="SAS",3])
pc2_SD_sas <-   sd(data[data$str=="SAS",4])   
pc1_mn_sas <- mean(data[data$str=="SAS",3])  
pc2_mn_sas <- mean(data[data$str=="SAS",4])     

n <- 4

x1_eur <- pc1_mn_eur-n*pc1_SD_eur
x2_eur <- pc1_mn_eur+n*pc1_SD_eur
y1_eur <- pc2_mn_eur-n*pc2_SD_eur
y2_eur <- pc2_mn_eur+n*pc2_SD_eur

x1_eas <- pc1_mn_eas-n*pc1_SD_eas
x2_eas <- pc1_mn_eas+n*pc1_SD_eas
y1_eas <- pc2_mn_eas-n*pc2_SD_eas
y2_eas <- pc2_mn_eas+n*pc2_SD_eas

x1_afr <- pc1_mn_afr-n*pc1_SD_afr
x2_afr <- pc1_mn_afr+n*pc1_SD_afr
y1_afr <- pc2_mn_afr-n*pc2_SD_afr
y2_afr <- pc2_mn_afr+n*pc2_SD_afr

x1_sas <- pc1_mn_sas-n*pc1_SD_sas
x2_sas <- pc1_mn_sas+n*pc1_SD_sas
y1_sas <- pc2_mn_sas-n*pc2_SD_sas
y2_sas <- pc2_mn_sas+n*pc2_SD_sas

# Create and save plot
new.fig <- fig +
 geom_vline(xintercept = x1_eur, color='royalblue', size=0.5)  +
 geom_vline(xintercept = x2_eur, color='royalblue', size=0.5)  +
 geom_hline(yintercept = y1_eur, color='royalblue', size=0.5)  +
 geom_hline(yintercept = y2_eur, color='royalblue', size=0.5)  +
 geom_vline(xintercept = x1_eas, color='seagreen3', size=0.5) +
 geom_vline(xintercept = x2_eas, color='seagreen3', size=0.5) +
 geom_hline(yintercept = y1_eas, color='seagreen3', size=0.5) +
 geom_hline(yintercept = y2_eas, color='seagreen3', size=0.5) +
 geom_vline(xintercept = x1_afr, color='orange3', size=0.5) +
 geom_vline(xintercept = x2_afr, color='orange3', size=0.5) +
 geom_hline(yintercept = y1_afr, color='orange3', size=0.5) +
 geom_hline(yintercept = y2_afr, color='orange3', size=0.5) +
 geom_vline(xintercept = x1_sas, color='purple', size=0.5) +
 geom_vline(xintercept = x2_sas, color='purple', size=0.5) +
 geom_hline(yintercept = y1_sas, color='purple', size=0.5) +
 geom_hline(yintercept = y2_sas, color='purple', size=0.5)

ggsave(paste(out, "pc_plot_ancestry.png", sep="/"), new.fig, height = 9, width = 9)


#------------------------------------------------------------------------------
# Population assignment 
#------------------------------------------------------------------------------

pc$population = "other"
pc$population[which(pc$PC1>x1_eur & pc$PC1<x2_eur  & pc$PC2>y1_eur & pc$PC2<y2_eur)] = "EUR"
pc$population[which(pc$PC1>x1_eas & pc$PC1<x2_eas  & pc$PC2>y1_eas & pc$PC2<y2_eas)] = "EAS"
pc$population[which(pc$PC1>x1_afr & pc$PC1<x2_afr  & pc$PC2>y1_afr & pc$PC2<y2_afr)] = "AFR"
pc$population[which(pc$PC1>x1_sas & pc$PC1<x2_sas  & pc$PC2>y1_sas & pc$PC2<y2_sas)] = "SAS"

# Write output table
write.table(pc, file = paste(out, "pc_and_population.txt", sep="/"), 
	row.names=F, col.names=T, quote=F, sep="\t")

