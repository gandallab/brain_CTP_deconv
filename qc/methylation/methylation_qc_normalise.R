#==============================================================================
#
# Methylation QC
# With thanks to Marta Nabais for scripts and pipeline
#
#==============================================================================

#------------------------------------------------------------------------------
# Directories
#------------------------------------------------------------------------------

qcdir 		<- commandArgs(trailingOnly = TRUE)[1]
filen		<- commandArgs(trailingOnly = TRUE)[2]
pc			<- commandArgs(trailingOnly = TRUE)[3]
ncores		<- commandArgs(trailingOnly = TRUE)[4]
# out_dir <- "/30days/uqcyap3/ASD/Data/4_methylation/qc"

#------------------------------------------------------------------------------
# Read-in data and libraries
#------------------------------------------------------------------------------

# Libraries
library(meffil)

# Cores
options(mc.cores = ncores)

# Paste together directories
qc_summary_dir <- paste(qcdir, "/", "qc.summary.clean_", filen, ".Robj", sep = "")
qc_object_dir <- paste(qcdir, "/", "qc.objects.clean_", filen, ".Robj", sep = "")

# Load .Robj
load(qc_summary_dir)
load(qc_object_dir)

# Ensure pc is in numeric form
pc <- as.numeric(pc)

#------------------------------------------------------------------------------
# Normalise
#------------------------------------------------------------------------------

print("Quantile normalisation")
norm.objects <- meffil.normalize.quantiles(qc.objects, fixed.effects = "plate", number.pcs = pc) 

print("Save quantile normalisation file")
save(norm.objects, file = paste(qcdir, "/", "norm.obj_", filen, ".Robj", sep = ""))

print("Normalise samples")
norm.beta <- meffil.normalize.samples(norm.objects, 
	cpglist.remove = qc.summary$bad.cpgs$name) 

print("Save sample normalisation file")
save(norm.beta, file = paste(qcdir, "/", "norm.beta_", filen, ".Robj", sep = ""))

#------------------------------------------------------------------------------
# Generate report
#------------------------------------------------------------------------------

print("Generate normalisation report")
str(norm.objects[[1]]$samplesheet)

#You change it by running a loop

for (i in 1:length(norm.objects)){
norm.objects[[i]]$samplesheet$plate<-as.factor(norm.objects[[i]]$samplesheet$plate)
norm.objects[[i]]$samplesheet$Slide<-as.factor(norm.objects[[i]]$samplesheet$Slide)
norm.objects[[i]]$samplesheet$Sex<-as.factor(norm.objects[[i]]$samplesheet$Sex)
norm.objects[[i]]$samplesheet$sentrix_row<-as.factor(norm.objects[[i]]$samplesheet$sentrix_row)
norm.objects[[i]]$samplesheet$sentrix_col<-as.factor(norm.objects[[i]]$samplesheet$sentrix_col)
}

batch_var <- c("plate", "Slide", "sentrix_row", "sentrix_col", "Sex")
norm.parameters <- meffil.normalization.parameters(
	norm.objects,
	variables = batch_var,
	control.pcs = 1:10,
	batch.pcs = 1:10,
	batch.threshold = 0.01
)

print("Obtain methylation PCs")
pcs <- meffil.methylation.pcs(norm.beta, probe.range = 20000)
save(pcs, file = paste(qcdir, "/", "pcs.norm.beta_", filen, ".Robj", sep = ""))

print("Generate normalisation summary report")
norm.summary <- meffil.normalization.summary(norm.objects, 
	pcs = pcs, parameters = norm.parameters)
meffil.normalization.report(norm.summary, 
	output.file = paste(qcdir, "/", "normalization_report_", filen, ".html", sep = ""),
	author = "Chloe Yap", study = paste("Methylation:", filen))
