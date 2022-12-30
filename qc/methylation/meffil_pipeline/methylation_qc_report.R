#==============================================================================
#
# Methylation QC
# With thanks to Marta Nabais for scripts and pipeline
#
#==============================================================================

#------------------------------------------------------------------------------
# Directories
#------------------------------------------------------------------------------

qc_dir 		<- commandArgs(trailingOnly = TRUE)[1]
sampsheetn	<- commandArgs(trailingOnly = TRUE)[2]
filen		<- commandArgs(trailingOnly = TRUE)[3]
geno		<- commandArgs(trailingOnly = TRUE)[4]
ncores		<- commandArgs(trailingOnly = TRUE)[5]

#------------------------------------------------------------------------------
# Read-in data and libraries
#------------------------------------------------------------------------------

# Libraries
library(meffil)

# Cores
options(mc.cores = ncores)

# Samplesheet
samplesheet <- read.csv(paste(qc_dir, sampsheetn, sep = "/"), header = T, as.is = T)

#------------------------------------------------------------------------------
# Background, dye bias correction, sex prediction and cell count estimates.
#------------------------------------------------------------------------------

# Create QC object
print("Create QC object: meffil.qc()")
# Ensure featureset is correct! Or lose extra EPIC probes
qc.objects <- meffil.qc(samplesheet,
	cell.type.reference = "guintivano dlpfc", verbose = TRUE)
save(qc.objects, file = paste(qc_dir, "/qc.objects_", filen, ".Robj", sep = ""))

#------------------------------------------------------------------------------
# Add genotyping
#------------------------------------------------------------------------------

if (geno != "no_geno") {

	print("Add genotyping")
	# Assign variables for readability
	snpnames <- paste(qc_dir, "/snpnames_", filen, ".txt", sep = "")
	genoraw <- paste(qc_dir, "/genotypes_", filen, sep = "")
	
	# Extract control SNP names - full dataset
	# annotation <- qc.objects[[1]]$annotation
	writeLines(meffil.snp.names(), con=snpnames)# Extract AAB genotypes
	system(paste("/90days/uqcyap3/software/plink2", "--bfile", geno, 
		"--extract", snpnames, "--recode A",
		"--keep", paste(qc_dir, "sample.id", sep = "/"),
		"--out", genoraw, sep = " "))
	genotypes <- meffil.extract.genotypes(paste(genoraw, ".raw", sep = ""))

}

#------------------------------------------------------------------------------
# Run QC report
#------------------------------------------------------------------------------

qc.parameters <- meffil.qc.parameters(
	beadnum.samples.threshold             = 0.1,
	detectionp.samples.threshold          = 0.1,
	detectionp.cpgs.threshold             = 0.1, 
	beadnum.cpgs.threshold                = 0.1,
	sex.outlier.sd                        = 5,
	snp.concordance.threshold             = 0.95,
	sample.genotype.concordance.threshold = 0.8
)

# Generate the QC summary
print("Generate QC summary: meffil.qc.summary()")
if (geno != "no_geno") {
	qc.summary <- meffil.qc.summary(qc.objects, parameters = qc.parameters, 
		genotypes = genotypes)
} else {
	qc.summary <- meffil.qc.summary(qc.objects, parameters = qc.parameters)
}

# Save the object
save(qc.summary, file = paste(qc_dir, "/qc.summary_", filen, ".Robj", sep = ""))

# Generate .html report
print("Generate QC report: meffil.qc.report()")
qc_report_dir <- paste(qc_dir, "/qc.report_", filen, sep = "")
mkdir_command <- paste("mkdir ", qc_report_dir, sep = " ")
system(mkdir_command)

meffil.qc.report(qc.summary, 
	output.file = paste(qc_report_dir, "/qc.report_", filen, ".html", sep = ""), 
	author = "Chloe Yap", study = paste("Methylation:", filen))


