#==============================================================================
#
# Methylation QC - remove outliers
#
#==============================================================================

#------------------------------------------------------------------------------
# Directories
#------------------------------------------------------------------------------

filen <- "Jaffe2018"
filen <- "Jaffe2018_age18"
filen <- "Jaffe2018_age18_plate"
#filen <- "Jaffe2018_fetal"
filen <- "Jaffe2018_fetal_plate"
filen <- "Jaffe2018_age0"
qcdir <- "~/shared-gandalm/brain_CTP/Data/methylation/Jaffe2018/processed"

filen <- "ROSMAP"
qcdir <- "~/shared-gandalm/brain_CTP/Data/methylation/ROSMAP/processed"

qc_summary_dir <- paste(qcdir, "/", "qc.summary_", filen, ".Robj", sep = "")
qc_object_dir <- paste(qcdir, "/", "qc.objects_", filen, ".Robj", sep = "")

#------------------------------------------------------------------------------
# Read-in data and libraries
#------------------------------------------------------------------------------

load(qc_summary_dir)
load(qc_object_dir)

library(meffil)

#------------------------------------------------------------------------------
# Identify outliers (SHOULD HAVE CHECKED MANUALLY BEFORE)
#------------------------------------------------------------------------------

outlier <- qc.summary$bad.samples
table(outlier$issue)

write.table(outlier, file = paste(qcdir, "/", "outlier_", filen, ".id", sep = ""), 
	col.names = T, row.names = F, quote = F, sep = "\t")

#------------------------------------------------------------------------------
# Remove outliers from qc.objects
#------------------------------------------------------------------------------

length(qc.objects) # check
qc.objects <- meffil.remove.samples(qc.objects, outlier$sample.name)
length(qc.objects) # check

save(qc.objects, file = paste(qcdir, "/", "qc.objects.clean_", filen, ".Robj", sep = ""))

#------------------------------------------------------------------------------
# Re-generate summary
#------------------------------------------------------------------------------

# Get QC parameters

qc.parameters <- meffil.qc.parameters(
	beadnum.samples.threshold             = 0.1,
	detectionp.samples.threshold          = 0.1,
	detectionp.cpgs.threshold             = 0.1, 
	beadnum.cpgs.threshold                = 0.1,
	sex.outlier.sd                        = 5,
	snp.concordance.threshold             = 0.95,
	sample.genotype.concordance.threshold = 0.8
)

qc.summary <- meffil.qc.summary(qc.objects, 
	parameters = qc.parameters)

save(qc.summary, file = paste(qcdir, "/", "qc.summary.clean_", filen, ".Robj", sep = ""))

#------------------------------------------------------------------------------
# Check PCs
#------------------------------------------------------------------------------

y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot,filename = paste(qcdir, "/pc.fit_", filen, ".pdf", sep = ""),
	height = 6, width = 6)

#------------------------------------------------------------------------------
# Extract cell counts
#------------------------------------------------------------------------------

cc <- t(sapply(qc.objects, function(obj) obj$cell.counts$counts))
cc <- data.frame(IID = row.names(cc), cc)
write.table(cc, paste(qcdir, "/cellcounts_", filen, ".txt", sep = ""),
	col.names = T,row.names = F, quote = F, sep = "\t")
