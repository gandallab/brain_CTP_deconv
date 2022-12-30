#==============================================================================
#
# Format files for OSCA
#
#==============================================================================

#------------------------------------------------------------------------------
# Directories
#------------------------------------------------------------------------------

norm_beta_dir	<- commandArgs(trailingOnly = TRUE)[1]
out_dir			<- commandArgs(trailingOnly = TRUE)[2]

#------------------------------------------------------------------------------
# Read-in data and libraries
#------------------------------------------------------------------------------

print("Loading normalised beta table (from meffil)")
load(norm_beta_dir)

print("Tidying table")
norm.beta.df <- data.frame(IID = rownames(norm.beta), norm.beta, row.names = NULL)
colnames(norm.beta.df) <- gsub("X", "", colnames(norm.beta.df))

print("Writing table")
write.table(norm.beta.df, out_dir, 
	col.names = T, row.names = F, quote = F, sep = "\t")