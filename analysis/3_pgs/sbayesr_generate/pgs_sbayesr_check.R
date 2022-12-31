#==============================================================================
#
# Check SBayesR output (script courtesy of Tian Lin)
#
#==============================================================================

#==============================================================================
# Arguments
#==============================================================================

gwas.file		<- commandArgs(trailingOnly = TRUE)[1]
result.prefix 	<- commandArgs(trailingOnly = TRUE)[2]
out_dir			<- commandArgs(trailingOnly = TRUE)[3]
trait			<- commandArgs(trailingOnly = TRUE)[4]

#==============================================================================
# Libraries and function (from Tian)
#==============================================================================

##### NOW: go to prs_analysis.Rmd file

options(bitmapType='cairo')
library(data.table)
library(ggplot2)

# Inputs (script from Tian:
# - trait:	for naming
# - gwas:	to compare effect sizes
# - prefix: to read in SBayesR output

gwas 		<- data.frame(fread(gwas.file, showProgress = F))

pred.file 	<- paste0(result.prefix, ".snpRes")

parRes.file <- read.table(paste(result.prefix,".parRes", sep=""), skip =2)
parRes.file[1:4,] <- round(parRes.file[1:4,], 3)
parRes.file[5,] <- round(parRes.file[5,])
parRes.file[7,] <- round(parRes.file[7,],3)
parRes.file[c(6,8,9),] <- format(parRes.file[c(6,8,9),], format = "e", digits = 2)

predictor = data.frame(fread(pred.file, showProgress = F))
predictor$b = gwas[match(predictor$Name, gwas$SNP) , "EFFECT"]
predictor$gwasA1 = gwas[match(predictor$Name , gwas$SNP),"A1"]
predictor$gwas.b = NA
predictor[which(predictor$A1 == predictor$gwasA1),]$gwas.b =  predictor[which(predictor$A1 == predictor$gwasA1),"b"]
predictor[which(predictor$A1 != predictor$gwasA1),]$gwas.b = (-predictor[which(predictor$A1 != predictor$gwasA1),"b"])

effect.plot <- ggplot(predictor, aes(x=gwas.b, y=A1Effect)) + 
  geom_point(size=0.5)  + 
  geom_abline(intercept=0, slope=1, color="blue")    +
  labs(title=paste0(trait, ", with ", 
  	nrow(predictor[which(predictor$A1Effect != 0),]), " SNPs != 0 , cor = ", 
  	signif(cor(x = predictor$b, y = predictor$A1Effect, use = "pairwise.complete.obs" ),2)),
  	x="GWAS marginal effect", y = "SBayesR Effect size") 

#==============================================================================
# Generate outputs
#==============================================================================

# Output the plot
ggsave(paste(out_dir, "/", trait, "_sbayesr.png", sep = ""), effect.plot, device = "png")

# Output the predictor file
write.table(predictor, paste(out_dir, "/", trait, "_sbayesr.predictor", sep = ""), 
	col.names = T, row.names = T, sep = "\t", quote = F)





