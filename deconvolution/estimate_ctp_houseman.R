#==============================================================================
#
# Estimate CTP using Houseman 2012 method + sequencing data (Luo et al. 2022)
#
#==============================================================================

#------------------------------------------------------------------------------
# FUNCTIONS: READ IN FIRST

# Houseman 2012 deconvolution implemented in minfi (Aryee et al. 2014)
# minfi::projectCellType within estimateCellCounts.R 
# (https://rdrr.io/bioc/minfi/src/R/estimateCellCounts.R)
# this function preferred here as requiring just methylation beta matrix 
# is more flexible

# Inputs:
# Y: bulk methylation beta matrix
# coefCellType: reference, matrix of m rows (CpG probes) x n columns (cell-types)

projectCellType <- function(Y, coefCellType, contrastCellType = NULL,
                            nonnegative = TRUE, lessThanOne = TRUE) {
    if (is.null(contrastCellType)) {
        Xmat <- coefCellType
    } else {
        Xmat <- tcrossprod(coefCellType, contrastCellType)
    }

    nCol <- dim(Xmat)[2]
    if (nCol == 2) {
        Dmat <- crossprod(Xmat)
        mixCoef <- t(
            apply(Y, 2, function(x) solve(Dmat, crossprod(Xmat, x))))
        colnames(mixCoef) <- colnames(Xmat)
        return(mixCoef)
    } else {
        nSubj <- dim(Y)[2]

        mixCoef <- matrix(0, nSubj, nCol)
        rownames(mixCoef) <- colnames(Y)
        colnames(mixCoef) <- colnames(Xmat)

        if (nonnegative) {
            if (lessThanOne) {
                Amat <- cbind(rep(-1, nCol), diag(nCol))
                b0vec <- c(-1, rep(0, nCol))
            } else {
                Amat <- diag(nCol)
                b0vec <- rep(0, nCol)
            }
            for (i in seq_len(nSubj)) {
                obs <- which(!is.na(Y[,i]))
                Dmat <- crossprod(Xmat[obs,])
                mixCoef[i,] <- solve.QP(
                    Dmat = Dmat,
                    dvec = crossprod(Xmat[obs,], Y[obs,i]),
                    Amat = Amat,
                    bvec = b0vec)$sol
            }
        } else {
            for (i in seq_len(nSubj)) {
                obs <- which(!is.na(Y[,i]))
                Dmat <- crossprod(Xmat[obs,])
                mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
            }
        }
        mixCoef
    }
}

# Seiler-Vellame et al. 2022 biorXiv https://www.biorxiv.org/content/10.1101/2022.06.15.496235v1

# Inputs:
# applyIndex: getErrorPerSample() implemented within sapply loop to iterate through each matrix row
# predictedIN: output cell-type proportions
# coefDataIN: reference cell-type profile matrix
# betasBulkIN: bulk methylation matrix

getErrorPerSample = function(applyIndex,
                                 predictedIN = counts,
                                 coefDataIN = coefs,
                                 betasBulkIN = coefdat){
      
    trueBulk = matrix(ncol = 1, nrow = nrow(coefDataIN), data = 0)
      
    RMSE = function(m, o){
        sqrt(mean((m - o)^2))
    }
      
    for (i in 1:ncol(coefDataIN)){
        
        trueBulk[,1] = trueBulk[,1] + coefDataIN[,i]*predictedIN[applyIndex,i]
    }
      
    betasBulkIN = t(apply(betasBulkIN, 1, function(x){x[is.na(x)] = 0; return(x)}))
      
    error = RMSE(trueBulk, betasBulkIN[,applyIndex])
    return(error)
}
#------------------------------------------------------------------------------

#==============================================================================
# Deconvolution code
#==============================================================================

#------------------------------------------------------------------------------
# Libraries
#------------------------------------------------------------------------------

library(minfi)
library(quadprog)
library(csSAM)
library(data.table)
library(tidyverse)
library(ggplot2)
library(forcats)
library(dplyr)

#------------------------------------------------------------------------------
# Bulk data directory (DNA methylation beta matrix, batch corrected)
#------------------------------------------------------------------------------

input <- "txt"

# Jaffe
data_dir <- "~/shared-gandalm/brain_CTP/Data/methylation/Jaffe2018/analysis"
filen <- "Jaffe2018_age0_aut_mask_t"
meth_dir <- paste(data_dir, "/", filen, ".txt", sep = "")

# ASD brain
data_dir <- "~/shared-gandalm/brain_CTP/Data/methylation/ASD_methylation_brain/analysis"
filen <- "KCL_R01MH094714_ASD_Illumina450K_PFC"
meth_dir <- paste(data_dir, "/", filen, ".csv", sep = "")

# ROSMAP (own batch correction)
data_dir <- "~/shared-gandalm/brain_CTP/Data/methylation/ROSMAP/analysis"
filen <- "ROSMAP_aut_mask_t"
meth_dir <- paste(data_dir, "/", filen, ".txt", sep = "")

# BDR AZD
input <- "rdat"
data_dir <- "~/shared-gandalm/brain_CTP/Data/methylation/BDR_DNAm_AZD_HannonMills"
filen <- "BDR_DNAm_FINAL"
meth_dir <- paste(data_dir, "/", filen, ".rdat", sep = "")

#------------------------------------------------------------------------------
# Reference data directory
#------------------------------------------------------------------------------

# Coefficients from sequencing data
coefs_seq_dir <- "~/shared-gandalm/brain_CTP/Data/reference_cell_profile/Luo2020/Luo2020_extremes_dmr_ilmn450kepic_aggto100bp_c10_cell7_split6040all_beta.rds"

#------------------------------------------------------------------------------
# Read-in + QC of bulk + reference
#------------------------------------------------------------------------------

# Methylation data: read-in + format
if (input == "txt") {

    meth.tmp <- fread(meth_dir)
    meth <- as.matrix(meth.tmp[,2:ncol(meth.tmp)])
    rownames(meth) <- data.frame(meth.tmp)[,1]
    meth <- meth[order(rownames(meth)),]

} else if (input == "rdat") {

    load(meth_dir)
    meth <- betas
}

# Reference probes, order, and overlap with methylation matrix
# Formatted as a list here, as it gives flexibility to deconvolve 
# 1 methylation beta matrix with multiple different reference
# panels for comparison
coefs_seq <- readRDS(coefs_seq_dir)
coefs_seq.ls <- list(coefs_seq)
coefs_seq.ls <- lapply(coefs_seq.ls, function(x) x[order(rownames(x)),])
coefs_seq.ls <- lapply(coefs_seq.ls, function(x) x[which(rownames(x) %in% rownames(meth)),])

#==============================================================================
# DECONVOLUTION
#==============================================================================

# Deconvolve CTP
ctp_seq.ls <- lapply(coefs_seq.ls, function(x) projectCellType(meth[which(rownames(meth) %in% (rownames(x))), ], x, lessThanOne = TRUE))

# Add errors
ctp_seq_error.ls <- lapply(1:length(ctp_seq.ls), function(y) data.frame(IID = rownames(ctp_seq.ls[[y]]), ctp_seq.ls[[y]], error =
    sapply(1:nrow(ctp_seq.ls[[y]]), function(x) getErrorPerSample(x, 
        ctp_seq.ls[[y]], 
        coefs_seq.ls[[y]][which(rownames(coefs_seq.ls[[y]]) %in% rownames(meth)),], 
        meth[which(rownames(meth) %in% (rownames(coefs_seq.ls[[y]]))), ]))))
names(ctp_seq_error.ls) <- names(ctp_seq.ls)

# Write output
saveRDS(ctp_seq_error.ls, paste(data_dir, "/", filen, "_Luo2020_extremes_dmr_ilmn450kepic_aggto100bp_c10_cell7_split6040all_houseman_ls.rds", sep = ""))

# Plots
# seq
ctp_seq.long.ls <- lapply(ctp_seq_error.ls, function(x) melt(x, id.var = "IID", variable.name = "celltype", value.name = "CTP"))
# - add comparison identifier
compar <- names(ctp_seq.long.ls)
compar <- c("markers_Luo2020")
tmp.ls <- mapply(cbind, ctp_seq.long.ls, "comparison"=compar, SIMPLIFY=F)
# - append list elements
ctp_seq.long.df <- do.call(rbind, tmp.ls)
ctp_seq.gg <- ctp_seq.long.df %>% 
        ggplot(aes(x=celltype, y=CTP, colour=celltype)) + 
        geom_boxplot() + 
        theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
        facet_wrap(~comparison) +
        geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray20")
ggsave(paste(data_dir, "/", filen, "_Luo2020_extremes_dmr_ilmn450kepic_aggto100bp_c10_cell7_split6040all_houseman_ls.png", sep = ""), ctp_seq.gg)