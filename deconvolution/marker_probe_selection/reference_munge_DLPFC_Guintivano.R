#==============================================================================
#
# Obtain reference panel using DLPFC_450K_guintivano
# This is the standard option implemented in minfi
# Used here as a comparison reference panel
#
#==============================================================================

#==============================================================================
# Libraries
#==============================================================================

# http://bioconductor.org/packages/release/data/experiment/html/FlowSorted.DLPFC.450k.html
BiocManager::install("FlowSorted.DLPFC.450k")
library(FlowSorted.DLPFC.450k)

library(minfi)
library(genefilter) # required within pickCompProbes()

#==============================================================================
# Functions
#==============================================================================

# internal function from minfi::estimateCellCounts
# accessed at: https://rdrr.io/bioc/minfi/src/R/estimateCellCounts.R

pickCompProbes <- function(mSet, cellTypes = NULL, numProbes = 50,
                           compositeCellType = compositeCellType,
                           probeSelect = probeSelect) {
   #.isMatrixBackedOrStop(mSet)
    splitit <- function(x) {
        split(seq_along(x), x)
    }

    p <- getBeta(mSet)
    pd <- as.data.frame(colData(mSet))
    if (!is.null(cellTypes)) {
        if (!all(cellTypes %in% pd$CellType))
            stop("elements of argument 'cellTypes' is not part of ",
                 "'mSet$CellType'")
        keep <- which(pd$CellType %in% cellTypes)
        pd <- pd[keep,]
        p <- p[,keep]
    }
    # NOTE: Make cell type a factor
    pd$CellType <- factor(pd$CellType, levels = cellTypes)
    ffComp <- rowFtests(p, pd$CellType)
    prof <- vapply(
        X = splitit(pd$CellType),
        FUN = function(j) rowMeans2(p, cols = j),
        FUN.VALUE = numeric(nrow(p)))
    r <- rowRanges(p)
    compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
    names(compTable)[1] <- "Fstat"
    names(compTable)[c(-2, -1, 0) + ncol(compTable)] <-
        c("low", "high", "range")
    tIndexes <- splitit(pd$CellType)
    tstatList <- lapply(tIndexes, function(i) {
        x <- rep(0,ncol(p))
        x[i] <- 1
        return(rowttests(p, factor(x)))
    })

    if (probeSelect == "any") {
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[, "p.value"] < 1e-8, ]
            yAny <- y[order(abs(y[, "dm"]), decreasing = TRUE), ]
            c(rownames(yAny)[seq(numProbes * 2)])
        })
    } else {
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[, "p.value"] < 1e-8, ]
            yUp <- y[order(y[, "dm"], decreasing = TRUE), ]
            yDown <- y[order(y[, "dm"], decreasing = FALSE), ]
            c(rownames(yUp)[seq_len(numProbes)],
              rownames(yDown)[seq_len(numProbes)])
        })
    }

    trainingProbes <- unique(unlist(probeList))
    p <- p[trainingProbes,]

    pMeans <- colMeans2(p)
    names(pMeans) <- pd$CellType

    form <- as.formula(
        sprintf("y ~ %s - 1", paste(levels(pd$CellType), collapse = "+")))
    phenoDF <- as.data.frame(model.matrix(~ pd$CellType - 1))
    colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
    if (ncol(phenoDF) == 2) {
        # Two group solution
        X <- as.matrix(phenoDF)
        coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
    } else {
        # > 2 groups solution
        tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
        coefEsts <- tmp$coefEsts
    }

    list(
        coefEsts = coefEsts,
        compTable = compTable,
        sampleMeans = pMeans)
}

#==============================================================================
# Directories
#==============================================================================

# beta_dir <- USE 450K REFERENCE
probes_dir <- "~/shared-gandalm/brain_CTP/Data/methylation/reference/HM450.hg19.manifest_aut_mask.rds"
out_dir <- "~/shared-gandalm/brain_CTP/Data/reference_cell_profile/dlpfc_450k_guintivano"
filen <- "dlpfc_450k_guintivano_rowftest"

#==============================================================================
# Obtain coefficients
#==============================================================================

mask_rm <- rownames(readRDS(probes_dir))
dlpfc_exclude <- subsetByLoci(FlowSorted.DLPFC.450k, mask_rm)
dlpfc_mset <- preprocessIllumina(dlpfc_exclude)

compData <- pickCompProbes(
        mSet = dlpfc_mset,
        cellTypes = c("NeuN_pos", "NeuN_neg"),
        compositeCellType = "DPLFC",
        probeSelect = "any")

coefs <- compData$coefEsts

saveRDS(coefs, paste(out_dir, "/", filen, ".rds", sep = ""))
