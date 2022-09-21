#==============================================================================
#
# Find marker probes for sequencing data using "extremes" method
#
#==============================================================================

#==============================================================================
# 1. Set-up
#==============================================================================

# Libraries
library(data.table)
library(tidyverse)

# Directories
cgn_dmr_ilmn450k_all_methylC.ma_dir <- "~/shared-gandalm/brain_CTP/Data/reference_cell_profile/Luo2020/ilmn450kepic_aggto100bp_celltype_c10_M.txt"
cgn_dmr_ilmn450k_all_methylunmethylC.ma_dir <- "~/shared-gandalm/brain_CTP/Data/reference_cell_profile/Luo2020/ilmn450kepic_aggto100bp_celltype_c10_MU.txt"

# Read in data
cgn_dmr_ilmn450k_all_methylC <- fread(cgn_dmr_ilmn450k_all_methylC.ma_dir)
cgn_dmr_ilmn450k_all_methylunmethylC <- fread(cgn_dmr_ilmn450k_all_methylunmethylC.ma_dir)

# Convert to matrix format
cgn_dmr_ilmn450k_all_methylC.ma <- as.matrix(cgn_dmr_ilmn450k_all_methylC[,2:ncol(cgn_dmr_ilmn450k_all_methylC)])
rownames(cgn_dmr_ilmn450k_all_methylC.ma) <- cgn_dmr_ilmn450k_all_methylC$Probe
cgn_dmr_ilmn450k_all_methylunmethylC.ma <- as.matrix(cgn_dmr_ilmn450k_all_methylunmethylC[,2:ncol(cgn_dmr_ilmn450k_all_methylunmethylC)])
rownames(cgn_dmr_ilmn450k_all_methylunmethylC.ma) <- cgn_dmr_ilmn450k_all_methylunmethylC$Probe

# Extract common cell-types
cgn_dmr_ilmn450k_all_common_methylC.ma <- cgn_dmr_ilmn450k_all_methylC.ma
cgn_dmr_ilmn450k_all_common_methylunmethylC.ma <- cgn_dmr_ilmn450k_all_methylunmethylC.ma

# Calculate U values
cgn_dmr_ilmn450k_all_common_unmethylC.ma <- cgn_dmr_ilmn450k_all_common_methylunmethylC.ma - cgn_dmr_ilmn450k_all_common_methylC.ma

#==============================================================================
# 2. QC of data
#==============================================================================

#------------------------------------------------------------------------------
# a) Grouping into excitatory, inhibitory, astrocytes, microglial, endo, oligo +/- fetal
#------------------------------------------------------------------------------

# - sum across subtypes
exc_col <- c("Exc_L1-3_CUX2", "Exc_L2-4_RORB", "Exc_L4-5_FOXP2", "Exc_L4-6_LRRK1", "Exc_L4_PLCH1", "Exc_L4-5_TOX", "Exc_L5-6_PDZRN4", "Exc_L6_TLE4", "Exc_L6_TSHZ2")
inh_col <- c("Inh_CGE_LAMP5", "Inh_CGE_NDNF", "Inh_CGE_VIP", "Inh_CGE-MGE_CHST9", "Inh_MGE_B3GAT2", "Inh_MGE_CALB1", "Inh_MGE_PVALB", "Inh_MGE_UNC5B")
keep_col <- c("NonN_Astro_FGF3R", "NonN_Endo", "NonN_Micro", "NonN_Oligo_MBP", "NonN_OPC")

# exc_col <- c("Exc_L1-3_CUX2", "Exc_L2-4_RORB", "Exc_L4-5_FOXP2", "Exc_L4-6_LRRK1", "Exc_L4_PLCH1", "Exc_L5-6_PDZRN4", "Exc_L6_TLE4")
# inh_col <- c("Inh_CGE_LAMP5", "Inh_CGE_NDNF", "Inh_CGE_VIP", "Inh_MGE_B3GAT2", "Inh_MGE_CALB1", "Inh_MGE_PVALB")
# keep_col <- c("NonN_Astro_FGF3R", "NonN_Micro-Endo_TYROBP", "NonN_Oligo_MBP", "hs_fc_UMB_412")

# - methylC
cgn_dmr_ilmn450k_all_common_methylC_ctgroup.ma <- cbind(rowSums(cgn_dmr_ilmn450k_all_common_methylC.ma[,exc_col], na.rm = T), rowSums(cgn_dmr_ilmn450k_all_common_methylC.ma[,inh_col], na.rm = T), cgn_dmr_ilmn450k_all_common_methylC.ma[,keep_col])
colnames(cgn_dmr_ilmn450k_all_common_methylC_ctgroup.ma)[1:2] <- c("Exc", "Inh")

# - unmethylC
cgn_dmr_ilmn450k_all_common_unmethylC_ctgroup.ma <- cbind(rowSums(cgn_dmr_ilmn450k_all_common_unmethylC.ma[,exc_col], na.rm = T), rowSums(cgn_dmr_ilmn450k_all_common_unmethylC.ma[,inh_col], na.rm = T), cgn_dmr_ilmn450k_all_common_unmethylC.ma[,keep_col])
colnames(cgn_dmr_ilmn450k_all_common_unmethylC_ctgroup.ma)[1:2] <- c("Exc", "Inh")

# - methylunmethylC
cgn_dmr_ilmn450k_all_common_methylunmethylC_ctgroup.ma <- cbind(rowSums(cgn_dmr_ilmn450k_all_common_methylunmethylC.ma[,exc_col], na.rm = T), rowSums(cgn_dmr_ilmn450k_all_common_methylunmethylC.ma[,inh_col], na.rm = T), cgn_dmr_ilmn450k_all_common_methylunmethylC.ma[,keep_col])
colnames(cgn_dmr_ilmn450k_all_common_methylunmethylC_ctgroup.ma)[1:2] <- c("Exc", "Inh")

# - check NAs in methylC file
na_check_probe <- rowSums(cgn_dmr_ilmn450k_all_common_methylC_ctgroup.ma)
length(which(!is.na(na_check_probe)))
na_check_ctp <- apply(cgn_dmr_ilmn450k_all_common_methylC_ctgroup.ma, 2, function(x) length(which(is.na(x))))
na_check_ctp

# - remove fetal to improve missingness
cgn_dmr_ilmn450k_all_common_methylC_ctgroup_rmfetal.ma <- cgn_dmr_ilmn450k_all_common_methylC_ctgroup.ma[,c("Exc", "Inh", keep_col)]
cgn_dmr_ilmn450k_all_common_unmethylC_ctgroup_rmfetal.ma <- cgn_dmr_ilmn450k_all_common_unmethylC_ctgroup.ma[,c("Exc", "Inh", keep_col)]
cgn_dmr_ilmn450k_all_common_methylunmethylC_ctgroup_rmfetal.ma <- cgn_dmr_ilmn450k_all_common_methylunmethylC_ctgroup.ma[,c("Exc", "Inh", keep_col)]

na_check_probe_rmfetal <- rowSums(cgn_dmr_ilmn450k_all_common_methylC_ctgroup_rmfetal.ma[,c("Exc", "Inh", keep_col)])
length(which(!is.na(na_check_probe_rmfetal)))

# - keep NA0 only
na_rm <- which(is.na(na_check_probe_rmfetal))
cgn_dmr_ilmn450k_all_common_methylC_ctgroup_rmfetal_na0.ma <- cgn_dmr_ilmn450k_all_common_methylC_ctgroup_rmfetal.ma[-na_rm,]
cgn_dmr_ilmn450k_all_common_unmethylC_ctgroup_rmfetal_na0.ma <- cgn_dmr_ilmn450k_all_common_unmethylC_ctgroup_rmfetal.ma[-na_rm,]
cgn_dmr_ilmn450k_all_common_methylunmethylC_ctgroup_rmfetal_na0.ma <- cgn_dmr_ilmn450k_all_common_methylunmethylC_ctgroup_rmfetal.ma[-na_rm,]
# - remove rows where methylunmethylC = 0, as this is equivalent to 0 when adding with na.rm = T command
na_sum_rm <- which(cgn_dmr_ilmn450k_all_common_methylunmethylC_ctgroup_rmfetal_na0.ma == 0, arr.ind = T)
cgn_dmr_ilmn450k_all_common_methylC_ctgroup_rmfetal_na0.ma <- cgn_dmr_ilmn450k_all_common_methylC_ctgroup_rmfetal_na0.ma[-na_sum_rm[,1],]
cgn_dmr_ilmn450k_all_common_unmethylC_ctgroup_rmfetal_na0.ma <- cgn_dmr_ilmn450k_all_common_unmethylC_ctgroup_rmfetal_na0.ma[-na_sum_rm[,1],]
cgn_dmr_ilmn450k_all_common_methylunmethylC_ctgroup_rmfetal_na0.ma <- cgn_dmr_ilmn450k_all_common_methylunmethylC_ctgroup_rmfetal_na0.ma[-na_sum_rm[,1],]

#==============================================================================
# 3. Identify marker probes
#==============================================================================

methyl_pc <- cgn_dmr_ilmn450k_all_common_methylC_ctgroup_rmfetal_na0.ma/cgn_dmr_ilmn450k_all_common_methylunmethylC_ctgroup_rmfetal_na0.ma

# Find how many have %>60 and %<40
ndiff <- (ncol(methyl_pc)-1)
methyl_pc_diff <- data.frame(methyl_pc)
methyl_pc_diff$high <- apply(data.frame(methyl_pc), 1, function(x) length(which(x >= 0.6)))
methyl_pc_diff$low <- apply(data.frame(methyl_pc), 1, function(x) length(which(x <= 0.4)))
methyl_pc_diff_marker <- methyl_pc_diff %>% filter((high == ndiff & low == 1) | (high == 1 & low == ndiff))
methyl_pc_diff_marker$direction <- ifelse(methyl_pc_diff_marker$high == (ncol(methyl_pc)-1), "down", "up")
methyl_pc_diff_marker.tmp <- data.frame(Probe = rownames(methyl_pc_diff_marker), methyl_pc_diff_marker[,-which(colnames(methyl_pc_diff_marker) %in% c("high", "low"))])
methyl_pc_diff_marker.long <- melt(methyl_pc_diff_marker.tmp, id.vars = c("Probe", "direction"), variable.name = "celltype", value.name = "prop")
methyl_pc_diff_marker.long <- methyl_pc_diff_marker.long %>% mutate(marker = ifelse((direction == "down" & prop <= 0.4) | (direction == "up" & prop >= 0.6), as.character(celltype), NA))
methyl_pc_diff_marker.long <- methyl_pc_diff_marker.long %>% 
    mutate(binary = case_when((!is.na(marker) & direction == "down") ~ 0,
                                (!is.na(marker) & direction == "up") ~ 1,
                                (is.na(marker) & direction == "down") ~ 1,
                                (is.na(marker) & direction == "up") ~ 0))

methyl_pc_diff_marker.long.out <- methyl_pc_diff_marker.long %>% filter(!is.na(marker))
saveRDS(methyl_pc_diff_marker.long.out, "~/shared-gandalm/brain_CTP/Data/reference_cell_profile/Luo2020/Luo2020_extremes_dmr_ilmn450kepic_aggto100bp_c10_cell7_split6040all_markerprobes.rds")

methyl_pc_diff_marker.beta <- as.matrix(methyl_pc_diff_marker.tmp[,-which(colnames(methyl_pc_diff_marker.tmp) %in% c("Probe", "direction"))])
saveRDS(methyl_pc_diff_marker.beta, "~/shared-gandalm/brain_CTP/Data/reference_cell_profile/Luo2020/Luo2020_extremes_dmr_ilmn450kepic_aggto100bp_c10_cell7_split6040all_beta.rds")

cell_type_markers_extremes_6040 <- unique(methyl_pc_diff_marker.long$Probe)

# - take the U, M, U+M matrices corresponding, and sort rownames to correspond with cell_type_markers_chisq_top200
M_out <- cgn_dmr_ilmn450k_all_common_methylC_ctgroup_rmfetal_na0.ma[which(rownames(cgn_dmr_ilmn450k_all_common_methylC_ctgroup_rmfetal_na0.ma) %in% cell_type_markers_extremes_6040),]
M_out <- M_out[order(match(rownames(M_out), cell_type_markers_extremes_6040)),]
U_out <- cgn_dmr_ilmn450k_all_common_unmethylC_ctgroup_rmfetal_na0.ma[which(rownames(cgn_dmr_ilmn450k_all_common_unmethylC_ctgroup_rmfetal_na0.ma) %in% cell_type_markers_extremes_6040),]
U_out <- U_out[order(match(rownames(U_out), cell_type_markers_extremes_6040)),]
MU_out <- cgn_dmr_ilmn450k_all_common_methylunmethylC_ctgroup_rmfetal_na0.ma[which(rownames(cgn_dmr_ilmn450k_all_common_methylunmethylC_ctgroup_rmfetal_na0.ma) %in% cell_type_markers_extremes_6040),]
MU_out <- MU_out[order(match(rownames(MU_out), cell_type_markers_extremes_6040)),]

identical(rownames(M_out), cell_type_markers_extremes_6040)
identical(rownames(MU_out), cell_type_markers_extremes_6040)

# Write output
saveRDS(M_out, "~/shared-gandalm/brain_CTP/Data/reference_cell_profile/Luo2020/ilmn450kepic_aggto100bp_c10_extremes_split6040all_celltype7_M.rds")
saveRDS(U_out, "~/shared-gandalm/brain_CTP/Data/reference_cell_profile/Luo2020/ilmn450kepic_aggto100bp_c10_extremes_split6040all_celltype7_U.rds")
saveRDS(MU_out, "~/shared-gandalm/brain_CTP/Data/reference_cell_profile/Luo2020/ilmn450kepic_aggto100bp_c10_extremes_split6040all_celltype7_MU.rds")

#==============================================================================
# 4. Heatmap
#==============================================================================

M_MU_out <- M_out/MU_out
#zmat <- readRDS("~/shared-gandalm/brain_CTP/Data/reference_cell_profile/Luo2020/Luo2020_chisq_dmr_cell5_bonf001_contr50pc_top200.rds")$zmat
cell_type_markers_extreme_markers_stats <- readRDS("~/shared-gandalm/brain_CTP/Data/reference_cell_profile/Luo2020/Luo2020_extremes_dmr_ilmn450kepic_aggto100bp_c10_cell7_split6040all_markerprobes.rds")
#cell_type_markers_chisq_top200_stats <- readRDS("~/shared-gandalm/brain_CTP/Data/reference_cell_profile/Luo2020/Luo2020_chisq_dmr_cell7_bonf001_contr50pc_top200_stats.rds")
M_MU_out_order <- M_MU_out[order(match(rownames(M_MU_out), cell_type_markers_extreme_markers_stats$Probe)),]
identical(rownames(M_MU_out_order), cell_type_markers_extreme_markers_stats$Probe)
M_MU_out_order_t <- t(M_MU_out_order)
rownames(M_MU_out_order_t) <- ifelse(rownames(M_MU_out_order_t) %in% c("NonN_Astro_FGF3R", "NonN_Endo", "NonN_Micro", "NonN_Oligo_MBP", "NonN_OPC"),
    unlist(lapply(strsplit(rownames(M_MU_out_order_t), "_"), function(x) x[2])), rownames(M_MU_out_order_t))
celltype <- gsub("NonN_", "", cell_type_markers_extreme_markers_stats$marker)
celltype <- unlist(lapply(strsplit(celltype, "_"), function(x) x[1]))
direction <- cell_type_markers_extreme_markers_stats$direction
#png("~/shared-gandalm/brain_CTP/Data/reference_cell_profile/heatmap_Luo2020_chisq_dmr_cell7_bonf001_contr50pc_top200_MdivMUprop_probedendro_cellcluster.png", height = 750, width = 1000)
png("~/shared-gandalm/brain_CTP/Data/reference_cell_profile/heatmap_LLuo2020_extremes_dmr_ilmn450kepic_aggto100bp_c10_cell7_split6040all_MdivMUprop_probedendro_cellcluster.png", height = 750, width = 1000)
ha <- HeatmapAnnotation(celltype = celltype, direction = direction)
Heatmap(M_MU_out_order_t, name = "mat", top_annotation = ha, show_row_names = TRUE, show_column_names = FALSE, column_split = celltype, cluster_rows = FALSE)
dev.off()

table(methyl_pc_diff_marker.long$marker)
#
#             Exc              Inh NonN_Astro_FGF3R        NonN_Endo 
#             983               99              499              682 
#      NonN_Micro   NonN_Oligo_MBP         NonN_OPC 
#             763              423              838 

#==============================================================================
# 5. Get into format for deconvolution with essential parts:
#==============================================================================

# - regions_all: GRanges object (I think you only really need the positional info)
# - zmat: 1/0 matrix per CpG, where 1 is up-methylated vs 0 is down-methylated

# regions_all
ilmn450k.gr_dir <- "~/shared-gandalm/brain_CTP/Data/methylation/reference/HM450.hg19.manifest.rds"
ilmn450k.gr <- readRDS(ilmn450k.gr_dir)
mcols(ilmn450k.gr) <- NULL

#------------------------------------------------------------------------------
# a) Collapsed
#------------------------------------------------------------------------------

# zmat
zmat.tmp <- dcast(methyl_pc_diff_marker_addNoise.long, Probe + direction ~ celltype, value.var = "binary") %>%
    mutate(Noise = case_when((is.na(Noise) & direction == "down") ~ 1,
                                (is.na(Noise) & direction == "up") ~ 0,
                                (!is.na(Noise) & direction == "down") ~ 0,
                                (!is.na(Noise) & direction == "up") ~ 1))
zmat.tmp <- dcast(methyl_pc_diff_marker_addGABA.long, Probe + direction ~ celltype, value.var = "binary")
zmat <- zmat.tmp %>% dplyr::select(-c("Probe", "direction"))
rownames(zmat) <- zmat.tmp$Probe

# regions_all
regions_all <- ilmn450k.gr[which(names(ilmn450k.gr) %in% zmat.tmp$Probe),]
# - sort to be in order
regions_all <- regions_all[order(match(names(regions_all), zmat.tmp$Probe)),]
identical(names(regions_all), zmat.tmp$Probe)

out <- list(regions_all = regions_all, zmat = zmat)
saveRDS(out, "~/shared-gandalm/brain_CTP/Data/reference_cell_profile/Luo2020/Luo2020_extremes_dmr_ilmn450kepic_aggto100bp_c10_cell7_split6040all.rds")
