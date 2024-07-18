###################################### Figure 7 ######################################

################## ------------------ Figure 7B ------------------ ##################
# 1. library
library(readr)
library(tidyverse)

# 2. params
dir_for_HALO_FCRL4 <- "data/FCRL4_induction_experiment.csv"
dir_for_result <- "figures"
celltype_color_panel <- c("#C3423F", "#5BC0EB", "#9BC53D")

# 3. load data
distance_df <- read_tsv(dir_for_HALO_FCRL4)
distance_df %>%
    group_by(from, to) %>%
    summarise(median_distance = median(distance)) %>%
    print()

# 4. plot
plot_df <- distance_df %>%
    dplyr::filter(from == "CD20+FCRL4+")
plot_df$to[plot_df$to == "CD4+ CD68-"] <- "CD4+ T"
plot_df$to[plot_df$to == "CD8+"] <- "CD8+ T"
plot_df$to[plot_df$to == "CD68+"] <- "CD68+ myeloid"
plot_df$to <- factor(plot_df$to, levels = rev(c("CD4+ T", "CD8+ T", "CD68+ myeloid")))

ggplot(plot_df, aes(x = to, y = distance, color = to, fill = to)) +
    geom_violin(scale = "width") +
    scale_color_manual(values = rev(celltype_color_panel)) +
    scale_fill_manual(values = rev(celltype_color_panel)) +
    stat_summary(
        fun = median, fun.min = median, fun.max = median,
        geom = "crossbar", width = 0.8, color = "black", size = 0.1
    ) +
    cowplot::theme_cowplot() +
    theme(
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        text = element_text(size = 7, family = "ArialMT"),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        plot.title = element_text(hjust = 0.5, size = 8),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
    ) +
    ggpubr::stat_compare_means(paired = FALSE, comparisons = list(c("CD4+ T", "CD8+ T"), c("CD4+ T", "CD68+ myeloid")), size = 0.35 * 8, tip.length = 0, label = "p.signif", label.y = c(50, 60)) +
    ylab("Closest distance (μm)") +
    xlab("") +
    ylim(0, 100) +
    coord_flip()

ggsave(file.path(dir_for_result, "7B.HALO_CD4_distance.pdf"),
    width = 2.7,
    height = 1.6
)

################## ------------------ Figure 7C ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(Seurat)

# 2. params
dir_for_B_obj <- "data/All_obj.rds"
dir_for_result <- "figures"
source("code/functions.R")

# 3. load data
B_obj <- read_rds(dir_for_B_obj)
B_obj_tumor <- subset(B_obj, subset = Tissue_short == "T")
Bm_tumor <- subset(B_obj_tumor, subset = Annotation_major == "Bm")

# 4. plot
gene_list <- c("LGMN", "CTSB", "IFI30", "CD74", "HLA-DMA", "HLA-DMB", "RFX5", "RFXANK", "NFYC", "CIITA")
ht <- gene_expression_heatmap(Bm_tumor, genes = gene_list, group.by = "Annotation", tile_size = 0.18)
pdf(file.path(dir_for_result, "7C.TAAB_APC_signatures.pdf"),
    height = 3.5,
    width = 4.5
)
draw(ht)
dev.off()

################## ------------------ Figures 7E and S7C ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(Seurat)
library(CellChat)

# 2. params
dir_for_CD45_obj <- "data/obj_CD45_20240708.rds"
dir_for_result <- "figures"
dir_for_data <- "data/CellChat"

# 3. cellchat function
run_cellchat <- function(expr_data,
                         label,
                         do_parrallel = TRUE,
                         workers = 10) {
    meta <- data.frame(
        label = label,
        stringsAsFactors = FALSE
    )
    require(CellChat)

    rownames(meta) <- colnames(expr_data)
    meta$label <- droplevels(meta$label, exclude = setdiff(levels(meta$label), unique(meta$label)))

    cellchat <- createCellChat(object = expr_data, meta = meta, group.by = "label")

    # Set the ligand-receptor interaction database
    CellChatDB <- CellChatDB.human
    # use all CellChatDB for cell-cell communication analysis
    CellChatDB.use <- CellChatDB # simply use the default CellChatDB
    # set the used database in the object
    cellchat@DB <- CellChatDB.use

    # subset the expression data of signaling genes for saving computation cost
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

    if (do_parrallel) {
        future::plan("multicore", workers = workers) # do parallel

        options(future.globals.maxSize = 5 * 1024^3)
    }

    print(paste(Sys.time(), "identifyOverExpressedGenes"))
    cellchat <- identifyOverExpressedGenes(cellchat)
    print(paste(Sys.time(), "identifyOverExpressedInteractions"))
    cellchat <- identifyOverExpressedInteractions(cellchat)

    # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
    print(paste(Sys.time(), "projectData"))
    cellchat <- projectData(cellchat, PPI.human)

    # cell-cell communication calculation
    # Compute the communication probability and infer cellular communication network
    print(paste(Sys.time(), "computeCommunProb"))
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    print(paste(Sys.time(), "filterCommunication"))
    cellchat <- filterCommunication(cellchat, min.cells = 5)

    # Infer the cell-cell communication at a signaling pathway level
    print(paste(Sys.time(), "computeCommunProbPathway"))
    cellchat <- computeCommunProbPathway(cellchat)

    # Calculate the aggregated cell-cell communication network
    print(paste(Sys.time(), "aggregateNet"))
    cellchat <- aggregateNet(cellchat)

    # Compute the network centrality scores
    print(paste(Sys.time(), "netAnalysis_computeCentrality"))
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

    return(cellchat)
}

# 4. run cellchat
CD45_obj <- read_rds(dir_for_CD45_obj)
CD45_tumor <- subset(CD45_obj, subset = Tissue_short == "T")

data <- list(
    expr_data = CD45_tumor@assays$RNA@data,
    cancer = CD45_tumor$Cancer,
    label = CD45_tumor$Annotation_CD45_2
)
data$label[data$label %in% c("Bn", "FCRL4- Bm", "FCRL4+ Bm", "Bgc", "PC_IGHA", "PC_IGHG", "ASC_other")] <- CD45_tumor$Annotation_CD45[data$label %in% c("Bn", "FCRL4- Bm", "FCRL4+ Bm", "Bgc", "PC_IGHA", "PC_IGHG", "ASC_other")]
data$label[startsWith(data$label, "CD4")] <- "CD4 T"
data$label[startsWith(data$label, "CD8")] <- "CD8 T"
data$label[startsWith(data$label, "cDC")] <- "cDC"
data$label <- factor(data$label,
    levels = c(
        "c01_Bn_TCL1A", "c02_Bn_NR4A2", "c03_Bn_IFN-response",
        "c04_classical-Bm_TXNIP", "c05_classical-Bm_GPR183", "c06_Bm_stress-response", "c07_Bm_IFN-response", "c08_ABC_FCRL4", "c09_ABC_FGR", "c10_Bm_TCL1A", "c11_pre-GC",
        "c12_Bgc_LZ-like",
        "c13_Bgc_DZ-like", "c14_Bm_activated-cycling", "c15_cycling_ASC",
        "c16_PC_IGHG", "c17_PC_IGHA", "c18_early-PC_MS4A1low", "c19_early-PC_LTB", "c20_early-PC_RGS13",
        "CD4 T",
        "CD8 T", "NK",
        "cDC", "pDC", "monoMacro", "Mast"
    )
)

rm(CD45_obj)
rm(CD45_tumor)
gc()

cellchat <- run_cellchat(data$expr_data,
    data$label,
    do_parrallel = TRUE,
    workers = 10
)
write_rds(cellchat, file.path(dir_for_data, "pan_cancer_CD45_cellchat_result.rds"), compress = "gz")

# 5. LR count lollipop plot
B_celltypes <- c(
    "c01_Bn_TCL1A", "c02_Bn_NR4A2", "c03_Bn_IFN-response",
    "c04_classical-Bm_TXNIP", "c05_classical-Bm_GPR183", "c06_Bm_stress-response", "c07_Bm_IFN-response", "c08_ABC_FCRL4", "c09_ABC_FGR", "c10_Bm_TCL1A", "c11_pre-GC",
    "c12_Bgc_LZ-like",
    "c13_Bgc_DZ-like", "c14_Bm_activated-cycling", "c15_cycling_ASC",
    "c16_PC_IGHG", "c17_PC_IGHA", "c18_early-PC_MS4A1low", "c19_early-PC_LTB", "c20_early-PC_RGS13"
)

## LR count with CD4 T cells
count <- list()
for (i_celltype_1 in B_celltypes) {
    i_celltype_2 <- "CD4 T"
    count <- append(count, list(data.frame(
        celltype_1 = i_celltype_1,
        celltype_2 = i_celltype_2,
        count = cellchat@net$count[i_celltype_1, i_celltype_2] + cellchat@net$count[i_celltype_2, i_celltype_1]
    )))
}
count <- bind_rows(count)

## plot
plot_df <- count %>%
    dplyr::filter(!str_detect(celltype_1, "^c1[345]")) # remove cycling clusters
plot_df$celltype_1 <- with(plot_df, reorder(celltype_1, count))
ggplot(plot_df, aes(x = celltype_1, y = count)) +
    geom_segment(aes(
        x = celltype_1,
        xend = celltype_1,
        y = 0,
        yend = count,
    ), color = "grey", size = 2) +
    geom_point(aes(color = count), size = 4.5) +
    scale_color_distiller(name = "Interaction count", palette = "YlOrRd", direction = 1) +
    geom_text(aes(label = count), color = "black", size = 7 * 0.35) +
    theme_light() +
    theme(
        text = element_text(size = 7, family = "ArialMT"),
        axis.text.x = element_text(size = 7, colour = "black", family = "ArialMT"),
        axis.text.y = element_text(size = 7, colour = "black", family = "ArialMT"),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        # legend.position = "none",
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        legend.key.width = unit(0.2, "inch"),
        legend.key.height = unit(0.15, "inch"),
        legend.title = element_text(size = 8)
    ) +
    coord_flip() +
    xlab("") +
    ylab("")

ggsave(file.path(dir_for_result, "S7C.B_and_CD4T_LR_count_lollipop.pdf"),
    height = 3.7,
    width = 4
)

# 6. LR pairs between TAABs and CD4 T cells
Bm_celltypes <- c("c04_classical-Bm_TXNIP", "c05_classical-Bm_GPR183", "c06_Bm_stress-response", "c07_Bm_IFN-response", "c08_ABC_FCRL4", "c09_ABC_FGR", "c10_Bm_TCL1A", "c11_pre-GC")

## ligand from B, receptor from T
lr_1 <- netVisual_bubble(cellchat,
    targets.use = "CD4 T",
    sources.use = Bm_celltypes,
    remove.isolate = FALSE,
    return.data = TRUE
)[[1]]
### filter differential LR pair
lr_statistics <- lr_1 %>%
    group_by(interaction_name) %>%
    summarise(n = n()) %>%
    dplyr::filter(n < n_distinct(lr_1$source))
lr_1 <- lr_1 %>%
    dplyr::filter(interaction_name %in% lr_statistics$interaction_name)
### format
lr_1 <- lr_1 %>%
    mutate(
        with = source,
        interaction_new = paste(ligand, "->", receptor)
    )

## ligand from T, receptor from B
lr_2 <- netVisual_bubble(cellchat,
    sources.use = "CD4 T",
    targets.use = Bm_celltypes,
    remove.isolate = FALSE,
    return.data = TRUE
)[[1]]
### filter differential LR pair
lr_statistics <- lr_2 %>%
    group_by(interaction_name) %>%
    summarise(n = n()) %>%
    dplyr::filter(n < n_distinct(lr_2$target))
lr_2 <- lr_2 %>%
    dplyr::filter(interaction_name %in% lr_statistics$interaction_name)
### format
lr_2 <- lr_2 %>%
    mutate(
        with = target,
        interaction_new = paste(receptor, "<-", ligand)
    )

## combine
lr <- rbind(lr_1, lr_2)

## plot
plot_df <- lr %>%
    dplyr::select(with, interaction_new, prob, pval)

plot_df$with <- factor(plot_df$with, levels = Bm_celltypes)

plot_df$interaction_new[plot_df$interaction_new %in% c("CD99 -> CD99", "CD99 <- CD99")] <- "CD99 <-> CD99"
plot_df$interaction_new[plot_df$interaction_new %in% c("CD22 -> PTPRC", "CD22 <- PTPRC")] <- "CD22 <-> PTPRC"
plot_df <- plot_df %>% distinct()
plot_df$interaction_new <- factor(plot_df$interaction_new, levels = sort(unique(plot_df$interaction_new), decreasing = TRUE))

values <- c(1, 2, 3)
names(values) <- c("P > 0.05", "0.01 < P < 0.05", "P < 0.01")

ggplot(plot_df, aes(x = with, y = interaction_new)) +
    geom_point(aes(color = prob, size = pval)) +
    scale_color_distiller(
        palette = "YlOrRd", direction = 1,
        breaks = c(quantile(plot_df$prob, 0, na.rm = T), quantile(plot_df$prob, 1, na.rm = T)), labels = c("Min", "Max")
    ) +
    scale_radius(
        range = c(min(plot_df$pval), max(plot_df$pval)),
        breaks = sort(unique(plot_df$pval)), labels = names(values)[values %in%
            sort(unique(plot_df$pval))], name = "p-value"
    ) +
    cowplot::theme_cowplot() +
    guides(size = guide_legend(title = "P-value")) +
    guides(color = guide_colorbar(title = "Communication\nprobility")) +
    theme(
        axis.text.x = element_text(
            size = 7,
            angle = 45,
            hjust = 1
        ),
        axis.text.y = element_text(size = 7),
        text = element_text(size = 7),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        panel.background = element_rect(colour = "black", fill = "white")
    ) +
    labs(x = "", y = "") +
    ggtitle("Interaction with CD4 T cells")

ggsave(file.path(dir_for_result, "7E.LR_TAAB_and_CD4_T.pdf"),
    height = 4.5,
    width = 4.2
)

################## ------------------ Figure 7G ------------------ ##################
# 1. library
library(readr)
library(tidyverse)

# 2. params
dir_for_HALO_CD69 <- "data/CD4_CD69_pct.csv"
dir_for_result <- "figures"

# 3. load data
CD69_pct_df <- read_csv(dir_for_HALO_CD69)
plot_df <- CD69_pct_df %>%
    gather(`0-20um`:`100-200um`, key = "distance", value = "CD69_pct") %>%
    mutate(CD69_pct = CD69_pct * 100) %>%
    dplyr::filter(distance %in% c("0-20um", "50-200um"))
plot_df$cancer <- str_extract(plot_df$slideID, "^[^-]*(?=-)")

# 4. plot
plot_df$distance <- factor(as.character(plot_df$distance), levels = c("0-20um", "50-200um"))
ggplot(plot_df, aes(x = distance, y = CD69_pct)) +
    geom_boxplot(outlier.color = NA, lwd = 0.3) +
    geom_point(size = 0.5, aes(color = cancer)) +
    scale_color_manual(values = c("#FF7F00FF", "#6A3D9AFF", "#B2DF8AFF", "#6DCCDAFF")) +
    geom_line(aes(group = slideID), colour = "grey", linetype = "11", size = 0.3) +
    cowplot::theme_cowplot() +
    theme(
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        text = element_text(size = 7, family = "ArialMT"),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        plot.title = element_text(hjust = 0.5, size = 8),
        # legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
    ) +
    ggpubr::stat_compare_means(paired = TRUE, method = "t.test", label.y = 75, size = 0.35 * 6, label = "p.signif") +
    ylab("Percentage of CD69+ in CD4+ T cells (%)") +
    xlab("Distance to TAABs")

ggsave(file.path(dir_for_result, "7G.HALO_CD4_CD69_pct.pdf"),
    width = 1.8,
    height = 1.8
)

################## ------------------ Figures 7H and S7F ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(Seurat)
library(DropletUtils)

# 2. params
dir_for_B_obj <- "data/All_obj.rds"
dir_for_result <- "figures"
dir_for_data <- "data/CytoSig"
dir_for_CytoSig_output <- "data/CytoSig/02.output.Zscore"
Cancer_color_panel <- c(
    AML = "#5A7B8F",
    BCC = "#EE4C97",
    BRCA = "#3D806F",
    CRC = "#A08634",
    CTCL = "#F37C95",
    ESCA = "#608541",
    FTC = "#7D4E57",
    HCC = "#BC3C29",
    HNSCC = "#958056",
    ICC = "#9FAFA3",
    MB = "#6F99AD",
    MELA = "#0072B5",
    NPC = "#CFC59A",
    NSCLC = "#E18727",
    OV = "#FFDC91",
    PACA = "#718DAE",
    RC = "#7876B1",
    STAD = "#F9AC93",
    THCA = "#3E6086",
    UCEC = "#20854E",
    UM = "#7581AF",
    cSCC = "#4A7985"
)

# 3. load data
B_obj <- read_rds(dir_for_B_obj)
Bm_tumor <- subset(B_obj, subset = Tissue_short == "T" & Annotation_major == "Bm")

# 4. run CytoSig
DropletUtils::write10xCounts(
    x = Bm_tumor@assays$RNA@counts,
    path = file.path(dir_for_data, "01.cellranger_output"),
    version = "3"
)
dir.create(file.path(dir_for_data, "02.output"), recursive = TRUE)
system("python3 CytoSig_run.py -i data/CytoSig/01.cellranger_output/ -o data/CytoSig/02.output -c 0 -z 1", wait = FALSE)

# 5. load CytoSig output
zscore_all <- data.table::fread(dir_for_CytoSig_output)
zscore_all <- zscore_all %>%
    column_to_rownames("V1") %>%
    t()

# 6. plot IL21 signaling in each Bm subset per cancer
cytokine <- "IL21"
plot_df <- data.frame(
    score = zscore_all[, cytokine],
    Annotation = Bm_tumor$Annotation[match(rownames(zscore_all), colnames(Bm_tumor))],
    Cancer = Bm_tumor$Cancer[match(rownames(zscore_all), colnames(Bm_tumor))]
)
plot_df <- plot_df %>%
    group_by(Cancer, Annotation) %>%
    summarise(
        mean_score = mean(score)
    )
plot_df$Annotation <- with(plot_df, reorder(Annotation, mean_score, function(x) -mean(x)))

my_comparisons <- list(
    c("c08_ABC_FCRL4", "c07_Bm_IFN-response"),
    c("c08_ABC_FCRL4", "c04_classical-Bm_TXNIP")
)

ggplot(plot_df, aes(x = Annotation, y = mean_score, color = Cancer)) +
    geom_jitter(size = 0.4, width = 0.2) +
    stat_summary(
        fun.data = mean_sdl, fun.args = list(mult = 1),
        geom = "pointrange", color = "black", fatten = 1.5
    ) +
    scale_color_manual(values = Cancer_color_panel) +
    cowplot::theme_cowplot() +
    theme(
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        text = element_text(size = 7, family = "ArialMT"),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        plot.title = element_text(hjust = 0.5, size = 8, face = "plain"),
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
    ) +
    xlab("") +
    ylab("Average activity") +
    ggtitle(paste0(cytokine, " signaling by CytoSig")) +
    ggpubr::stat_compare_means(comparisons = my_comparisons, label.y = c(2.9, 2.6), tip.length = 0, label = "p.signif", size = 2) +
    coord_cartesian(ylim = c(-3.4, 3.3))

ggsave(file.path(dir_for_result, "7H.IL21_signaling_CytoSig.pdf"),
    width = 3,
    height = 2.7
)

# 7. plot top 10 cytokines in c08
annotation_B <- Bm_tumor$Annotation[match(rownames(zscore_all), colnames(Bm_tumor))]
idx <- annotation_B == "c08_ABC_FCRL4"
mean_zscore <- colMeans(zscore_all[idx, ])
plot_df <- data.frame(
    cytokine = names(mean_zscore),
    activity = mean_zscore,
    stringsAsFactors = FALSE
)
plot_df <- plot_df %>%
    arrange(desc(activity)) %>%
    slice_head(n = 10)
plot_df$cytokine <- factor(plot_df$cytokine, levels = rev(plot_df$cytokine))

ggplot(plot_df, aes(x = cytokine, y = activity)) +
    geom_bar(stat = "identity", fill = "#6F99AD") +
    cowplot::theme_cowplot() +
    theme(
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        text = element_text(size = 7, family = "ArialMT"),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        plot.title = element_text(hjust = 0.5, size = 8, face = "plain"),
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
    ) +
    xlab("") +
    ylab("Average activity") +
    ggtitle("CytoSig top 10 cytokine signaling in TAAB") +
    ylim(-0.05, 1.55)

ggsave(file.path(dir_for_result, "S7F.TAAB_top10_cytokines.pdf"),
    width = 3,
    height = 2
)

# 8. plot top 10 cytokines in c04
annotation_B <- Bm_tumor$Annotation[match(rownames(zscore_all), colnames(Bm_tumor))]
idx <- annotation_B == "c04_classical-Bm_TXNIP"
mean_zscore <- colMeans(zscore_all[idx, ])
plot_df <- data.frame(
    cytokine = names(mean_zscore),
    activity = mean_zscore,
    stringsAsFactors = FALSE
)
plot_df <- plot_df %>%
    arrange(desc(activity)) %>%
    slice_head(n = 10)
plot_df$cytokine <- factor(plot_df$cytokine, levels = rev(plot_df$cytokine))

ggplot(plot_df, aes(x = cytokine, y = activity)) +
    geom_bar(stat = "identity", fill = "#6F99AD") +
    cowplot::theme_cowplot() +
    theme(
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        text = element_text(size = 7, family = "ArialMT"),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        plot.title = element_text(hjust = 0.5, size = 8, face = "plain"),
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
    ) +
    xlab("") +
    ylab("Average activity") +
    ggtitle("CytoSig top 10 cytokine signaling in c04_classical-Bm_TXNIP") +
    ylim(-0.05, 1.55)

ggsave(file.path(dir_for_result, "S7F.c04_top10_cytokines.pdf"),
    width = 3,
    height = 2
)

################## ------------------ Figure 7I ------------------ ##################
# 1. library
library(readr)
library(tidyverse)

# 2. params
dir_for_BCR <- "data/panB_BCR_20240709.rds"
dir_for_result <- "figures"

## palette
celltype_color_panel <- c(
    # Bm
    "c04_classical-Bm_TXNIP" = "#a82d06",
    "c05_classical-Bm_GPR183" = "#4592bf",
    "c06_Bm_stress-response" = "#d38219",
    "c07_Bm_IFN-response" = "#74a764",
    "c08_ABC_FCRL4" = "#8ca2b4",
    "c09_ABC_FGR" = "#cbb190",
    "c10_Bm_TCL1A" = "#e7ca8d",
    "c11_pre-GC" = "#9d9ec3"
)

# 3. load data
clone <- read_rds(dir_for_BCR)
clone_subset <- clone[clone$Tissue_short == "T", ]

# 4. clones containing ASCs
keep_clone_id <- clone_subset %>%
    group_by(clone_id) %>%
    summarise(target_cellType_count = sum(Annotation_major_2 == "ASC")) %>%
    dplyr::filter(target_cellType_count >= 1)
keep_clone_id <- keep_clone_id$clone_id

# 5. calculate percentage of each Bm subset belonging to these clones
plot_df <- clone_subset %>%
    dplyr::filter(Annotation_major == "Bm") %>%
    mutate(sharing_with_target_cellType = clone_id %in% keep_clone_id) %>%
    group_by(Annotation, SampleID) %>%
    summarise(
        sharing_count = sum(sharing_with_target_cellType),
        not_sharing_count = sum(!sharing_with_target_cellType),
        sharing_percent = 100 * sum(sharing_with_target_cellType) / n(),
        n = n()
    ) %>%
    dplyr::filter(n > 2) %>%
    arrange(desc(sharing_percent))

# 6. plot
my_comparisons <- list(
    c("c08_ABC_FCRL4", "c07_Bm_IFN-response"),
    c("c08_ABC_FCRL4", "c04_classical-Bm_TXNIP")
)
plot_df$Annotation <- with(plot_df, reorder(Annotation, sharing_percent, function(x) -median(x)))
ggplot(plot_df, aes(x = Annotation, y = sharing_percent, color = Annotation)) +
    geom_boxplot(outlier.color = NA, lwd = 0.3) +
    geom_jitter(size = 0.2, width = 0.2) +
    scale_color_manual(values = celltype_color_panel) +
    cowplot::theme_cowplot() +
    theme(
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        text = element_text(size = 7, family = "ArialMT"),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        plot.title = element_text(hjust = 0.5, size = 8),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
    ) +
    ggtitle("") +
    coord_cartesian(ylim = c(0, 100)) +
    ggpubr::stat_compare_means(comparisons = my_comparisons, label.y = c(80, 90), tip.length = 0, label = "p.signif", size = 2) +
    xlab("") +
    ylab("Proportion of cells sharing BCR with tumor ASCs")

ggsave(file.path(dir_for_result, "7I.Bm_BCR_sharing_with_ASCs.pdf"),
    width = 2.7,
    height = 2.5
)

################## ------------------ Figure 7J ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(Seurat)

# 2. params
dir_for_CD45_obj <- "data/obj_CD45_20240708.rds"
dir_for_result <- "figures"
source("code/functions.R")

# 3. load data
CD45_obj <- read_rds(dir_for_CD45_obj)
CD45_tumor <- subset(CD45_obj, subset = Tissue_short == "T")

# 4. plot
CD45_tumor$Annotation_CD45_2 <- factor(CD45_tumor$Annotation_CD45_2, levels = c(
    "Bn", "Bgc", "FCRL4- Bm", "FCRL4+ Bm", "PC_IGHA", "PC_IGHG", "ASC_other",
    "CD4_Tn", "CD4_Th_Other", "CD4_CXCL13", "CD4_Treg",
    "CD8_Tn", "CD8_Tex", "CD8_T_Other", "NK",
    "cDC_LAMP3", "cDC1", "cDC2", "pDC", "monoMacro", "Mast"
))
custom_dotplot(CD45_tumor, "IL21", group.by = "Annotation_CD45_2", coord_flip = TRUE)
ggsave(file.path(dir_for_result, "7J.IL21_expression.pdf"),
    width = 5.5,
    height = 1.5
)

################## ------------------ Figures 7L and S7J ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(AUCell)

# 2. params
dir_for_TAAB_signature <- "data/bulk_analysis/01.B_subset_signature_in_CD45/c08_ABC_FCRL4_wilcoxon_DE_genes.tsv"
dir_for_ICB_datasets <- "data/bulk_analysis/03.ICB_bulk_datasets"
dir_for_result <- "figures"

# 3. load TAAB signature
top_n_DE_genes <- 20
markers <- read_tsv(dir_for_TAAB_signature)
markers <- markers %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = top_n_DE_genes)
markers <- markers$gene

# 4. plot
aucMaxRank_perc <- 0.05

datasetIDs <- c(
    "NSCLC.Cho,J.-2020-ExpMolMed", "MELA.Gide,T.N.-2019-CancerCell", "MELA.Liu,D.-2019-NatMed", "MELA.VanAllen,E.M.-2015-Science",
    "MELA.Lauss,M.-2017-NatCommun", "MELA.Auslander,N.-2018-NatMed", "NSCLC.Prat,A.-2017-CancerRes"
)
for (datasetID in datasetIDs) {
    # load data
    data_obj <- read_rds(file.path(dir_for_ICB_datasets, paste0(datasetID, ".rds")))
    # subset
    if (datasetID %in% c("NSCLC.Cho,J.-2020-ExpMolMed", "MELA.Gide,T.N.-2019-CancerCell", "MELA.Liu,D.-2019-NatMed", "MELA.VanAllen,E.M.-2015-Science")) {
        subset <- data_obj$meta$biopsy_time == "pre-treatment" & data_obj$meta$responder %in% c("NR", "R")
    } else if (datasetID %in% c("MELA.Lauss,M.-2017-NatCommun", "MELA.Auslander,N.-2018-NatMed", "NSCLC.Prat,A.-2017-CancerRes")) {
        subset <- data_obj$meta$Treatment %in% c("PRE", "Pre", "pre") & data_obj$meta$response_NR %in% c("N", "R")
    }
    subset[is.na(subset)] <- FALSE
    tpm <- data_obj$tpm[, subset]
    meta <- data_obj$meta[subset, ]
    # AUCell scoring
    cells_rankings <- AUCell::AUCell_buildRankings(as.matrix(tpm), nCores = 10, plotStats = FALSE)
    cells_AUC <- AUCell::AUCell_calcAUC(list("1" = markers),
        cells_rankings,
        nCores = 10,
        aucMaxRank = ceiling(aucMaxRank_perc * nrow(cells_rankings))
    )
    abundance <- t(cells_AUC@assays@data$AUC)
    # plot
    if (datasetID %in% c("NSCLC.Cho,J.-2020-ExpMolMed", "MELA.Gide,T.N.-2019-CancerCell", "MELA.Liu,D.-2019-NatMed", "MELA.VanAllen,E.M.-2015-Science")) {
        plot_df <- data.frame(
            abundance = abundance,
            response = meta$responder
        )
    } else if (datasetID %in% c("MELA.Lauss,M.-2017-NatCommun", "MELA.Auslander,N.-2018-NatMed", "NSCLC.Prat,A.-2017-CancerRes")) {
        plot_df <- data.frame(
            abundance = abundance,
            response = meta$response_NR
        )
    }
    ggplot(plot_df, aes(x = response, y = abundance, color = response)) +
        geom_boxplot(outlier.shape = NA, lwd = 0.3) +
        geom_jitter(size = 0.4, width = 0.2) +
        scale_color_manual(values = c("#0173C1", "#EFBF02")) +
        cowplot::theme_cowplot() +
        theme(
            axis.text.x = element_text(size = 7),
            axis.text.y = element_text(size = 7),
            text = element_text(size = 7, family = "ArialMT"),
            plot.title = element_text(hjust = 0.5, size = 8, face = "plain"),
            legend.position = "none",
            axis.line = element_line(linetype = 1, color = "black", size = 0.3),
            axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
        ) +
        ggpubr::stat_compare_means(size = 0.35 * 7, label = "p.format") +
        ggtitle(paste0(datasetID, "\n", unique(meta$treatment))) +
        xlab("") +
        ylab("TAAB score")

    ggsave(file.path(dir_for_result, paste0("7L_S7J.ICB_", datasetID, ".pdf")),
        width = 1.5,
        height = 2
    )
}

################## ------------------ Figures 7K and S7H ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(Seurat)

# 2. params
dir_for_CD45_obj <- "data/obj_CD45_20240708.rds"
dir_for_result <- "figures"
dir_for_data <- "data"
source("code/functions.R")

# 3. load data
CD45_obj <- read_rds(dir_for_CD45_obj)

# 4. divide samples based on the median level of TAAB abundance
sample_df <- CD45_obj@meta.data %>%
    dplyr::filter(Tissue_short == "T") %>%
    group_by(Cancer, SampleID) %>%
    summarise(
        CD4_CXCL13 = sum(Annotation_CD45_2 == "CD4_CXCL13"),
        CD4 = sum(Annotation_CD45_2 %in% c("CD4_CXCL13", "CD4_Th_Other", "CD4_Tn", "CD4_Treg")),
        TAAB = sum(Annotation_CD45_2 == "FCRL4+ Bm"),
        B = sum(Annotation_CD45_major == "B"),
        n = n()
    )

sample_df$CD4_CXCL13_pct <- sample_df$CD4_CXCL13 / sample_df$CD4 * 100
sample_df$TAAB_pct <- sample_df$TAAB / sample_df$B * 100

sample_df$TAAB_group <- ifelse(sample_df$TAAB_pct > median(sample_df$TAAB_pct), "TAAB-high", "TAAB-low")
sample_df$TAAB_group <- factor(sample_df$TAAB_group, levels = c("TAAB-low", "TAAB-high"))

# 5. plot CD4_CXCL13 percentage in two groups
ggplot(sample_df, aes(x = TAAB_group, y = CD4_CXCL13_pct, color = TAAB_group)) +
    geom_boxplot(outlier.shape = NA, lwd = 0.3) +
    geom_jitter(size = 0.4, width = 0.2) +
    scale_color_manual(values = c("#0173C1", "#EFBF02")) +
    cowplot::theme_cowplot() +
    theme(
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        text = element_text(size = 7, family = "ArialMT"),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
    ) +
    ggpubr::stat_compare_means(size = 0.35 * 7, label = "p.format") +
    xlab("") +
    ylab("Percentage of CD4_CXCL13 cells in all CD4 T cells (%)")

ggsave(file.path(dir_for_result, "7K.CD4_CXCL13_percent_in_TAAB_high_vs_low_groups.pdf"),
    width = 1.2,
    height = 2.1
)

# 6. differential expression between the CD4 T cells from TAAB-high/low groups
CD4_obj <- subset(CD45_obj, subset = Tissue_short == "T" & Annotation_CD45_2 %in% c("CD4_CXCL13", "CD4_Th_Other", "CD4_Tn", "CD4_Treg"))

## divide sample groups
sample_df <- CD45_obj@meta.data %>%
    dplyr::filter(Tissue_short == "T") %>%
    group_by(Cancer, SampleID) %>%
    summarise(TAAB_pct = sum(Annotation_CD45_2 == "FCRL4+ Bm") / sum(Annotation_CD45_2 %in% c("FCRL4+ Bm", "FCRL4- Bm")))
sample_df$TAAB_group <- ifelse(sample_df$TAAB_pct > median(sample_df$TAAB_pct), "TAAB-high", "TAAB-low")
sample_df$TAAB_group <- factor(sample_df$TAAB_group, levels = c("TAAB-low", "TAAB-high"))
CD4_obj$TAAB_group <- sample_df$TAAB_group[match(CD4_obj$SampleID, sample_df$SampleID)]

## remove IG genes
CD4_obj <- subset(CD4_obj, features = rownames(CD4_obj)[!str_detect(rownames(CD4_obj), "^IG[HLK]")])
CD4_obj <- Seurat::NormalizeData(CD4_obj)

## DE
markers <- FindMarkers(CD4_obj, ident.1 = "TAAB-high", ident.2 = "TAAB-low", group.by = "TAAB_group", logfc.threshold = 0, test.use = "wilcox")
markers <- markers %>%
    rownames_to_column("gene") %>%
    arrange(desc(avg_log2FC))
write_tsv(markers, file.path(dir_for_data, "DEGs_in_CD4T_from_TAAB_high_vs_low_groups.tsv"))

## plot
FC_cutoff <- 0.3
p_cutoff <- 0.05

plot_df <- markers
plot_df$DE <- "Nonsignificant"
plot_df$DE[plot_df$avg_log2FC > FC_cutoff & plot_df$p_val_adj < p_cutoff] <- "TAAB-high group"
plot_df$DE[plot_df$avg_log2FC < -FC_cutoff & plot_df$p_val_adj < p_cutoff] <- "TAAB-low group"

plot_df$label <- NA
plot_df$label[plot_df$DE != "Nonsignificant"] <- plot_df$gene[plot_df$DE != "Nonsignificant"]

mycolors <- c("#0F7B9F", "#C3423F", "grey")
names(mycolors) <- c("TAAB-low group", "TAAB-high group", "Nonsignificant")

ggplot(data = plot_df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = DE, label = label)) +
    geom_point(size = 0.2) +
    scale_colour_manual(values = mycolors, name = "Enriched in tumor CD4 T from") +
    geom_vline(xintercept = c(FC_cutoff, -FC_cutoff), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black") +
    ggrepel::geom_text_repel(size = 7 * 0.35, show.legend = FALSE) +
    cowplot::theme_cowplot() +
    theme(
        text = element_text(size = 7, family = "ArialMT"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, size = 8),
        # legend.position = "none",
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
    ) +
    xlab("Log2(fold change)") +
    ylab("-Log10(q-value)")

ggsave(file.path(dir_for_result, "S7H.CD4_differential_expression.pdf"),
    width = 3.8,
    height = 2.3
)

ggsave(file.path(dir_for_result, "S7H.CD4_differential_expression_large.pdf"),
    width = 6,
    height = 5
)

################## ------------------ Figure S7K ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(Seurat)
library(harmony)
library(future)

library(Nebulosa)

# set multiple core
plan("multicore", workers = 10)
options(future.globals.maxSize = 50 * 1024^3)

# 2. params
dir_for_obj_Ma <- "data/Ma_et_al_reanalysis/panB_scRNA_processed_data.rds"
dir_for_black_list_genes <- "data/black_list_genes.tsv"
dir_for_result <- "figures"
dir_for_data <- "data/Ma_et_al_reanalysis"
source("code/functions.R")
source("code/harmony_functions.R")

# 3. load data
obj_Ma <- read_rds(dir_for_obj_Ma)
AtM_Ma <- subset(obj_Ma, subset = celltype == "B.09.DUSP4+AtM")

# 4. run Harmony
## remove batches with too few cells
AtM_Ma <- subset(AtM_Ma, subset = dataid %in% names(table(AtM_Ma$dataid)[table(AtM_Ma$dataid) > 50]))
## run Harmony
AtM_Ma <- yy_harmony_batch_correction(
    obj = AtM_Ma,
    path = dir_for_data,
    batchID = "dataid",
    find_commmon_hvg = TRUE,
    min_cell_in_batch = 50,
    nhvg = 2000,

    # black list genes
    use_black_list = TRUE,
    dir_for_black_list_genes = dir_for_black_list_genes,
    remove_tom_gene = TRUE,
    remove_yang_gene = TRUE,

    # scale data
    scale_seperately = TRUE,
    regress_vector = NULL,

    # harmony batchID
    try_theta = FALSE,
    harmony_batchID = "dataid",
    harmony_theta = 2,
    harmony_lambda = 1,

    # calculate DEGs
    calculate_deg = FALSE
)

# 5. load data after Harmony
obj <- read_rds(file.path(dir_for_data, "obj_after_harmony.rds"))

# 6. plot tissue distribution
obj$type <- factor(obj$type, levels = c("Cancer", "Adjacent", "Blood", "LN_Met"))
custom_dimplot(obj = obj, group.by = "type", pt.size = 0.3, label = FALSE) +
    scale_color_manual(values = c(
        "Cancer" = "#5BC0EB", "Adjacent" = "#9BC53D",
        "Blood" = "#C3423F", "LN_Met" = "#FDE74C"
    ))
ggsave(file.path(dir_for_result, "S7K.Ma_et_al_AtM_tissue_distribution.pdf"), width = 4, height = 2.8)

# 7. plot FCRL4 and FGR gene expression
Nebulosa::plot_density(obj, features = c("FCRL4", "FGR"), size = 0.3) &
    # scale_color_viridis_c(limits = c(0, 0.07), oob = scales::squish) &
    theme_bw() &
    theme(
        aspect.ratio = 1,
        plot.margin = margin(t = 0.2, r = 0, b = 0, l = 0),
        plot.title = element_text(size = 7, hjust = 0.5, vjust = -3),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()
    )
ggsave(file.path(dir_for_result, "S7K.Ma_et_al_AtM_FCRL4_FGR_expression.pdf"), width = 5, height = 3)
