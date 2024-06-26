###################################### Figure 7 ######################################

################## ------------------ Figure 7B ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(Seurat)

# 2. params
dir_for_HALO_FCRL4 <- ""
dir_for_result <- ""
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

ggsave(file.path(dir_for_result_part, "7B.HALO_CD4_distance.pdf"),
    width = 2.7,
    height = 1.6
)

################## ------------------ Figure 7C ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(Seurat)

# 2. params
dir_for_B_obj <- ""
dir_for_result <- ""
source("./functions.R")

# 3. load data
B_obj <- read_rds(dir_for_B_obj)
B_obj_tumor <- subset(B_obj, subset = Tissue_short == "T")
Bm_tumor <- subset(B_obj_tumor, subset = Annotation_major == "Bm")

# 4. plot
gene_list <- c("LGMN", "CTSB", "IFI30", "CD74", "HLA-DMA", "HLA-DMB", "RFX5", "RFXANK", "NFYC", "CIITA")
ht <- gene_expression_heatmap(Bm_tumor, genes = gene_list, group.by = "annotation", tile_size = 0.18)
pdf(file.path(dir_for_result, "7C.TAAB_APC_signatures.pdf"),
    height = 3.5,
    width = 4.5
)
draw(ht)
dev.off()
