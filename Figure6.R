###################################### Figure 6 ######################################

################## ------------------ Figure 6A ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(Seurat)

# 2. params
dir_for_B_obj <- ""
dir_for_celltypist_result <- ""
dir_for_result <- ""
Cancer_color_panel <- read_rds("")

# 3. load data
B_obj <- read_rds(dir_for_B_obj)
B_obj_tumor <- subset(B_obj, subset = Tissue_short == "T")

# add celltypist prediction of cycling Bm cells
celltypist_prediction <- read.csv(dir_for_celltypist_result)
celltypist_prediction$celltypist_prediction <- ifelse(celltypist_prediction$predicted_labels == 7, "FCRL4+", "FCRL4-")
B_obj_tumor$celltypist_prediction <- NA
B_obj_tumor$celltypist_prediction[B_obj_tumor$Annotation_major == "Bcycling"] <- celltypist_prediction$celltypist_prediction[match(B_obj_tumor$CellID_final[B_obj_tumor$Annotation_major == "Bcycling"], celltypist_prediction$CellID_final)]
B_obj_tumor$Annotation_major_3 <- as.character(B_obj_tumor$Annotation_major_2)
B_obj_tumor$Annotation_major_3[B_obj_tumor$Annotation_major_2 == "Bm"] <- "FCRL4- Bm"
B_obj_tumor$Annotation_major_3[B_obj_tumor$Annotation == "c08_ABC_FCRL4"] <- "FCRL4+ Bm"
B_obj_tumor$Annotation_major_3[B_obj_tumor$Annotation == "c14_Bm_activated−cycling" & B_obj_tumor$celltypist_prediction == "FCRL4+"] <- "FCRL4+ Bm"

# 4. plot
plot_df <- B_obj_tumor@meta.data %>%
    dplyr::filter(treatment_status %in% c("pre", "Pre", "Pre-treatment")) %>%
    group_by(SampleID, Cancer) %>%
    summarise(
        FCRL4 = sum(Annotation_major_3 == "FCRL4+ Bm"),
        n = n(),
        FCRL4_pct = FCRL4 / n * 100
    ) %>%
    dplyr::filter(n > 50)

plot_df$Cancer <- with(plot_df, reorder(Cancer, FCRL4_pct, median))

ggplot(plot_df, aes(x = Cancer, y = FCRL4_pct, color = Cancer)) +
    geom_boxplot(outlier.shape = NA, lwd = 0.4, fatten = 1) +
    geom_jitter(shape = 16, size = 0.3, stroke = 0.5, width = 0.2) +
    scale_color_manual(values = Cancer_color_panel) +
    cowplot::theme_cowplot() +
    theme(
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, size = 8, family = "ArialMT", face = "plain"),
        text = element_text(size = 7, family = "ArialMT"),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
    ) +
    xlab("") +
    ylab("% among TIBs") +
    ggtitle("FCRL4+ Bm abundance") +
    coord_flip(ylim = c(0, 25))

ggsave(file.path(dir_for_result, "6A.FCLR4_Bm_abundance.pdf"),
    device = "pdf",
    width = 1.8,
    height = 2.5
)

################## ------------------ Figure 6B ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(Matrix)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)

# 2. params
dir_for_peak_counts <- ".../filtered_peak_bc_matrix.h5"
dir_for_fragments <- ".../fragments.tsv.gz"
dir_for_meta <- ".../singlecell.csv"
dir_for_B_obj <- ""
dir_for_result <- ""

# 3. create Seurat object
datasetID <- "BC.Zhang,Y.-2021-CancerCell.10X5"
counts <- Seurat::Read10X_h5(dir_for_peak_counts)
chrom_assay <- Signac::CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = "hg38",
    fragments = dir_for_fragments,
    min.cells = 10,
    min.features = 200
)
sample_ref <- c(
    "1" = "Pre_P012_t",
    "2" = "Post_P012_t",
    "3" = "Pre_P020_t",
    "4" = "Pre_P023_t",
    "5" = "Prog_P013_t"
)
metadata <- data.frame(
    barcode = colnames(chrom_assay),
    datasetID = datasetID,
    sampleID = sample_ref[str_extract(colnames(chrom_assay), "(?<=-).*$")],
    stringsAsFactors = FALSE
)
metadata$treatment_status <- str_extract(metadata$sampleID, "^[A-Za-z]*(?=_)")
metadata_tmp <- read_csv(dir_for_meta)
metadata <- metadata %>%
    left_join(metadata_tmp %>% dplyr::select(barcode:duplicate), by = "barcode")
rownames(metadata) <- metadata$barcode
ATAC_obj <- SeuratObject::CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
)

# 4. add gene annotations
## extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
## change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
## add the gene information to the object
Annotation(ATAC_obj) <- annotations

# 5. QC
## compute nucleosome signal score per cell
ATAC_obj <- NucleosomeSignal(object = ATAC_obj)
## compute TSS enrichment score per cell
ATAC_obj <- TSSEnrichment(object = ATAC_obj, fast = FALSE)
## add blacklist ratio and fraction of reads in peaks
ATAC_obj$pct_reads_in_peaks <- ATAC_obj$peak_region_fragments / ATAC_obj$passed_filters * 100
ATAC_obj$blacklist_ratio <- ATAC_obj$blacklist_region_fragments / ATAC_obj$peak_region_fragments
## filter
ATAC_obj <- subset(
    x = ATAC_obj,
    subset = peak_region_fragments > 3000 &
        peak_region_fragments < 20000 &
        pct_reads_in_peaks > 15 &
        blacklist_ratio < 0.05 &
        nucleosome_signal < 4 &
        TSS.enrichment > 2
)

# 6. normalization and linear dimensional reduction
ATAC_obj <- RunTFIDF(ATAC_obj)
ATAC_obj <- FindTopFeatures(ATAC_obj, min.cutoff = "q75")
ATAC_obj <- RunSVD(ATAC_obj, features = VariableFeatures(ATAC_obj))

# 7. non-linear dimension reduction and clustering
ATAC_obj <- RunUMAP(object = ATAC_obj, reduction = "lsi", dims = 2:30)
ATAC_obj <- RunTSNE(object = ATAC_obj, reduction = "lsi", dims = 2:30)
ATAC_obj <- FindNeighbors(object = ATAC_obj, reduction = "lsi", dims = 2:30)
ATAC_obj <- FindClusters(object = ATAC_obj, verbose = FALSE, algorithm = 3, resolution = 0.8)

# 8. gene activity
plan("multisession", workers = 20)
plan()
options(future.globals.maxSize = 10 * 1024^3) # for 10 Gb RAM
gene.activities <- GeneActivity(ATAC_obj)
plan("sequential")
## add the gene activity matrix to the Seurat object as a new assay and normalize it
ATAC_obj[["RNA"]] <- CreateAssayObject(counts = gene.activities)
ATAC_obj <- NormalizeData(
    object = ATAC_obj,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = median(ATAC_obj$nCount_RNA)
)

# 9. select B cells
ATAC_obj <- subset(ATAC_obj, subset = peaks_snn_res.0.8 %in% c(0, 12, 16, 25)) # 0,12,25: B; 16: plasma

# 10. re- normalization and clustering
DefaultAssay(ATAC_obj) <- "peaks"
ATAC_obj <- RunTFIDF(ATAC_obj)
ATAC_obj <- FindTopFeatures(ATAC_obj)
ATAC_obj <- RunSVD(ATAC_obj)
ATAC_obj <- RunUMAP(object = ATAC_obj, reduction = "lsi", dims = 2:30)
ATAC_obj <- RunTSNE(object = ATAC_obj, reduction = "lsi", dims = 2:30)
ATAC_obj <- FindNeighbors(object = ATAC_obj, reduction = "lsi", dims = 2:30)
ATAC_obj <- FindClusters(object = ATAC_obj, verbose = FALSE, algorithm = 3, resolution = 0.8)
ATAC_obj <- NormalizeData(
    object = ATAC_obj,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = median(ATAC_obj$nCount_RNA)
)
## remove doublets
ATAC_obj <- subset(ATAC_obj, subset = peaks_snn_res.0.8 != 7)
## re-normalization and clustering
DefaultAssay(ATAC_obj) <- "peaks"
ATAC_obj <- RunTFIDF(ATAC_obj)
ATAC_obj <- FindTopFeatures(ATAC_obj)
ATAC_obj <- RunSVD(ATAC_obj)
ATAC_obj <- RunUMAP(object = ATAC_obj, reduction = "lsi", dims = 2:30)
ATAC_obj <- RunTSNE(object = ATAC_obj, reduction = "lsi", dims = 2:30)
ATAC_obj <- FindNeighbors(object = ATAC_obj, reduction = "lsi", dims = 2:30)
ATAC_obj <- FindClusters(object = ATAC_obj, verbose = FALSE, algorithm = 3, resolution = 0.8)
ATAC_obj <- FindClusters(object = ATAC_obj, verbose = FALSE, algorithm = 3, resolution = 1.5)
ATAC_obj <- FindClusters(object = ATAC_obj, verbose = FALSE, algorithm = 3, resolution = 2)
ATAC_obj <- NormalizeData(
    object = ATAC_obj,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = median(ATAC_obj$nCount_RNA)
)

# 11. prepare scRNA-seq data for label transfer
B_obj <- read_rds(dir_for_B_obj)
RNA_obj <- subset(B_obj, subset = sampleID %in% unique(ATAC_obj$sampleID))
RNA_obj$annotation_tmp <- as.character(RNA_obj$Annotation_major_2)
RNA_obj$annotation_tmp[as.character(RNA_obj$Annotation) == "c08_ABC_FCRL4"] <- "FCRL4+ Bm"
RNA_obj$annotation_tmp[RNA_obj$annotation_tmp == "Bm"] <- "FCRL4- Bm"

# 12. label transfer
transfer.anchors <- FindTransferAnchors(
    reference = RNA_obj,
    query = ATAC_obj,
    reduction = "cca"
)
predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    refdata = RNA_obj$annotation,
    weight.reduction = ATAC_obj[["lsi"]],
    dims = 2:30
)
ATAC_obj <- AddMetaData(object = ATAC_obj, metadata = predicted.labels)
Idents(ATAC_obj) <- ATAC_obj$peaks_snn_res.0.8
ATAC_obj <- RenameIdents(
    object = ATAC_obj,
    "0" = "FCRL4- Bm",
    "1" = "FCRL4- Bm",
    "2" = "FCRL4- Bm",
    "3" = "ASC",
    "4" = "FCRL4- Bm",
    "5" = "Bn",
    "6" = "Bn",
    "7" = "FCRL4- Bm",
    "8" = "FCRL4- Bm",
    "9" = "Bgc",
    "10" = "FCRL4+ Bm"
)
Idents(ATAC_obj) <- factor(Idents(ATAC_obj), levels = c("Bn", "Bgc", "FCRL4- Bm", "FCRL4+ Bm", "ASC"))
ATAC_obj$annotation <- Idents(ATAC_obj)

# 13. plot tSNE
color <- c(
    "Bn" = "#E18727FF",
    "Bgc" = "#BC3C29FF",
    "FCRL4- Bm" = "#0072B5FF",
    "FCRL4+ Bm" = "#8ca2b4",
    "ASC" = "#20854EFF"
)
DimPlot(
    object = ATAC_obj,
    reduction = "tsne",
    label = FALSE
) +
    ggtitle("scATAC-seq") +
    coord_fixed() +
    scale_color_manual(values = color) +
    theme(
        text = element_text(size = 7, family = "ArialMT"),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        legend.key.width = unit(0.2, "inch"),
        legend.key.height = unit(0.15, "inch"),
        legend.title = element_text(size = 8)
    )
ggsave(file.path(dir_for_result, "6B.ATAC_seq_tSNE.pdf"),
    device = "pdf",
    width = 4,
    height = 3
)

# 14. plot gene peaks
genes <- c("FCRL4", "ITGAX", "CXCR3", "CCR1")
for (i_gene in genes) {
    p <- CoveragePlot(
        object = ATAC_obj,
        assay = "peaks",
        region = i_gene,
        annotation = TRUE,
        peaks = TRUE,
        features = i_gene,
        expression.assay = "RNA",
        extend.upstream = 2000,
        extend.downstream = 2000
    )
    ggsave(file.path(dir_for_result, paste0(i_gene, ".pdf")),
        device = "pdf",
        width = 3,
        height = 3
    )
}

################## ------------------ Figure 6C ------------------ ##################
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
Bm_tumor_tissue_enriched <- subset(Bm_tumor, subset = Annotation %in% c("c05_classical-Bm_GPR183", "c06_Bm_stress-response", "c07_Bm_IFN-response", "c08_ABC_FCRL4", "c11_pre-GC"))

# 4. plot
gene_list <- list(
    "ABC markers" = c("FCRL4", "FCRL5", "ITGAX", "TBX21", "CR2"),
    "Activation" = c("CD80", "CD86", "FAS"),
    "Transcription factors" = c("BHLHE40", "RUNX2", "RUNX3", "SPIB", "ZBTB32", "ZEB2"),
    "Inflammatory cytokine receptors" = c("IFNGR1", "IL2RB", "IL12RB1", "TNFRSF13B", "TNFRSF17"),
    "Chemokine receptors" = c("CXCR3", "CCR1", "CCR5", "CCR6"),
    "MHC-II" = c("HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DQA2", "HLA-DQB2", "HLA-DRA", "HLA-DRB1", "HLA-DOA", "HLA-DOB")
)
ht <- gene_expression_heatmap(Bm_tumor_tissue_enriched, genes = gene_list, group.by = "Annotation", tile_size = 0.12)
pdf(file.path(dir_for_result, "6B.TAAB_signatures.pdf"),
    height = 8,
    width = 4
)
draw(ht)
dev.off()

################## ------------------ Figure 6D ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(Seurat)
library(clusterProfiler)
library(future)

# 2. params
dir_for_B_obj <- ""
dir_for_celltypist_result <- ""
dir_for_result <- ""
source("./functions.R")

# 3. load data
B_obj <- read_rds(dir_for_B_obj)
B_obj_tumor <- subset(B_obj, subset = Tissue_short == "T")

# add celltypist prediction of cycling Bm cells
celltypist_prediction <- read.csv(dir_for_celltypist_result)
celltypist_prediction$celltypist_prediction <- ifelse(celltypist_prediction$predicted_labels == 7, "FCRL4+", "FCRL4-")
B_obj_tumor$celltypist_prediction <- NA
B_obj_tumor$celltypist_prediction[B_obj_tumor$Annotation_major == "Bcycling"] <- celltypist_prediction$celltypist_prediction[match(B_obj_tumor$CellID_final[B_obj_tumor$Annotation_major == "Bcycling"], celltypist_prediction$CellID_final)]
B_obj_tumor$Annotation_major_3 <- as.character(B_obj_tumor$Annotation_major_2)
B_obj_tumor$Annotation_major_3[B_obj_tumor$Annotation_major_2 == "Bm"] <- "FCRL4- Bm"
B_obj_tumor$Annotation_major_3[B_obj_tumor$Annotation == "c08_ABC_FCRL4"] <- "FCRL4+ Bm"
B_obj_tumor$Annotation_major_3[B_obj_tumor$Annotation == "c14_Bm_activated−cycling" & B_obj_tumor$celltypist_prediction == "FCRL4+"] <- "FCRL4+ Bm"

# 4. differential expression (FCRL4+ vs FCRL4- Bm)
plan("multisession", workers = 5)
plan()
options(future.globals.maxSize = 10 * 1024^3) # 10 GB
markers <- FindMarkers(B_obj_tumor, ident.1 = "FCRL4+ Bm", ident.2 = "FCRL4- Bm", group.by = "Annotation_major_3", only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox")
markers <- markers %>%
    rownames_to_column("gene")
plan("sequential")

# 5. GO enrichment
GO_result_up <- GO_enrichment(markers$gene, rownames(B_obj_tumor))

# 6. plot
plot_df <- GO_result_up$go_res@result %>%
    mutate(log10q = -log10(qvalue)) %>%
    arrange(desc(log10q)) %>%
    slice_head(n = 17)
plot_df$Description <- with(plot_df, reorder(Description, log10q))
ggplot(plot_df, aes(x = Description, y = log10q)) +
    geom_bar(stat = "identity", fill = "#FFADAD") +
    theme_minimal() +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(family = "ArialMT", size = 6, color = "black"),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_line(colour = "black"),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        text = element_text(family = "ArialMT", size = 6),
        plot.title = element_text(hjust = 0.5, size = 8),
        panel.grid = element_blank(),
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
    ) +
    ylim(-35, 16) +
    geom_text(data = plot_df, aes(x = Description, y = -0.1, label = Description), hjust = 1, color = "black", size = 2) +
    coord_flip() +
    ylab("-Log10(q-value)") +
    ggtitle("TAAB-enriched pathways")
ggsave(file.path(dir_for_result, "6D.TAAB_pathways.pdf"),
    device = "pdf",
    width = 4.5,
    height = 2.5
)

################## ------------------ Figure 6E ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(Seurat)
library(AUCell)
library(SCENIC)

# 2. params
dir_for_B_obj <- ""
dir_for_celltypist_result <- ""
dir_for_SCENIC <- ".../scenic_auc_mtx.txt"
dir_for_result <- ""

# 3. load data
B_obj <- read_rds(dir_for_B_obj)
B_obj_tumor <- subset(B_obj, subset = Tissue_short == "T")

# add celltypist prediction of cycling Bm cells
celltypist_prediction <- read.csv(dir_for_celltypist_result)
celltypist_prediction$celltypist_prediction <- ifelse(celltypist_prediction$predicted_labels == 7, "FCRL4+", "FCRL4-")
B_obj_tumor$celltypist_prediction <- NA
B_obj_tumor$celltypist_prediction[B_obj_tumor$Annotation_major == "Bcycling"] <- celltypist_prediction$celltypist_prediction[match(B_obj_tumor$CellID_final[B_obj_tumor$Annotation_major == "Bcycling"], celltypist_prediction$CellID_final)]
B_obj_tumor$Annotation_major_3 <- as.character(B_obj_tumor$Annotation_major_2)
B_obj_tumor$Annotation_major_3[B_obj_tumor$Annotation_major_2 == "Bm"] <- "FCRL4- Bm"
B_obj_tumor$Annotation_major_3[B_obj_tumor$Annotation == "c08_ABC_FCRL4"] <- "FCRL4+ Bm"
B_obj_tumor$Annotation_major_3[B_obj_tumor$Annotation == "c14_Bm_activated−cycling" & B_obj_tumor$celltypist_prediction == "FCRL4+"] <- "FCRL4+ Bm"

# 4. run SCENIC
# yangyu

# 5. load SCENIC result
regulonAUC <- data.table::fread(dir_for_SCENIC)
regulonAUC <- regulonAUC[regulonAUC$V1 %in% B_obj_tumor$CellID, ]
regulonAUC <- regulonAUC %>%
    column_to_rownames("V1") %>%
    as.matrix() %>%
    t()
meta <- data.frame(
    cellID = colnames(regulonAUC),
    Annotation = as.character(B_obj_tumor$Annotation[match(colnames(regulonAUC), B_obj_tumor$CellID)]),
    Annotation_tmp = as.character(B_obj_tumor$Annotation_major_3[match(colnames(regulonAUC), B_obj_tumor$CellID)])
)
# subset Bm
subset <- meta$Annotation_tmp %in% c("FCRL4- Bm,", "FCRL4+ Bm")
regulonAUC_subset <- regulonAUC[, subset]
meta_subset <- meta[subset, ]
meta_subset$Annotation[meta_subset$Annotation_tmp == "FCRL4+ Bm"] <- "FCRL4+ Bm"

# 6. calculate RSS
rss <- calcRSS(AUC = regulonAUC_subset, cellAnnotation = meta_subset$Annotation)

# 7. plot
setName <- "FCRL4+ Bm"
n <- 5
# sort
rssThisType <- sort(rss[, setName], decreasing = TRUE)
# label
thisRss <- data.frame(
    regulon = names(rssThisType), rank = seq_along(rssThisType),
    rss = rssThisType
)
thisRss$regulon[(n + 1):nrow(thisRss)] <- NA
# color
thisRss$color <- FALSE
thisRss$color[(n + 1):nrow(thisRss)] <- TRUE
# plot
ggplot(thisRss, aes(x = rank, y = rss, color = color)) +
    geom_point(size = 0.5, stroke = 0, shape = 16) +
    scale_color_manual(values = c("#C3423F", "#807F7F")) +
    ggrepel::geom_text_repel(aes(label = regulon), segment.color = "grey40", size = 0.35 * 5) +
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
    xlab("Rank") +
    ylab("Regulon specificity score") +
    ggtitle("TAAB-specific regulons")
ggsave(file.path(dir_for_result, "6E.TAAB_SCENIC.pdf"),
    width = 2.5,
    height = 2.5
)
