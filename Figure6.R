###################################### Figure 6 ######################################

################## ------------------ Figure 6A ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(Seurat)

# 2. params
dir_for_B_obj <- "data/All_obj.rds"
dir_for_celltypist_result <- "data/predict_pro_memory.csv"
dir_for_result <- "figures"
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
B_obj_tumor <- subset(B_obj, subset = Tissue_short == "T")

# add celltypist prediction of cycling Bm cells
celltypist_prediction <- read.csv(dir_for_celltypist_result)
celltypist_prediction$celltypist_prediction <- ifelse(celltypist_prediction$predicted_labels == 7, "FCRL4+", "FCRL4-")
B_obj_tumor$celltypist_prediction <- NA
B_obj_tumor$celltypist_prediction[B_obj_tumor$Annotation_major == "Bcycling"] <- celltypist_prediction$celltypist_prediction[match(B_obj_tumor$CellID_final[B_obj_tumor$Annotation_major == "Bcycling"], celltypist_prediction$CellID_final)]
B_obj_tumor$Annotation_major_3 <- as.character(B_obj_tumor$Annotation_major_2)
B_obj_tumor$Annotation_major_3[B_obj_tumor$Annotation_major_2 == "Bm"] <- "FCRL4- Bm"
B_obj_tumor$Annotation_major_3[B_obj_tumor$Annotation == "c08_ABC_FCRL4"] <- "FCRL4+ Bm"
B_obj_tumor$Annotation_major_3[B_obj_tumor$Annotation == "c14_Bm_activated-cycling" & B_obj_tumor$celltypist_prediction == "FCRL4+"] <- "FCRL4+ Bm"

# 4. plot
plot_df <- B_obj_tumor@meta.data %>%
    dplyr::filter(Treatment_status == "treatment naïve") %>%
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
library(future)

# 2. params
dir_for_peak_counts <- "data/ATAC_seq/filtered_peak_bc_matrix.h5"
dir_for_fragments <- "data/ATAC_seq/fragments.tsv.gz"
dir_for_meta <- "data/ATAC_seq/singlecell.csv"
dir_for_B_obj <- "data/All_obj.rds"
dir_for_data <- "data"
dir_for_result <- "figures"

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
ATAC_obj <- subset(ATAC_obj, subset = peaks_snn_res.0.8 %in% c(0, 12, 16, 25)) # 0,12,25: CD20+ B; 16: ASC

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
# write_rds(ATAC_obj, file.path(dir_for_data, "BRCA_ATAC_seq_obj_before_label_transfer.rds"), compress = "gz")

# 11. prepare scRNA-seq data for label transfer
B_obj <- read_rds(dir_for_B_obj)
RNA_obj <- subset(B_obj, subset = sampleID %in% unique(ATAC_obj$sampleID))

# 12. label transfer
DefaultAssay(ATAC_obj) <- "RNA"
transfer.anchors <- FindTransferAnchors(
    reference = RNA_obj,
    query = ATAC_obj,
    reduction = "cca"
)
predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    refdata = RNA_obj$Annotation,
    weight.reduction = ATAC_obj[["lsi"]],
    dims = 2:30
)
ATAC_obj <- AddMetaData(object = ATAC_obj, metadata = predicted.labels)
Idents(ATAC_obj) <- ATAC_obj$peaks_snn_res.0.8

## confusion matrix
res <- 0.8
Var1 <- ATAC_obj@meta.data[[paste0("peaks_snn_res.", res)]]
Var2 <- ATAC_obj@meta.data[["predicted.id"]]
plot_df <- as.data.frame(table(Var1, Var2), stringsAsFactors = FALSE)
tmp <- table(Var1, Var2)
tmp <- tmp / rowSums(tmp)
plot_df_2 <- as.data.frame(tmp, stringsAsFactors = FALSE)
plot_df$percentage <- plot_df_2$Freq
plot_df$Var1 <- factor(plot_df$Var1, levels = levels(Var1))
ggplot(plot_df, aes(Var1, Var2, fill = percentage)) +
    geom_tile() +
    geom_text(aes(label = ifelse(percentage > 0.1, scales::percent(percentage, accuracy = 0.1), "")), size = 2, family = "Arial") +
    scale_fill_gradient(low = "white", high = "#3575b5") +
    theme_minimal() +
    coord_equal() +
    theme(
        text = element_text(size = 10, family = "Arial"),
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, size = 10),
        plot.margin = unit(c(1, 1, 1, 1), "char")
    ) +
    labs(x = paste0("peaks_snn_res.", res), y = "Predicted cell type")

## set annotation
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
# write_rds(ATAC_obj, file.path(dir_for_data, "BRCA_ATAC_seq_obj.rds"), compress = "gz")

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
    ggsave(file.path(dir_for_result, paste0("6B.ATAC_seq_", i_gene, ".pdf")),
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
dir_for_B_obj <- "data/All_obj.rds"
dir_for_result <- "figures"
source("code/functions.R")

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
pdf(file.path(dir_for_result, "6C.TAAB_signatures.pdf"),
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
dir_for_B_obj <- "data/All_obj.rds"
dir_for_celltypist_result <- "data/predict_pro_memory.csv"
dir_for_result <- "figures"
source("code/functions.R")

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
B_obj_tumor$Annotation_major_3[B_obj_tumor$Annotation == "c14_Bm_activated-cycling" & B_obj_tumor$celltypist_prediction == "FCRL4+"] <- "FCRL4+ Bm"

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
    slice_head(n = 18)
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
    height = 2.6
)

################## ------------------ Figure 6E ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(Seurat)
library(AUCell)
library(SCENIC)

# 2. params
dir_for_B_obj <- "data/All_obj.rds"
dir_for_celltypist_result <- "data/predict_pro_memory.csv"
dir_for_SCENIC <- "data/SCENIC/scenic_auc_mtx.txt"
dir_for_result <- "figures"

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
B_obj_tumor$Annotation_major_3[B_obj_tumor$Annotation == "c14_Bm_activated-cycling" & B_obj_tumor$celltypist_prediction == "FCRL4+"] <- "FCRL4+ Bm"

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
subset <- meta$Annotation_tmp %in% c("FCRL4- Bm", "FCRL4+ Bm")
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

################## ------------------ Figure 6G ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library("Startrac")

# 2. params
dir_for_BCR <- "data/panB_BCR_20240709.rds"
dir_for_result <- "figures"

# 3. load data
clone <- read_rds(dir_for_BCR)
clone_subset <- clone[clone$Tissue_short == "T", ]

# 4. run STARTRAC
startrac_input <- data.frame(
    Cell_Name = clone_subset$sequence_id,
    clone.id = clone_subset$clone_id,
    patient = clone_subset$PatientID,
    majorCluster = clone_subset$Annotation,
    loc = clone_subset$Tissue_short
)
startrac_output <- Startrac.run(startrac_input, proj = "panB", verbose = T)

# 5. plot
target <- "c08_ABC_FCRL4"
plot_df <- startrac_output@pIndex.sig.tran[startrac_output@pIndex.sig.tran$aid == "panB", ]
plot_df <- plot_df %>%
    dplyr::filter(majorCluster == target, index != target)
plot_df$index <- with(plot_df, reorder(index, value, function(x) -x))
plot_df$color <- plot_df$value > 0.06
# only consider CD20+ B cells, remove ASCs for simplicity of visualization
plot_df <- plot_df %>%
    dplyr::filter(!index %in% c(
        "c15_cycling_ASC",
        "c16_PC_IGHG",
        "c17_PC_IGHA",
        "c18_early-PC_MS4A1low",
        "c19_early-PC_LTB",
        "c20_early-PC_RGS13"
    ))

ggplot(plot_df, aes(x = index, y = value, fill = color)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("#5BC0EB", "#C3423F")) +
    cowplot::theme_cowplot() +
    theme(
        axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 7),
        text = element_text(size = 7, family = "ArialMT"),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        plot.title = element_text(hjust = 0.5, size = 8, family = "ArialMT", face = "plain"),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
    ) +
    ylab("pTrans of TAABs") +
    xlab("")

ggsave(file.path(dir_for_result, "6G.TAAB_pTrans.pdf"),
    width = 2.2,
    height = 2.5
)

################## ------------------ Figures 6I, S6D and S6E ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(Seurat)

# 2. params
dir_for_B_obj <- "data/All_obj.rds"
dir_for_diff <- "data/trajectory_analysis/diffMap.txt"
dir_for_pseudotime <- "data/trajectory_analysis/dpt_pseudotime.txt"
dir_for_result <- "figures"

# palette
celltype_color_panel <- c(
    "c04_classical-Bm_TXNIP" = "#a82d06",
    "c07_Bm_IFN-response" = "#74a764",
    "c08_ABC_FCRL4" = "#8ca2b4"
)

# 3. load data
B_obj <- read_rds(dir_for_B_obj)
Bm_obj <- subset(B_obj, subset = Annotation_major == "Bm")
Bm_obj <- Seurat::ScaleData(Bm_obj, features = rownames(Bm_obj))
Bm_subset <- subset(Bm_obj, subset = Reference == "thisstudy" & Tissue == "Tumor" & Annotation %in% c("c04_classical-Bm_TXNIP", "c07_Bm_IFN-response", "c08_ABC_FCRL4"))

# 4. add diffusion components and pseudotime
# add diffusion components
diffmap <- data.table::fread(dir_for_diff, header = T, stringsAsFactors = F, sep = "\t")
diffmap <- diffmap[match(Bm_subset$CellID_final, diffmap$CellID_final), ]
diffmap <- diffmap[, c("1", "2")] %>% as.data.frame()
colnames(diffmap) <- c("DC_1", "DC_2")
rownames(diffmap) <- colnames(Bm_subset)
Bm_subset@reductions$diffmap <- CreateDimReducObject(embeddings = as.matrix(diffmap), key = "DC_", assay = "RNA")
# add diffusion pseudotime
dpt_pseudotime <- data.table::fread(dir_for_pseudotime, header = T, stringsAsFactors = F, sep = "\t")
Bm_subset$dpt_pseudotime <- dpt_pseudotime$dpt_pseudotime[match(Bm_subset$CellID_final, dpt_pseudotime$CellID_final)]
## convert dpt to percentile
convert_pseudotime_to_percentile <- function(pseudotime) {
    return(rank(pseudotime) / length(pseudotime))
}
Bm_subset$dpt_order_perc <- convert_pseudotime_to_percentile(Bm_subset$dpt_pseudotime)

# 4. plot pseudotime density
plot_df <- data.frame(
    Annotation = Bm_subset$Annotation,
    dpt_pseudotime = Bm_subset$dpt_pseudotime,
    dpt_order_perc = Bm_subset$dpt_order_perc
)

plot_df <- plot_df %>%
    dplyr::filter(!is.na(dpt_pseudotime))

ggplot(plot_df, aes(x = dpt_order_perc, fill = Annotation, color = Annotation)) +
    geom_density(adjust = 1.5, trim = F, alpha = 0.6) +
    scale_fill_manual(values = celltype_color_panel) +
    scale_color_manual(values = celltype_color_panel) +
    xlim(0, 1) +
    cowplot::theme_cowplot() +
    theme(
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        text = element_text(size = 7, family = "ArialMT"),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        plot.title = element_text(hjust = 0.5, size = 8, family = "ArialMT", face = "plain"),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
    )

ggsave(file.path(dir_for_result, "S6D.pseudotime_density.pdf"),
    width = 3,
    height = 2
)

# 5. plot gene expression dynamics along the pseudotime
## 6I
gene_list <- list(
    "Activation & costimulation" = c("CD86", "FAS", "LGALS9", "ICAM1"),
    "Antigen presentation" = c("HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1")
)

for (i_name in names(gene_list)) {
    genes <- gene_list[[i_name]]
    plot_df <- Bm_subset@assays$RNA@scale.data[genes, ]
    colnames(plot_df) <- Bm_subset$CellID_final
    plot_df <- reshape2::melt(plot_df)
    plot_df$pseudotime <- Bm_subset$dpt_order_perc[match(plot_df$Var2, Bm_subset$CellID_final)]
    plot_df <- plot_df[!is.na(plot_df$value), ]

    ggplot(data = plot_df, mapping = aes(x = pseudotime, y = value, color = Var1, fill = Var1)) +
        # geom_point() +
        ggsci::scale_color_npg(name = "") +
        ggsci::scale_fill_npg(name = "") +
        geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), alpha = 0.2, size = 0.5) +
        cowplot::theme_cowplot() +
        theme(
            axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 7),
            text = element_text(size = 7, family = "ArialMT"),
            plot.margin = unit(c(1, 1, 1, 1), "char"),
            plot.title = element_text(hjust = 0.5, size = 8, family = "ArialMT", face = "plain"),
            axis.line = element_line(linetype = 1, color = "black", size = 0.3),
            axis.ticks = element_line(linetype = 1, color = "black", size = 0.3),
            legend.position = "none"
        ) +
        ylab("Scaled expression") +
        xlab("Pseudo-time") +
        ggtitle(i_name)
    ggsave(file.path(dir_for_result, paste0("6I.", i_name, "_gene_trends.pdf")),
        width = 1.8,
        height = 2
    )
}

## S6E
gene_list <- list(
    "ABC markers" = c("FCRL4", "FCRL5", "ITGAX", "TBX21", "CR2"),
    "Transcription factors" = c("TBX21", "BHLHE40", "SPIB", "RUNX2", "RUNX3", "CREM", "ZBTB32"),
    "Inflammatory cytokine receptors" = c("CCR1", "CXCR3", "CXCR5", "IFNGR1", "IL21R")
)

for (i_name in names(gene_list)) {
    genes <- gene_list[[i_name]]
    plot_df <- Bm_subset@assays$RNA@scale.data[genes, ]
    colnames(plot_df) <- Bm_subset$CellID_final
    plot_df <- reshape2::melt(plot_df)
    plot_df$pseudotime <- Bm_subset$dpt_order_perc[match(plot_df$Var2, Bm_subset$CellID_final)]
    plot_df <- plot_df[!is.na(plot_df$value), ]

    ggplot(data = plot_df, mapping = aes(x = pseudotime, y = value, color = Var1, fill = Var1)) +
        # geom_point() +
        ggsci::scale_color_npg(name = "") +
        ggsci::scale_fill_npg(name = "") +
        geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), alpha = 0.2, size = 0.5) +
        cowplot::theme_cowplot() +
        theme(
            axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 7),
            text = element_text(size = 7, family = "ArialMT"),
            plot.margin = unit(c(1, 1, 1, 1), "char"),
            plot.title = element_text(hjust = 0.5, size = 8, family = "ArialMT", face = "plain"),
            axis.line = element_line(linetype = 1, color = "black", size = 0.3),
            axis.ticks = element_line(linetype = 1, color = "black", size = 0.3),
            legend.position = "none"
        ) +
        ylab("Scaled expression") +
        xlab("Pseudo-time") +
        ggtitle(i_name)
    ggsave(file.path(dir_for_result, paste0("S6E.", i_name, "_gene_trends.pdf")),
        width = 1.8,
        height = 2
    )
}

################## ------------------ Figures 6J, S6F-I ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(Seurat)

# 2. params
dir_for_non_cancer_B_obj <- "data/non_cancer_B_obj_20240628.rds"
dir_for_cancer_B_obj <- "data/All_obj.rds"
dir_for_result <- "figures"
source("code/functions.R")

## palette
celltype_color_panel <- c(
    "Bn" = "#E18727FF",
    "ISG+ Bn" = "#e78071",
    "Bm" = "#0072B5FF",
    "atypical-like B" = "#8ca2b4",
    "Bgc" = "#BC3C29FF",
    "Bcycling" = "#7876B1FF",
    "ASC" = "#20854EFF"
)

# 3. load data
non_cancer_B <- read_rds(dir_for_non_cancer_B_obj)
cancer_B <- read_rds(dir_for_cancer_B_obj)

# 4. plot annotations
non_cancer_B$annotation <- factor(non_cancer_B$annotation, levels = c("Bn", "ISG+ Bn", "Bm", "atypical-like B", "Bgc", "Bcycling", "ASC"))
custom_dimplot(non_cancer_B, group.by = "annotation", label = FALSE) +
    scale_color_manual(values = celltype_color_panel)
ggsave(file.path(dir_for_result, "S6F.non_cancer_B_annotation.png"),
    height = 4, width = 8
)

# 5. plot datasets
## add sample number for each dataset
sample_df <- non_cancer_B@meta.data %>%
    group_by(datasource) %>%
    summarise(n_sample = n_distinct(SampleID))
non_cancer_B$Datasets <- paste(non_cancer_B$datasource, " N=", sample_df$n_sample[match(non_cancer_B$datasource, sample_df$datasource)])
## plot
custom_dimplot(non_cancer_B, group.by = "Datasets", label = FALSE)
ggsave(file.path(dir_for_result, "S6H.non_cancer_B_datasets.png"),
    height = 4, width = 8
)

# 6. plot atypical B cell markers
custom_dotplot(non_cancer_B, genes = c("FCRL4", "FCRL5", "ITGAX", "TBX21"), group.by = "annotation") +
    scale_radius(range = c(0.2, 4))
ggsave(file.path(dir_for_result, "S6G.non_cancer_B_ABC_markers.pdf"),
    width = 3,
    height = 2
)

# 7. atypical B cell marker fold change (atypical Bm cells vs other Bm cells)
genes <- c("FCRL4", "FCRL5", "ITGAX", "TBX21")

## calculate non-cancer
fc_df <- list()
for (i_tissue in unique(non_cancer_B$tissue)) {
    data_tmp <- subset(non_cancer_B, subset = tissue == i_tissue)
    fc_tmp <- seurat_foldchange(data_tmp, features = genes, ident.1 = "atypical-like B", ident.2 = "Bm", group.by = "annotation")
    fc_tmp$tissue <- i_tissue
    fc_tmp$tissue_group <- unique(data_tmp$tissue_group)
    fc_df <- append(fc_df, list(fc_tmp))
}
## calculate cancer
data_tmp <- subset(cancer_B, subset = Tissue == "Tumor" & Treatment_status == "treatment naïve" & Annotation_major == "Bm")
data_tmp$annotation <- ifelse(data_tmp$Annotation == "c08_ABC_FCRL4", "c08_ABC_FCRL4", "Bm")
fc_tmp <- seurat_foldchange(data_tmp, features = genes, ident.1 = "c08_ABC_FCRL4", ident.2 = "Bm", group.by = "annotation")
fc_tmp$tissue <- "Tumor tissue"
fc_tmp$tissue_group <- "Tumor tissue"
fc_df <- append(fc_df, list(fc_tmp))

fc_df <- bind_rows(fc_df)

## plot
plot_df <- fc_df[fc_df$tissue != "PA synovial tissue", ] # only have one atypical B cell
plot_df$tissue <- factor(plot_df$tissue,
    levels = c(
        "Tumor tissue",
        "SLE kidney", "RA synovial tissue", "Healthy tonsil", "SSc-ILD lung", "Influenza vaccination lymph node", "Healthy lymphoid tissue", "Healthy non-lymphoid tissue",
        "HIV blood", "Malaria blood", "SLE blood", "SARS-CoV-2 blood", "Healthy blood", "Influenza vaccination blood"
    )
)

ggplot(plot_df, aes(x = tissue, y = avg_log2FC)) +
    geom_bar(stat = "identity", fill = "#5BC0EB") +
    facet_wrap(~gene, nrow = 4) +
    cowplot::theme_cowplot() +
    theme(
        axis.text.x = element_text(size = 7, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 7),
        text = element_text(size = 7, family = "ArialMT"),
        plot.title = element_text(hjust = 0.5, size = 8, face = "plain"),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.3)
    ) +
    xlab("Tissue origins of atypical-like B cells") +
    ylab("Average log2(fold change) vs other Bm cells")

ggsave(file.path(dir_for_result, "6J_S6I.non_cancer_B_ABC_markers_fold_change.pdf"),
    width = 3.5,
    height = 5
)

# 8. FCRL4 fold change (atypical Bm cells vs other Bm cells) per cancer
## calculate
cancer_B_subset <- subset(cancer_B, subset = Tissue == "Tumor" & Treatment_status == "treatment naïve" & Annotation_major == "Bm")
fc_df <- list()
for (i_cancer in unique(cancer_B_subset$Cancer)) {
    data_tmp <- subset(cancer_B_subset, subset = Cancer == i_cancer)
    data_tmp$annotation <- ifelse(data_tmp$Annotation == "c08_ABC_FCRL4", "c08_ABC_FCRL4", "Bm")
    fc_tmp <- seurat_foldchange(data_tmp, features = genes, ident.1 = "c08_ABC_FCRL4", ident.2 = "Bm", group.by = "annotation")
    fc_tmp$Cancer <- i_cancer
    fc_df <- append(fc_df, list(fc_tmp))
}
fc_df <- bind_rows(fc_df)

## only cancers with at least five TAABs
select_cancers <- as.data.frame(table(cancer_B_subset$Annotation, cancer_B_subset$Cancer)) %>%
    dplyr::filter(
        Var1 == "c08_ABC_FCRL4",
        Freq >= 5
    )
select_cancers <- as.character(select_cancers$Var2)
plot_df <- fc_df %>%
    dplyr::filter(Cancer %in% select_cancers)

## plot
order_df <- plot_df %>%
    dplyr::filter(gene == "FCRL4") %>%
    arrange(desc(avg_log2FC))
plot_df$Cancer <- factor(plot_df$Cancer, levels = order_df$Cancer)

ggplot(plot_df, aes(x = Cancer, y = avg_log2FC)) +
    geom_bar(stat = "identity", fill = "#5BC0EB") +
    facet_wrap(~gene, nrow = 4) +
    cowplot::theme_cowplot() +
    theme(
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        text = element_text(size = 7, family = "ArialMT"),
        plot.title = element_text(hjust = 0.5, size = 8, face = "plain"),
        legend.position = "none",
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.3)
    ) +
    xlab("Cancer") +
    ylab("Average log2(fold change) vs other Bm cells")

ggsave(file.path(dir_for_result, "6J.cancer_B_ABC_markers_fold_change.pdf"),
    width = 3.58,
    height = 4.5
)

################## ------------------ Figure 6L ------------------ ##################
# 1. library
library(readr)
library(tidyverse)

# 2. params
dir_for_induction_results <- "data/FCRL4_induction_experiment.csv"
dir_for_result <- "figures"

# 3. load data
data <- read_csv(dir_for_induction_results)
data <- data %>%
    tidyr::gather(., key = "condition", value = "FCRL4_pct") %>%
    dplyr::filter(!is.na(FCRL4_pct))

# 4. calculate significance
conditions <- unique(data$condition)
significance_df <- as.data.frame(t(combn(conditions, 2)))
significance_df$p <- lapply(1:nrow(significance_df), function(x) {
    cor_tmp <- t.test(
        data$FCRL4_pct[data$condition == significance_df$V1[x]],
        data$FCRL4_pct[data$condition == significance_df$V2[x]]
    )
    return(cor_tmp$p.value)
})
significance_df$significance <- "ns"
significance_df$significance[significance_df$p < 0.05] <- "*"
significance_df$significance[significance_df$p < 0.01] <- "**"
significance_df$significance[significance_df$p < 0.001] <- "***"
significance_df$significance[significance_df$p < 0.0001] <- "****"

# 5. plot
data_summary <- function(data, varname, groupnames) {
    require(plyr)
    summary_func <- function(x, col) {
        c(
            mean = mean(x[[col]], na.rm = TRUE),
            sem = sd(x[[col]], na.rm = TRUE) / sqrt(length(x[[col]]))
        )
    }
    data_sum <- ddply(data, groupnames,
        .fun = summary_func,
        varname
    )
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
}
plot_df <- data_summary(data, varname = "FCRL4_pct", groupnames = "condition")
plot_df$condition <- factor(plot_df$condition, levels = c("Medium", "SW480-SN", "Huh7-SN", "Poly I:C", "LPS", "CEA", "anti-CD40", "CpG", "CpG+anti-CD40", "CpG+anti-CD40+CEA"))

ggplot(plot_df, aes(x = condition, y = FCRL4_pct, fill = condition)) +
    geom_errorbar(aes(ymin = FCRL4_pct - sem, ymax = FCRL4_pct + sem),
        width = .2,
        position = position_dodge(.9),
        size = 0.3
    ) +
    geom_bar(
        stat = "identity",
        position = position_dodge()
    ) +
    paletteer::scale_fill_paletteer_d("impressionist.colors::chanteuse_de_cafe_concert") +
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
    xlab("") +
    ylab("FCRL4+ in B cells (%)")
ggsave(file.path(dir_for_result, "6L.FCRL4_induction_experiment.pdf"),
    width = 2.8,
    height = 2.5
)
