fdir = ""
## I. load packages
library(readr)
library(dplyr)
library(ggplot2)
library(scales)

## II. load data
obj = read_rds("obj/All_obj.rds")
meta = obj@meta.data
BCR = read_rds("BCR.rds")
IHC_info = read_tsv("IHC_panel.tsv")

## III. load color
Tissue_color_panel = c(
  "Tumor" = "#5BC0EB",
  "Adjacent non-tumor" = "#9BC53D",
  "Blood" = "#C3423F",
  "Other tissue" = "#FDE74C"
)

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

Major_color_panel = c(
  "Bn" =  "#E18727FF",
  "Bm" = "#0072B5FF",
  "Bgc" = "#BC3C29FF",
  "ASC" = "#20854EFF",
  "Bcycling" = "#7876B1FF"
)

Subset_color_panel <- c(
  # Bn
  "c01_Bn_TCL1A" = "#e5d25b",
  "c02_Bn_NR4A2" = "#599014",
  "c03_Bn_IFN-response" = "#e78071",
  # Bm
  "c04_classical-Bm_TXNIP" = "#a82d06",
  "c05_classical-Bm_GPR183" = "#4592bf",
  "c06_Bm_stress-response" = "#d38219",
  "c07_Bm_IFN-response" = "#74a764",
  "c08_ABC_FCRL4" = "#8ca2b4",
  "c09_ABC_FGR" = "#cbb190",
  "c10_Bm_TCL1A" = "#e7ca8d",
  "c11_pre-GC" = "#9d9ec3",
  # Bgc
  "c12_Bgc_LZ-like" = "#593202",
  # ASC
  "c16_PC_IGHG" = "#ebafa4",
  "c17_PC_IGHA" = "#5e8a89",
  "c18_early-PC_MS4A1low" = "#ecd577",
  "c19_early-PC_LTB" = "#7c606c",
  "c20_early-PC_RGS13" = "#5c6489",
  # Bcycling
  "c13_Bgc_DZ-like" = "#ECE4B7",
  "c15_cycling_ASC" = "#D36135",
  "c14_Bm_activated-cycling" = "#467599"
)


## IV. load code
yy_Dotplot <- function(seuratObj,
                       genes, # list (grouped) or vector (ungrouped)
                       group.by,
                       coord_flip = FALSE,
                       scale = TRUE,
                       dot.scale = 4,
                       gene_expr_cutoff = -Inf,
                       gene_pct_cutoff = 0,
                       cell_expr_cutoff = -Inf,
                       cell_pct_cutoff = 0,
                       panel.spacing_distance = 0.5,
                       return_data = FALSE) {
  # Required libraries
  require(dplyr)
  require(Seurat)
  require(RColorBrewer)
  require(ggplot2)
  require(cowplot)
  require(tidyverse)
  
  # Parameters
  col.min <- -2.5
  col.max <- 2.5
  dot.min <- 0
  
  # Reformat genes
  if (is.list(genes)) {
    if (is.null(names(genes))) {
      names(genes) <- 1:length(genes)
    }
    genes <- stack(genes)
    features <- genes$values
    features_group <- genes$ind
  } else {
    features <- genes
    features_group <- NULL
  }
  
  # Check for missing genes
  subset <- features %in% rownames(seuratObj)
  if (sum(!subset) > 0) {
    cat(paste0(paste(features[!subset], collapse = ", "), " is missing in gene list"))
  }
  if (length(features) == 0) {
    stop("No intersecting genes, please check gene name format.\n")
  }
  
  # Prepare plot input
  data.features <- FetchData(seuratObj, cells = colnames(seuratObj), vars = features, slot = "data")
  data.features$id <- seuratObj@meta.data[[group.by]]
  
  data.plot <- lapply(unique(data.features$id), function(ident) {
    data.use <- data.features[data.features$id == ident, 1:(ncol(data.features) - 1), drop = FALSE]
    avg.exp <- apply(data.use, 2, function(x) mean(expm1(x)))
    pct.exp <- apply(data.use, 2, function(x) sum(x > 0) / length(x))
    list(avg.exp = avg.exp, pct.exp = pct.exp)
  })
  names(data.plot) <- unique(data.features$id)
  data.plot <- lapply(names(data.plot), function(x) {
    data.use <- as.data.frame(data.plot[[x]])
    data.use$features.plot <- rownames(data.use)
    data.use$features.plot_show <- features
    data.use$id <- x
    data.use
  })
  data.plot <- do.call(rbind, data.plot)
  
  avg.exp.scaled <- sapply(unique(data.plot$features.plot), function(x) {
    data.use <- data.plot[data.plot$features.plot == x, "avg.exp"]
    if (scale) {
      data.use <- scale(data.use)
      MinMax(data.use, min = col.min, max = col.max)
    } else {
      log1p(data.use)
    }
  })
  avg.exp.scaled <- as.vector(t(avg.exp.scaled))
  
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(data.plot$features.plot, levels = unique(data.plot$features.plot))
  data.plot$features.plot_show <- factor(data.plot$features.plot_show, levels = unique(features))
  
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  
  if (is.factor(seuratObj@meta.data[[group.by]])) {
    data.plot$id <- factor(data.plot$id, levels = levels(seuratObj@meta.data[[group.by]]))
  }
  
  if (!is.null(features_group)) {
    data.plot <- data.plot %>%
      left_join(data.frame(features = features, features_group = features_group), by = c("features.plot_show" = "features")) %>%
      mutate(features_group = factor(features_group, levels = unique(features_group)))
  }
  
  # Filter genes/cells
  filter_cell <- data.plot %>%
    group_by(id) %>%
    summarise(max_exp = max(avg.exp.scaled), max_pct = max(pct.exp)) %>%
    filter(max_exp >= cell_expr_cutoff & max_pct >= cell_pct_cutoff)
  filter_gene <- data.plot %>%
    group_by(features.plot) %>%
    summarise(max_exp = max(avg.exp.scaled), max_pct = max(pct.exp)) %>%
    filter(max_exp >= gene_expr_cutoff & max_pct >= gene_pct_cutoff)
  
  data.plot <- data.plot %>%
    filter(features.plot %in% filter_gene$features.plot) %>%
    filter(id %in% filter_cell$id)
  
  # Plot
  plot <- ggplot(data.plot, aes(x = features.plot, y = id)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_radius(range = c(0, dot.scale)) +
    scale_color_gradientn(colors = brewer.pal(9, "RdBu")[9:1]) +
    guides(size = guide_legend(title = "Percent expressed"), color = guide_colorbar(title = "Average expression")) +
    theme(
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      text = element_text(size = 8),
      plot.margin = unit(c(1, 1, 1, 1), "char"),
      panel.background = element_rect(colour = "black", fill = "white"),
      panel.grid = element_line(colour = "grey", linetype = "dashed", size = 0.2)
    ) +
    labs(x = "", y = "")
  
  if (coord_flip) {
    plot <- plot + coord_flip()
  }
  
  if (!is.null(features_group)) {
    plot <- plot +
      facet_grid(. ~ features_group, scales = "free_x", space = "free_x", switch = "y") +
      theme(panel.spacing = unit(panel.spacing_distance, "lines"), strip.background = element_blank())
  }
  
  if (return_data) {
    return(list(plot = plot, data = data.plot))
  } else {
    return(plot)
  }
}

#----- 01. Constitution of the atlas -----#

## Figure 1A; Schematics overview of our atlas
### add label: This study; In-house; Published
meta$new = ifelse(
  meta$Reference == "thisstudy",
  yes = "This study",
  no = ifelse(
    test = meta$Reference %in% c(
      "Kang,B.-2022-GenomeBio",
      "Zhang,L.-2020-Cell",
      "Zhang,Y.-2021-CancerCell",
      "Zhang,Q.-2019-Cell",
      "Liu,Y.-2021-NatCommun"
    ),
    yes = "In-house",
    no = "Published"
  )
)
meta$new = factor(meta$new, levels = c("This study","In-house","Published"))

### 01. Number of B cells
table(meta$new)/nrow(meta)
plot_df = meta %>% group_by(Cancer) %>% summarise(Cell_count = n()) %>% ungroup()

ggplot(plot_df) +
  geom_bar(aes(x = Cancer,
               y = Cell_count),
           stat = "identity",
           fill = "#93B5C6") +
  xlab("") + 
  ylab("Number of B cells") + 
  cowplot::theme_cowplot() +
  scale_y_continuous(
    trans = "sqrt",
    breaks = c(1000, 10000, 50000, 100000, 200000),
    labels = comma
  ) +
  theme(
    axis.text.x = element_text(size = 8,angle = 45,hjust = 1),
    axis.text.y = element_text(size = 8),
    text = element_text(size = 9),
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3)
  )

ggsave(
  filename = paste0(
    fdir,"Figure1A_Cell_count.pdf"
  ),
  device = "pdf",
  width = 4,
  height = 2.5
)

### 02. Number of BCR seuqneces
tmp = meta %>% filter(CellID %in% BCR$sequence_id)
table(tmp$new)/nrow(tmp)

plot_df = meta %>% filter(CellID %in% BCR$sequence_id) %>% 
  group_by(Cancer) %>% summarise(cell_count = n()) %>% ungroup()

ggplot(plot_df) +
  geom_bar(aes(x = Cancer,
               y = cell_count),
           stat = "identity",
           fill = "#91C788") +
  xlab("") + ylab("Number of BCR seuqneces") + cowplot::theme_cowplot() +
  scale_y_continuous(trans = "sqrt",
                     breaks = c(300, 1000, 5000, 10000, 30000),
                     labels = comma) +
  theme(
    axis.text.x = element_text(size = 8,angle = 45,hjust = 1),
    axis.text.y = element_text(size = 8),
    text = element_text(size = 9),
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    title = element_text(size = 9),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3)
  )

ggsave(
  filename = paste0(fdir,"Figure1A_BCR_source.pdf"),
  device = "pdf",
  width = 3,
  height = 2.5
)

### 03. IHC pie chart

plot_df <- IHC_info %>%
  group_by(`Cancer type`) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = lag(ymax, default = 0),
    labelPosition = (ymax + ymin) / 2,
    label = paste0(`Cancer type`, " (", count, ")")
  )

ggplot(plot_df, aes(
  ymax = ymax,
  ymin = ymin,
  xmax = 4,
  xmin = 3,
  fill = `Cancer type`
)) +
  geom_rect(color = "white",linewidth = 0.2) +
  geom_text(x = 3.5,
            aes(y = labelPosition, label = label),
            size = 1.06) +
  scale_fill_manual(values = Cancer_color_panel) +
  coord_polar(theta = "y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

ggsave(filename = paste0(fdir, "Figure1A_IHC.pdf"),width = 1, height = 1)


## Figure S1A; Number of patients
tmp = meta %>% group_by(Cancer,PatientID,new) %>% summarise()
table(tmp$new)/nrow(tmp)

plot_df = meta %>% group_by(Cancer) %>% summarise(Patient_count = length(unique(PatientID)))

ggplot(plot_df, aes(x = Cancer,
                    y = Patient_count)) +
  geom_bar(stat = "identity",
           fill = "#EFBBCF") +
  xlab("") + ylab("Number of patients") +
  scale_y_continuous(trans = "sqrt",
                     breaks = c(10, 20, 50, 100, 150),
                     labels = comma) +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(size = 5,angle = 45,hjust = 1),
    axis.text.y = element_text(size = 5),
    text = element_text(size = 6),
    axis.line = element_line(size = 0.2),
    axis.ticks = element_line(size = 0.2)
  ) 

ggsave(
  filename = paste0(
    fdir,"FigureS1A_Patient_count.pdf"
  ),
  device = "pdf",
  width = 2.3,
  height = 1.6
)

## Figure S1A; The composition of tissue sources

plot_df = meta %>% group_by(Cancer, Tissue) %>% summarise(n = n())

ggplot(plot_df, aes(Cancer, n, fill = Tissue)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = Tissue_color_panel) +
  xlab("") + ylab("Proportion") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5),
    text = element_text(size = 6),
    axis.line = element_line(size = 0.2),
    axis.ticks = element_line(size = 0.2),
    panel.border = element_rect(colour = "black", linewidth = 0.2),
    legend.position = "None"
  ) +
  labs(fill = 'Tissue') +
  Seurat::RotatedAxis()

ggsave(
  filename = paste0(
    fdir,"FigureS1A_Tissue_source.pdf"
  ),
  device = "pdf",
  width = 2.3,
  height = 1.5
)


## Figure S1B; staining panels heatmap
library(ComplexHeatmap)

m = IHC_info[,c(5:10)]
m[is.na(m)] = 0
m = t(m)

ha = HeatmapAnnotation(
  Cancer = IHC_info$`Cancer type`,
  Tissue = IHC_info$Tissue,
  annotation_name_gp = grid::gpar(fontsize = 6),
  col = list(
    Cancer = Cancer_color_panel,
    Tissue = c("Tumor" = "#5BC0EB", "Adjacent non-tumor" = "#9BC53D")
  ),
  gp = gpar(col = "white", lwd = 0.5)
)

colors = c("0" = "#d1d3d4", "1" = "#fbb040")

pdf(
  paste0(fdir,"FigureS1B.pdf"),
  width = 6,
  height = 3
)

ComplexHeatmap::Heatmap(
  m,
  col = colors,
  top_annotation = ha,
  cluster_rows = F,
  column_split  = IHC_info$`Cancer type` ,
  row_split = c(rep("IHC", 2), rep("mIHC", 4)),
  cluster_columns = F,
  column_names_gp = grid::gpar(fontsize = 5),
  row_names_gp = grid::gpar(fontsize = 5),
  column_title_gp =  grid::gpar(fontsize = 5),
  row_title_gp =  grid::gpar(fontsize = 5),
  width = ncol(m) * unit(0.03, "inch"),
  height = nrow(m) * unit(0.16, "inch"),
  rect_gp = gpar(col = "white", lwd = 0.5),
)

dev.off()

#----- 02. batch correction and annotation -----# 
## Figure S1C; UMAP - The distribution of datasets in the integrated B cell atlas

# batchID color
color = sample(get_palette(palette = c(get_palette("npg", 8),get_palette("Dark2", 8)),k = 54), 54)
color = read_rds("batchID_color.rds")

ggplot(meta, aes(x = UMAP1, y = UMAP2, color = batchID)) +
  geom_point(size = 0.03,
             shape = 16,
             stroke = 0) +
  theme_void()+
  scale_color_manual(values = color, name = '')+
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(
  paste0(fdir, "FigureS1C.png"),
  width = 2,
  height = 2,
  dpi = 300
)

ggplot(meta, aes(x = UMAP1, y = UMAP2, color = batchID)) +
  geom_point(size = 0,
             shape = 16,
             stroke = 0) +
  theme_void()+
  scale_color_manual(values = color, name = '') +
  theme(
    aspect.ratio = 1,
    legend.position = 'bottom',
    plot.margin = margin(0,0,0,0)
  ) +
  guides(color = guide_legend(
    ncol = 2,
    override.aes = list(size = 2, alpha = 1)
  )) +
  theme(legend.text = element_text(size = 5),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.height = unit(0,"cm"),
        legend.box.spacing = unit(0, 'cm'))

ggsave(
  paste0(fdir, "FigureS1C_legend.pdf"),
  width = 3.5,
  height = 6,
  dpi = 500
)

## Figure 1B; UMAP - The major clusters and subsets of B cells
ggplot(meta, aes(x = UMAP1, y = UMAP2, color = Annotation_major)) +
  geom_point(size = 0.01,
             shape = 16,
             stroke = 0) +
  theme_void()+
  scale_color_manual(values = Major_color_panel, name = '')+
  theme(aspect.ratio = 1,
        legend.position = "none")

ggsave(
  paste0(fdir, "Figure1B_major.png"),
  width = 2,
  height = 2,
  dpi = 300
)

## Figure S1D; UMAP - The expression patterns of major cluster marker genes 
color = RColorBrewer::brewer.pal(9, "Blues")

markers = c("MS4A1","IGHD","IGHM","TCL1A",
            "CD27","TNFRSF13B","FCRL4",
            "AICDA","BCL6","RGS13","MKI67",
            "MZB1","XBP1","CD38",
            "IGHG1","IGHA1")

plot_df = data.frame(
  UMAP1 = obj$UMAP1,
  UMAP2 = obj$UMAP2,
  FetchData(object = obj, vars = markers)
) %>%
  tidyr::pivot_longer(!c(UMAP1, UMAP2), names_to = 'Markers', values_to = 'Expr')

plot_df$Markers = factor(plot_df$Markers, level = markers)

ggplot(plot_df) + geom_point(
  aes(x = UMAP1, y = UMAP2, color = Expr),
  size = .05,
  stroke = 0,
  shape = 16
) +
  theme_void() +
  theme(
    aspect.ratio = 1,
    plot.margin = unit(c(-2,-2,-2,-2), "char"),
    text = element_text(family = "Arial")
  ) +
  facet_wrap( ~ Markers, scale = 'free', ncol = 4) +
  theme(strip.text = element_text(size = 7))+
  labs(color = "Exp") +
  patchwork::plot_layout(guides = 'collect') &
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.position = 'bottom',
    legend.key.size = unit(0.15, "inch"),
    legend.box.margin = margin(-10, 0, 0, 0)
  ) &
  scale_color_gradientn(colours = color[c(1, 2, 5, 6 , 7, 8, 9)],
                        limits = c(0, 4),
                        oob = scales::squish)

ggsave(filename = paste0(fdir,"FigureS1D.png"),
       device = "png",
       width = 4.2,
       height = 5.5)



#----- 03. BCR associated analysis -----#
BCR$SHM_rate = (1-BCR$v_identity) * 100
BCR$CellID = BCR$sequence_id
tmp = BCR[,c("CellID","clone_id","clonal","c_call","SHM_rate")]
BCR_df = tmp %>% left_join(meta)

##--- Figure 1D; Distribution of IgH isotypes (left) across B cell subsets
plot_df = BCR_df  %>%
  filter(!is.na(c_call),
         Tissue == "Tumor") %>%
  count(Annotation, c_call) %>%
  mutate(Annotation = forcats::fct_rev(as.factor(Annotation)))

plot_df$c_call = factor(
  plot_df$c_call,
  levels = c(
    "IGHE",
    "IGHM",
    "IGHD",
    "IGHA1",
    "IGHA2",
    "IGHG1",
    "IGHG2",
    "IGHG3",
    "IGHG4"
  )
)

p1 = ggplot(plot_df, aes(Annotation, n, fill = c_call)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(
    values =
      c(
        "IGHE" = "#CED6C3",
        "IGHM" = "#98C9DD",
        "IGHD" = "#207CB5",
        "IGHA1" = "#A6D38E",
        "IGHA2" = "#37A849",
        "IGHG1" = "#F69595",
        "IGHG2" = "#EB2A2A",
        "IGHG3" = "#FCBA71",
        "IGHG4" = "#f78200"
      )
  ) +
  labs(x = "", y = "Proportion", fill = 'BCR isotype') +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1,vjust = 0.5),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 8),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    panel.border = element_rect(colour = "black", linewidth = 0.3),
    legend.position = "none"
  ) + coord_flip()
p1

plot_df = BCR_df %>% 
  filter(Cancer %in% c("CRC","HCC","ICC") == FALSE, Tissue == "Tumor")%>%
  mutate(Annotation = forcats::fct_rev(as.factor(Annotation)))

plot_df$SHM_levels = ifelse(plot_df$SHM_rate < 1,"Low","Median")
plot_df$SHM_levels = ifelse(plot_df$SHM_rate > 5,"High",plot_df$SHM_levels)
plot_df$SHM_levels = factor(plot_df$SHM_levels,levels = c("Low","Median","High"))

plot_df = plot_df %>% group_by(Annotation,SHM_levels) %>% summarise(n = n()) %>%
  group_by(Annotation) %>% mutate(sum_n = sum(n)) %>% 
  ungroup() %>% mutate(percent = n/sum_n)

p2 = ggplot(plot_df,
           aes(x = Annotation, y = percent, fill = SHM_levels)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "High" = "#78290f",
    "Median" = "#ff7d00",
    "Low" = "#ffecd1"
  )) +
  labs(x = "", y = "Proportion", fill = 'SHM level') +
  theme_bw() +
  coord_flip()+
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 8),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    panel.border = element_rect(colour = "black", linewidth = 0.3),
    legend.position = "none"
  ) 
p2

##--- 
plot_df = BCR_df %>% 
  filter(Cancer %in% c("CRC","HCC","ICC") == FALSE, Tissue == "Tumor") %>%
  mutate(Annotation = forcats::fct_rev(as.factor(Annotation)))

plot_df = plot_df %>% group_by(Annotation) %>% 
  summarise(median_SHM = median(SHM_rate))

p3 = ggplot(plot_df,
       aes(x = Annotation, y = median_SHM)) +
  geom_point(shape = 16, stroke = 0, size = 1) +
  labs(x = "", y = "SHM rate") +
  theme_bw() +
  coord_flip()+
  theme(
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_blank(),
    text = element_text(size = 8),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    axis.ticks.y  = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.3),
    legend.position = "none"
  ) 

p3
cowplot::plot_grid(p1, p2, p3,nrow = 1,align = "h", rel_widths = c(1, 1, 0.4))

ggsave(paste0(fdir, "Figure1D.pdf"),width = 5,height = 2.7)


##--- Figure S2B; Heatmap showing the median SHM rate for each TIB major lineage 
tmp = BCR_df %>%
  filter(Cancer %in% c("CRC", "HCC", "ICC") == FALSE) %>%
  filter(Tissue == "Tumor") %>%
  group_by(Cancer) %>%
  mutate(sample_n = length(unique(SampleID))) %>%
  mutate(Cancer_new = paste0(Cancer," (N = ", sample_n,")")) %>%
  group_by(Cancer_new, Annotation_major_2) %>%
  dplyr::summarise(avg_SHM = median(SHM_rate), B_n = n())

colnames(tmp)[1] = "Cancer"  

tmp = tidyr::spread(data = tmp[, c("Cancer", "Annotation_major_2", "avg_SHM")], key = Annotation_major_2, value = avg_SHM)

matrix = tmp %>% tibble::column_to_rownames("Cancer")
matrix = matrix[,c("Bn","Bgc","Bm","ASC")]
col_fun <- circlize::colorRamp2(c(0, 5, 10),
                                c("#0F7B9F", "white", "#D83215"))
library(ComplexHeatmap)
pdf(
  paste0(fdir,"FigureS2B.pdf"),
  width = 3,
  height = 3
)

ComplexHeatmap::Heatmap(
  matrix,
  col = col_fun,
  na_col = "grey",
  cluster_rows = F,
  cluster_columns = F,
  row_names_side = "left",
  column_names_gp = grid::gpar(fontsize = 7),
  row_names_gp = grid::gpar(fontsize = 7),
  rect_gp = gpar(col = "white", lwd = 1),
  width = ncol(matrix) * unit(0.2, "inch"),
  height = nrow(matrix) * unit(0.2, "inch"),
  name = "Median SHM rate (%)",
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 7),
    title_gp = gpar(fontsize = 8)
  )
)

dev.off()


##--- Figure S2C; Differentially expressed genes between the SHM-high (red) and SHM-low (blue) groups
plot_df = BCR_df %>%
  filter(Cancer %in% c("CRC", "HCC", "ICC") == FALSE) %>% 
  group_by(Annotation_major_2) %>%
  mutate(median_SHM = median(SHM_rate)) %>% 
  mutate(group = ifelse(SHM_rate > median_SHM, "High","Low")) %>%
  as.data.frame()

sub_obj = subset(obj,CellID %in% plot_df$CellID)

rownames(plot_df) = plot_df$CellID
plot_df = plot_df[colnames(sub_obj),]
sub_obj$group = plot_df$group
Idents(sub_obj) = sub_obj$group

# perform DE
i = sort(unique(sub_obj$Annotation))[1]

DE_sum = data.frame()

for (i in sort(unique(sub_obj$Annotation_major_2))) {
  tmp_obj = subset(sub_obj, Annotation_major_2 == i)
  DE = Seurat::FindMarkers(
    tmp_obj,
    ident.1 = "High",
    ident.2 = "Low",
    logfc.threshold = 0
  )
  DE$gene = rownames(DE)
  DE$annotation = i
  DE_sum = rbind(DE_sum, DE)
}

FC_cutoff <- 0.3
p_cutoff <- 0.05

plot_df <- DE_sum %>% dplyr::filter(annotation != "Bn")
plot_df$pct_value = abs(plot_df$pct.1-plot_df$pct.2)

plot_df$DE <- "Not significant"
plot_df$DE[plot_df$avg_log2FC > FC_cutoff & plot_df$p_val_adj < p_cutoff ] <- "SHM-high group"
plot_df$DE[plot_df$avg_log2FC < -FC_cutoff & plot_df$p_val_adj < p_cutoff] <- "SHM-low group"

plot_df$label <-
  ifelse(
    abs(plot_df$avg_log2FC) > 0.3 &
      plot_df$p_val_adj < 0.001 & plot_df$pct_value > 0.1,
    plot_df$gene,
    NA
  )

mycolors <- c("#0F7B9F", "#D83215", "grey")
names(mycolors) <- c("SHM-low group", "SHM-high group", "Not significant")
plot_df$annotation = factor(plot_df$annotation,levels = c("Bgc","Bm","ASC"))

ggplot(data=plot_df, aes(x=avg_log2FC, y=-log10(p_val_adj), color = DE, label = label)) +
  geom_point(size = 0.2) +
  scale_colour_manual(values = mycolors, name = paste0("Enriched in")) +
  geom_vline(xintercept = c(FC_cutoff, -FC_cutoff), linetype="dashed", color = "black") +
  geom_hline(yintercept = -log10(p_cutoff), linetype="dashed", color = "black") +
  ggrepel::geom_text_repel(size = 5 * 0.35, show.legend = FALSE, max.overlaps = 40) +
  cowplot::theme_cowplot() +
  theme(text = element_text(size = 8, family = "ArialMT"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = "right",
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)) +
  xlab("Log2(fold change)") +
  ylab("-Log10(q-value)") + 
  facet_wrap(~annotation, scales = "free")

ggsave(paste0(fdir,"FigureS2B.pdf"),
       width = 8.2,
       height = 2.5)


#----- 04. B cell expression analysis ------# 
##--- Figure S1F; Expression of representative signature genes across B cell subsets. IL35 expression was defined as the concurrent expression of its subunits IL12A and EBI3

obj$Annotation = as.factor(obj$Annotation)

obj$IL12A = FetchData(obj,"IL12A")
obj$EBI3= FetchData(obj,"EBI3")
obj$IL35 = obj$IL12A * obj$EBI3
obj$IL35 = ifelse(obj$IL35 > 0,1,0)

gene_list = list(
  c("FCER2", "TCL1A", "IL4R", "CD72", "BACH2", "IGHD", "IGHM"),
  c("NR4A1", 'NR4A2', "CREM", "CD83"),
  c("ISG15", "IFI44L", "IFI6", "IFIT3"),
  
  c("CD27", "TNFRSF13B","TXNIP", "GPR183"),
  c("HSPA1A", "HSPA1B","DNAJB1", "EGR1"),
  c("FCRL4","FCRL5", "ITGAX", "TBX21", "CR2"),
  c("CCR1", "CXCR3", "PDCD1", "HCK", "FCRL3", "FGR"),
  c("NME1", "APEX1", "POLD2", "POLE3","MYC"),
  
  c("BCL6", "RGS13", "AICDA","IL21R"),
  c("MKI67","STMN1","HMGB2","TOP2A"),
  
  c("CD38","MZB1","PRDM1","IRF4","XBP1"),
  c("MS4A1","LTB","HLA-DRA", "HLA-DRB1","HLA-DPA1","HLA-DQA1"),
  c("IGHG1", "IGHG2", "IGHG3", "IGHG4","IGHA1", "IGHA2"),
  c("IL10","IL12A","EBI3","TGFB1"),
  c("IL35")
)

p = yy_Dotplot(seuratObj = obj,
                genes = gene_list,
                group.by = "Annotation")
p +
  theme(
    panel.spacing = unit(x = 0.2, units = "lines"),
    strip.background = element_blank(),
    text = element_text(size = 0),
    panel.grid = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    )
  ) + theme(
    axis.text.x = element_text(
      size = 6,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 6),
  ) + xlab("") + ylab("") +
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.position = 'bottom',
    legend.key.size = unit(0.15, "inch"),
    plot.margin = unit(c(0, 1, 0, 0), "char")
  ) +
  scale_radius(limits = c(0, 100), range = c(0, 2.5))

ggsave(paste0(fdir,"FigureS1F.pdf"),width = 8.5,height = 4)

##--- Figure S1G; Heatmap showing the expression of genes mechanistically linked with CSR
gene_list = c("APEX1","XRCC5","XRCC6","POLD2","POLE3",
              "NCL","NME2","DDX21",
              "NPM1","SERBP1",
              "MIR155HG","HSP90AB1",
              "BATF","HIVEP3","BHLHE40","IRF4")

obj$annotation_short = str_extract(obj$Annotation,pattern = "c[0-9]{2}")

plot_df = FetchData(obj,vars = c(gene_list,"annotation_short"),slot = "data")
plot_df = plot_df %>% group_by(annotation_short) %>% summarise_all(function(x) mean(x = expm1(x = x)))
plot_df = plot_df %>% tibble::column_to_rownames(var = "annotation_short")
plot_df = scale(plot_df)
range(plot_df)

col_fun <- circlize::colorRamp2(
  c(-1.5, 0, 3.5),
  c("#0F7B9F", "white", "#D83215")
)


pdf(paste0(fdir,"FigureS1G.pdf"), width = 3,height = 4)

ComplexHeatmap::Heatmap(
  plot_df,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "bottom",
  column_names_rot = 90,
  name = "Z-score",
  heatmap_legend_param = list(
    legend_direction = "vertical",
    title_position = "topcenter",
    legend_width = unit(0.05, "inch"),
    title_gp = gpar(fontsize = 7),
    labels_gp = gpar(fontsize = 6)
  ),
  column_names_gp = grid::gpar(fontsize = 6),
  row_names_gp = grid::gpar(fontsize = 6),
  width = ncol(plot_df) * unit(0.08, "inch"),
  height = nrow(plot_df) * unit(0.08, "inch")
)

dev.off()

##--- Figure S1H; Proportion of immunoglobulin gene count across B cell subsets.
obj$annotation_short = as.factor(obj$annotation_short)
Idents(obj) = obj$annotation_short

features <- grep(pattern = "^IG", x = rownames(obj), value = TRUE) 
percent.featureset <- colSums(x = GetAssayData(object = obj, slot = "counts")[features, , drop = FALSE])/ 
  colSums(x = obj@assays$RNA@counts) * 100 

obj$percent.ig = percent.featureset

Subset_color_panel_short = Subset_color_panel
names(Subset_color_panel_short) = str_extract(names(Subset_color_panel),pattern = "c[0-9]{2}")

p = Seurat::VlnPlot(object = obj,
                features = "percent.ig",
                pt.size = 0) +
  scale_fill_manual(values = Subset_color_panel_short) +
  theme(legend.position = "none") +
  labs(title = "The proportion of immunoglobulin gene counts", x = "", y = "") +
  theme(
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 7),
    plot.title = element_text(size = 7,face = "plain"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    plot.margin = unit(c(0, 1, 0, 0), "char"),
  )

p$layers[[1]]$aes_params$size = 0.1
p

ggsave(paste0(fdir, "FigureS1H.pdf"),width = 2.5,height = 1.5)

##--- Figure 1D; Heatmap showing expression of gene signatures in B cell subsets.
obj_withoutIG = read_rds("obj/All_obj_withoutIG.rds")

library("msigdbr")
library("clusterProfiler")
packageVersion("msigdbr")

gene_sets = msigdbr(species = "Homo sapiens")
gene_sets$gs_GO = str_replace_all(gene_sets$gs_url,"http://amigo.geneontology.org/amigo/term/","")

gene_sets_df = gene_sets %>% group_by(gs_name,gs_description,gs_GO) %>% summarise()

df <- read_tsv("gene_signatures.tsv") %>% as.data.frame()
df <- df %>% tibble::column_to_rownames(var = "index")

addm_List = list()

for (i in colnames(df)){
  addm_List[[i]] = df[,i][!is.na(df[,i])]
}

obj_withoutIG = Seurat::AddModuleScore(object = obj_withoutIG,features = addm_List,search = TRUE)

n = length(addm_List)
plot_df = FetchData(obj,vars = paste("Cluster", 1:n, sep = ""))
plot_df = FetchData(obj_withoutIG,vars = paste("Cluster", 1:n, sep = ""))

colnames(plot_df) = names(addm_List)

plot_df$annotation = obj$Annotation
plot_df = plot_df %>% group_by(annotation) %>% summarise_all(.funs = mean)
plot_df = plot_df %>% column_to_rownames(var = "annotation")
plot_df = scale(plot_df)

col_fun <- circlize::colorRamp2(
  c(-2, 0, 3),
  c("#0F7B9F", "white", "#D83215")
)

plot_df = t(plot_df)

pdf(paste0(fdir,"Figure1D.pdf"), width = 10,height = 12)
ComplexHeatmap::Heatmap(
  plot_df,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "right",
  column_names_side = "bottom",
  column_names_rot = 90,
  name = "Z-score",
  heatmap_legend_param = list(
    legend_direction = "vertical",
    title_position = "topcenter",
    legend_width = unit(0.8, "inch"),
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 7)
  ),
  column_names_gp = grid::gpar(fontsize = 7),
  row_names_gp = grid::gpar(fontsize = 7),
  width = ncol(plot_df) * unit(0.12, "inch"),
  height = nrow(plot_df) * unit(0.12, "inch")
)
dev.off()

#----- 05. spatial -----#
df = read_tsv("TLS_mIHC.tsv")

plot_df =  tidyr::gather(data = df,key = annotation, value = Percent, 3:6)
plot_df$Percent = plot_df$Percent*100
plot_df$annotation = factor(plot_df$annotation,levels = c("Bn","Bm","Bgc","ASC"))

ggplot(plot_df, aes(x = annotation, y = Percent)) +
  geom_boxplot(outlier.color = NA,
               lwd = 0.3) +
  geom_jitter(aes(color = Cancer),size = 0.2,
              width = 0.2) +
  scale_color_manual(values = Cancer_color_panel) +
  xlab("") + ylab("Percentage of cells within TLS") +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    text = element_text(size = 7),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    legend.position = "bottom",
    plot.margin = unit(c(1, 1, 1, 1), "char")
  ) +
  RotatedAxis()

ggsave(paste0(fdir,"Figure1E.pdf"),width = 2, height = 2.5)
