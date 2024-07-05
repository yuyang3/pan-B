# Figure5.R

fdir = ""
## I. load packages
library(readr)
library(dplyr)
library(ggplot2)
library(scales)

## II. load data
obj = read_rds("ss2_data/obj.rds")
meta = obj@meta.data

#----- smart-seq2 B cell atlas -----#
## Figure 5D; UMAP - distribution of datasets

length(unique(meta$Dataset))

ggplot(meta, aes(x = UMAP1, y = UMAP2, color = Dataset)) +
  geom_point(size = 0.2,
             shape = 16,
             stroke = 0) +
  theme_void()+
  scale_color_manual(values = get_palette("nejm", 8) , name = 'Dataset')+
  theme(aspect.ratio = 1,
        legend.position = "right",
        plot.margin = margin(0,0,0,0))+
  guides(color = guide_legend(
    ncol = 1,
    override.aes = list(size = 2, alpha = 1)
  )) +
  theme(legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.height = unit(0.3,"cm"),
        legend.box.spacing = unit(0.05, 'cm'))

ggsave(
  paste0(fdir, "Figure5D.pdf"),
  width = 3,
  height = 1.3,
  dpi = 300
)

## Figure 5E; UMAP - distribution of major clusters
meta$Annotation %>% unique() %>% sort()

ggplot(meta, aes(x = UMAP1, y = UMAP2, color = Annotation)) +
  geom_point(size = 0.2,
             shape = 16,
             stroke = 0) +
  theme_void()+
  scale_color_manual(values = get_palette("nejm", 5), name = '')+
  theme(aspect.ratio = 1,
        legend.position = "right",
        plot.margin = margin(0,0,0,0))+
  guides(color = guide_legend(
    ncol = 1,
    override.aes = list(size = 2, alpha = 1)
  )) +
  theme(legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.height = unit(0.3,"cm"),
        legend.box.spacing = unit(0.05, 'cm'))

ggsave(
  paste0(fdir, "Figure5E.pdf"),
  width = 3,
  height = 1.3,
  dpi = 300
)

## Figure 5F; Expression of signature genes across B cell clusters in the SMRAT-seq2 atlas.
# Define gene list
gene_list <- list(
  c("MZB1", "XBP1", "CD38", "SDC1", "IGHG1", "IGHA1"),
  c("CD19", "MS4A1", "IGHM", "IGHD"),
  c("CD27", "TNFRSF13B", "TXNIP", "GPR183"),
  c("MKI67", "STMN1", "HMGB2", "BIRC5"),
  c("BCL6", "AICDA", "RGS13", "TCL1A"),
  c("IL10", "IL12A", "EBI3", "TGFB1", "CD274"),
  c("IL35")
)

obj$Annotation_res.0.3 = str_extract(string = obj$Annotation, pattern = "c[1-5]")

# Generate dot plot
p <- yy_Dotplot(
  seuratObj = obj,
  genes = gene_list,
  group.by = "Annotation_res.0.3"
)

# Customize plot
p <- p +
  theme(
    panel.spacing = unit(0.2, "lines"),
    strip.background = element_blank(),
    text = element_text(size = 0),
    panel.grid = element_line(colour = "grey", linetype = "dashed", size = 0.2),
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.position = 'right',
    legend.key.size = unit(0.15, "inch"),
    plot.margin = unit(c(0, 0, 0, 0), "char")
  ) +
  labs(x = "", y = "") +
  scale_radius(limits = c(0, 100), range = c(0, 2.5))

p

ggsave(paste0(fdir,"Figure5F.pdf"),width = 4.6,height = 1.8)


## Figure 5G; UMAP - the expression pattern of IL10
plot_df = data.frame(
  UMAP1 = meta$UMAP1,
  UMAP2 = meta$UMAP2,
  FetchData(object = obj, vars = 'IL10')
) %>%
  tidyr::pivot_longer(!c(UMAP1, UMAP2), names_to = 'Markers', values_to = 'Expr')

sum(plot_df$Expr > 0)/nrow(plot_df)

ggplot(plot_df) + geom_point(
  aes(x = UMAP1, y = UMAP2, color = Expr),
  size = .2,
  stroke = 0,
  shape = 16
) +
  theme_void() +
  theme(
    aspect.ratio = 1,
    plot.margin = unit(c(0,0,0,0), "char")) +
  labs(title = "IL10 (15.91%)", color = "Exp") +
  theme(
    plot.title = element_text(size = 8, vjust = 0, hjust = 0.5), 
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.position = 'right',
    legend.key.size = unit(0.15, "inch"),
    legend.box.margin = margin(-10, 0, 0, 0)
  ) +
  scale_color_gradient2(mid = 'lightgrey', high = 'blue')

ggsave(
  paste0(fdir, "Figure5G.pdf"),
  width = 3,
  height = 1.5,
  dpi = 300
)

## Figure 5H; IL10 expression in B cells from the SMART-seq2 atlas sekecred by in-silco FACS according to reported Breg surface markers
# CD19+CD24+CD38+
matrix <- t(as.matrix(obj@assays$RNA@data[c("CD19", "CD24", "CD38"), ]))
sub <- filter(as.data.frame(matrix), CD19 > 0 & CD24 > 0 & CD38 > 0)

obj$group <- ifelse(rownames(matrix) %in% rownames(sub), "CD19+CD24+CD38+", "else")
Idents(obj) <- obj$group

VlnPlot(obj, c("IL10"), group.by = "group", ncol = 1, pt.size = 0) +
  labs(title = "", x= "", y = "Expression of IL10")+
  theme(legend.position = "none",
        axis.text = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        plot.margin = unit(c(0,0,0,0), "char"),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(size = 0.3)) +
  scale_fill_manual(values = c("CD19+CD24+CD38+" = "#e7bf41", "else" = "#2e72ba"))
ggsave(filename = paste0(fdir,"Figure5H_L.pdf"),width = 1.2, height = 2.5)

# CD19+CD24+CD27+
matrix <- t(as.matrix(obj@assays$RNA@data[c("CD19", "CD24", "CD27"), ]))
sub <- filter(as.data.frame(matrix), CD19 > 0 & CD24 > 0 & CD27 > 0)

obj$group <- ifelse(rownames(matrix) %in% rownames(sub), "CD19+CD24+CD27+", "else")
Idents(obj) <- obj$group

VlnPlot(obj, c("IL10"), group.by = "group", ncol = 1, pt.size = 0) +
  labs(title = "", x= "", y = "Expression of IL10")+
  theme(legend.position = "none",
        axis.text = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        plot.margin = unit(c(0,0,0,0), "char"),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(size = 0.3)) +
  scale_fill_manual(values = c("CD19+CD24+CD27+" = "#e7bf41", "else" = "#2e72ba"))
ggsave(filename = paste0(fdir,"Figure5H_R.pdf"),width = 1.2, height = 2.5)


## Figure 5I; Heatmap showing major lineage markers expression on IL10+ B cells from the Smart-seq2 atlas. Each column represents an IL10+ B cell
library(ComplexHeatmap)

subset = obj@assays$RNA@counts["IL10",] > 0
sum(subset)
tmp_obj = obj[, subset]

gene_list = c(
  c("CD19", "MS4A1"),
  c("RGS13", "BCL6"),
  c("MKI67", "TOP2A", "BIRC5"),
  c("MZB1", "XBP1", "PRDM1")
)

row_split = c(rep("1", 2),
              rep("2", 2),
              rep("3", 3),
              rep("4", 3))

matrix <- scale(t(as.matrix(tmp_obj@assays$RNA@data[gene_list, tmp_obj$CellID])))

df <- tmp_obj@meta.data %>%
  select(CellID, Annotation) %>%
  left_join(data.frame(Annotation = c("c2_Bn", "c3_Bm", "c5_Bgc", "c4_Bcycling", "c1_ASC"), 
                       anno_order = 1:5), by = "Annotation") %>%
  arrange(anno_order) %>%
  mutate(col_split = case_when(
    Annotation %in% c("c2_Bn", "c3_Bm") ~ 1,
    Annotation == "c5_Bgc" ~ 2,
    Annotation == "c4_Bcycling" ~ 3,
    TRUE ~ 4
  ))

col_split <- df$col_split
matrix <- t(matrix[df$CellID, ])

pdf(
  paste0(fdir,"Figure5I.pdf"),
  width = 4,
  height = 2
)

col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#0F7B9F", "white", "#D83215"))

ha = HeatmapAnnotation(anno = df$col_split,
                       col = list(anno = c("1" = "#E18727FF", "2" = "#E60012", 
                                           "3" = "#009FB9","4" = "#009944")))

ComplexHeatmap::Heatmap(
  top_annotation = ha,
  matrix,
  col = col_fun,
  row_split = row_split,
  column_split  = col_split,
  cluster_rows = F,
  cluster_columns = F,
  show_column_names = F,
  row_names_gp = grid::gpar(fontsize = 7),
  width = ncol(matrix) * unit(0.003, "inch"),
  height = nrow(matrix) * unit(0.15, "inch"),
  name = "Expression"
)

dev.off()
