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

## III. load color
Tissue_color_panel = c(
  "Tumor" = "#5BC0EB",
  "Adjacent non-tumor tissue" = "#9BC53D",
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

## IV. load code


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
    Tissue = c("Tumor" = "#5BC0EB", "Adjacent non-tumor tissue" = "#9BC53D")
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
