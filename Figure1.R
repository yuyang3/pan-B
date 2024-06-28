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


## Figure S1B
