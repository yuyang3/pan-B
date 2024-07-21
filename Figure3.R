# Figure3.R

fdir = ""
## load packages
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

## load data
obj = read_rds("obj/All_obj.rds")
meta = obj@meta.data
BCR = read_rds("BCR.rds")

## load color 
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

#----- 01. sub type abundance data frame -----#
if(TRUE){
  df.sample_count <- meta %>%
    group_by(SampleID, PatientID, Cancer, Tissue, Treatment_status) %>%
    summarise(B_n = n(), .groups = 'drop')
  
  for (i in sort(unique(meta$Annotation))) {
    tmp = meta %>%
      filter(Annotation == i) %>%
      group_by(SampleID) %>%
      summarise(n = n())
    
    tmp[is.na(tmp)] = 0
    
    df.sample_count = left_join(df.sample_count, tmp, by = "SampleID")
    df.sample_count[is.na(df.sample_count)] = 0
    
    df.sample_count$percent =
      round(df.sample_count$n / df.sample_count$B_n * 100, 2)
    
    i_names = paste0(i, "_n")
    colnames(df.sample_count)[length(df.sample_count) - 1] = i_names
    
    i_names = paste0(i, "_percent")
    colnames(df.sample_count)[length(df.sample_count)] = i_names
  }
  
  df.sample_count_flt = df.sample_count %>% filter(B_n > 50, Treatment_status == "treatment naïve")
}

##--- Figure 3B; Pearson correlation between the abundance of c12 and c13 Bgc cells in tumor samples with B cells > 50.
plot_df = as.data.frame(df.sample_count_flt) %>% filter(Tissue == "Tumor")

plot_df = plot_df[, c("Cancer",
                      "c12_Bgc_LZ-like_percent",
                      "c13_Bgc_DZ-like_percent")]

colnames(plot_df)[2:3] = c("x", "y")
plot_df$Cancer = as.character(plot_df$Cancer)
max(plot_df$y)

ggplot(plot_df, aes(x = x, y = y)) +
  geom_point(aes(color = Cancer), size = 0.5, shape = 16, stroke = 0) +
  geom_smooth(
    method = 'lm',
    se = F,
    color = 'red',
    size = 0.3
  ) +
  ggpubr::stat_cor(method = "pearson",
                   output.type = "text",
                   size = 7 * 0.35) +
  scale_color_manual(values = Cancer_color_panel, name = "Cancer") +
  xlab("c12_Bgc_LZ-like") + ylab("c13_Bgc_DZ-like") +
  cowplot::theme_cowplot() +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    text = element_text(size = 8),
    plot.margin = unit(c(0, 0, 0, 0), "char"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    legend.position = "right"
  ) +
  xlim(-0.1, 42)+
  ylim(-0.1, 42)+
  guides(color = guide_legend(
    ncol = 2,
    override.aes = list(size = 2, alpha = 1)
  ))+
  theme(legend.text = element_text(size = 7),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.height = unit(0.4,"cm"),
        legend.box.spacing = unit(0, 'cm'))

ggsave(paste0(fdir, "Figure3B.pdf"),
       width = 2,
       height = 2)

##--- Figure S4A and S4B; Boxplots showing the proportions of c12 and c13 Bgc cells among TIBs across cancer types. 
i = "c12_Bgc_LZ-like"
i = "c13_Bgc_DZ-like"
if(TRUE){
  i_names = paste0(i, "_percent")
  tmp = df.sample_count_flt[, c("SampleID", "Cancer", "Tissue", eval(i_names))]
  colnames(tmp)[length(tmp)] = "percent"
  
  tmp2 = tmp %>% 
    filter(Tissue %in% "Tumor") %>%
    group_by(Cancer) %>% 
    mutate(order_value = median(percent))
  
  a = summary(aov(percent ~ Cancer, data = tmp2))
  p_value = a[[1]][["Pr(>F)"]][1]
  F_value = a[[1]][["F value"]][1]
  
  p1 = ggplot(tmp2, aes(
    x = reorder(Cancer, order_value),
    y = percent,
    color = Cancer
  )) +
    geom_boxplot(outlier.color = NA, lwd = 0.3) +
    geom_jitter(size = 0.2, width = 0.2) +
    scale_color_manual(values = Cancer_color_panel) +
    xlab("") + ylab(expression(paste("% among TIBs"))) +
    ggtitle(paste0(i),
            paste0("P = ", signif(p_value, 2), ", F = ", round(F_value, 2))) +
    cowplot::theme_cowplot() +
    theme(
      axis.text.y = element_text(size = 7),
      text = element_text(size = 8),
      plot.title = element_text(size = 8,hjust = 0.5,face = 'plain'),
      plot.subtitle = element_text(size = 7),
      axis.line = element_line(size = 0.3),
      axis.ticks = element_line(size = 0.3),
      legend.position = "None",
      plot.margin = unit(c(0,0,-1,0), "char")
    ) +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank()
    )
  
  p2 <- ggplot(tmp2 %>% group_by(Cancer, order_value) %>%
                 mutate(fill = ifelse(percent > 0 , "non-zero", "zero"))) +
    geom_bar(
      aes(x = reorder(Cancer, order_value), fill = fill),
      stat = "count",
      position = "fill",
      width = 0.8
    ) +
    scale_fill_manual(values = c("non-zero" = "#d83215", "zero" = "lightgrey")) +
    xlab("") + ylab("Proportion") +
    cowplot::theme_cowplot() +
    theme(
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      text = element_text(size = 8),
      axis.line = element_line(size = 0.3),
      axis.ticks = element_line(size = 0.3),
      legend.position = "none",
      plot.margin = unit(c(0,0,0,0), "char")
    )
  
  p2
  
  
  # Remove legend for the bar plot
  cowplot::plot_grid(p1, p2, ncol = 1,align = "v", rel_heights = c(1, 1))
  
  ggsave(
    filename = paste0(fdir, "/FigureS4AB_", i, ".pdf"),
    device = "pdf",
    width = 3.7,
    height = 2.2
  )
}

##--- Figure 3F; Boxplots comparing the abundances of Bn cells, Bm cells and ASCs between tumors with and without Bgc cells.
if(TRUE) {
  df.sample_count <- meta %>%
    group_by(SampleID, PatientID, Cancer, Tissue, Treatment_status) %>%
    summarise(B_n = n(), .groups = 'drop')
  
  for (i in sort(unique(meta$Annotation))) {
    tmp = meta %>%
      filter(Annotation == i) %>%
      group_by(SampleID) %>%
      summarise(n = n())
    
    tmp[is.na(tmp)] = 0
    
    df.sample_count = left_join(df.sample_count, tmp, by = "SampleID")
    df.sample_count[is.na(df.sample_count)] = 0
    
    df.sample_count$percent =
      round(df.sample_count$n / df.sample_count$B_n * 100, 2)
    
    i_names = paste0(i, "_n")
    colnames(df.sample_count)[length(df.sample_count) - 1] = i_names
    
    i_names = paste0(i, "_percent")
    colnames(df.sample_count)[length(df.sample_count)] = i_names
  }
  
  for (i in sort(unique(meta$Annotation_major_2))) {
    tmp = meta %>%
      filter(Annotation_major_2 == i) %>%
      group_by(SampleID) %>%
      summarise(n = n())
    
    tmp[is.na(tmp)] = 0
    
    df.sample_count = left_join(df.sample_count, tmp, by = "SampleID")
    df.sample_count[is.na(df.sample_count)] = 0
    
    df.sample_count$percent =
      round(df.sample_count$n / df.sample_count$B_n * 100, 2)
    
    i_names = paste0(i, "_n")
    colnames(df.sample_count)[length(df.sample_count) - 1] = i_names
    
    i_names = paste0(i, "_percent")
    colnames(df.sample_count)[length(df.sample_count)] = i_names
  }
}

all_meta = read_rds("all_meta_0707.rds")

CD45_meta = all_meta %>%
  filter(
    majortype %in% c("B", "Myeloid", "T&NK"),
    `Isolated_or_sorted_cell_population` %in% c(
      "CD45+",
      "CD45+/CD45-(1:1 mix)",
      "CD45+/CD45-(3:1 mix)",
      "CD45+&EPCAM-",
      "EPCAM-",
      "None",
      "None, CD45+"
    )
  )

df.sample_count_CD45 = CD45_meta %>% 
  group_by(SampleID) %>% 
  summarise(CD45_n=n())

plot_df = df.sample_count %>%
  left_join(df.sample_count_CD45, by = "SampleID") %>%
  filter(B_n > 50, Treatment_status == "treatment naïve") %>%
  filter(Tissue == "Tumor") %>%
  filter(is.na(CD45_n) == FALSE)

plot_df$N_percent = plot_df$Bn_n / plot_df$CD45_n * 100
plot_df$M_percent = plot_df$Bm_n / plot_df$CD45_n * 100
plot_df$A_percent = plot_df$ASC_n / plot_df$CD45_n * 100

plot_df$GCB_percent = plot_df$`c12_Bgc_LZ-like_percent` + plot_df$`c13_Bgc_DZ-like_percent`
plot_df$GCB_sample = ifelse(plot_df$GCB_percent > 0, "GCB", "no-GCB")

i = "A_percent"

for (i in c("N_percent","M_percent","A_percent")) {
  tmp = plot_df[, c("GCB_sample", i)]
  colnames(tmp)[ncol(tmp)] = "abundance"
  
  ggplot(tmp,
         aes(x = GCB_sample,
             y = abundance,
             color = GCB_sample)) +
    geom_boxplot(outlier.color = NA,
                 lwd = 0.3) +
    geom_jitter(size = 0.6,shape = 16, stroke = 0,
                width = 0.2) +
    ggsignif::geom_signif(
      color = "black",
      comparisons = list(c("GCB", "no-GCB")),
      test = wilcox.test,
      step_increase = -0.1,
      textsize = 7 * 0.35,
      size = 0.2
    ) +
    scale_color_manual(values = c("GCB" = "#efbf05", "no-GCB" = "#0173c1")) +
    xlab("") + ylab("Percentage") +
    cowplot::theme_cowplot() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(size = 8, hjust = 0.5, face = 'plain'),
      title = element_text(size = 8),
      plot.margin = unit(c(0, 0, 0, 0), "char"),
      axis.line = element_line(size = 0.3),
      axis.ticks = element_line(size = 0.3),
      strip.background = element_rect(color = NA, fill = NA),
      strip.text = element_text(size = 8),
      legend.position = "none"
    )
  
  ggsave(paste0(fdir, "Figure3F_", i, ".pdf"),
         width = 1.3,
         height = 2.2)
}

#----- 02. BCR data -----#

BCR$SHM_rate = (1-BCR$v_identity) * 100
BCR$CellID = BCR$sequence_id
tmp = BCR[,c("CellID","clone_id","clonal","c_call","SHM_rate")]
BCR_df = tmp %>% left_join(meta)

##--- Figure3C; Proportions of cells sharing BCR with intratumoral c12_Bgc_LZ-like cells from each B cell subset. 
plot_df = BCR_df %>%
  filter(Tissue == "Tumor")

tmp = data.frame()

for (i in sort(unique(plot_df$Annotation))) {
  #
  tmp1 = plot_df %>% filter(Annotation %in% c("c12_Bgc_LZ-like"))
  tmp2 = plot_df %>% filter(clone_id %in% tmp1$clone_id) %>% filter(Annotation == i) %>%
    group_by(Annotation, SampleID) %>% summarise(overlapped_n = n())
  #
  tmp = rbind(tmp2, tmp)
}

sum_df = plot_df %>% group_by(Annotation, SampleID) %>% summarise(celltype_n = n())
sum_df = left_join(sum_df, tmp) %>% filter(Annotation != "c12_Bgc_LZ-like", celltype_n > 2) 
sum_df$percent = sum_df$overlapped_n / sum_df$celltype_n * 100
sum_df[is.na(sum_df)] = 0

ggplot(sum_df, aes(x = Annotation, y = percent, color = Annotation)) +
  geom_boxplot(outlier.color = NA, lwd = 0.2) +
  geom_jitter(size = 0.5,shape = 16,stroke = 0,width = 0.2) +
  scale_color_manual(name = "", values = Subset_color_panel) +
  xlab("") + ylab(paste0("% cells sharing BCR with","\n", "c12_Bgc_LZ-like")) +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    text = element_text(size = 8),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    legend.position = "None",
    plot.margin = unit(c(0,0,0,0), "char")
  ) 

ggsave(
  paste0(fdir,"Figure3C.pdf"),
  width = 3.2,
  height = 2.7
)

##--- Figure S4C; Heatmap showing the pTrans index for every combination of TIB subsets.
library("Startrac")
library("tictoc")
library("ggpubr")
library("ggplot2")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")

# Run STARTRAC
in.dat <-
  BCR_df[, c("CellID", "clone_id", "PatientID", "Annotation", "Tissue")]
in.dat <- in.dat %>% filter(Tissue == "Tumor")
colnames(in.dat) <-
  c("Cell_Name", 'clone.id', 'patient', 'majorCluster', 'loc')
out <- Startrac.run(in.dat, proj = "panB", verbose = F)

# plot
tmp = out@pIndex.tran %>% filter(aid == "panB") %>% as.data.frame()
rownames(tmp) = tmp$majorCluster
tmp = subset(tmp, select = -c(majorCluster, aid))
m = as.matrix(tmp)

summary(m)
col_fun <- colorRamp2(c(0, 0.05, 0.1, 0.2),
                      c("#fffde7", "#ffe0b2", "#ff9800", "#e65100"))

pdf(paste0(fdir, "FigureS4C.pdf"),
    width = 5.5,
    height = 5)
Heatmap(
  m,
  col = col_fun,
  cluster_rows = T,
  cluster_columns = T,
  column_names_gp = grid::gpar(fontsize = 6),
  row_names_gp = grid::gpar(fontsize = 6),
  width = ncol(m) * unit(0.12, "inch"),
  height = nrow(m) * unit(0.12, "inch"),
  name = "Number of overlaped clones",
)
dev.off()


##--- Heatmaps showing B cell major lineage marker expression among B cells from the lineage tree
require(ComplexHeatmap)
BCR$cdr3_aa <- as.character(Biostrings::translate(Biostrings::DNAStringSet(BCR$cdr3)))

i_clone_id = "48280"

tmp_BCR = filter(BCR, clone_id == i_clone_id) %>% as.data.frame()
tmp_obj = subset(obj, CellID %in% tmp_BCR$sequence_id)

matrix = tmp_obj@assays$RNA@data[c(
  c("BCL6", "AICDA", "RGS13", "NEIL1"),
  c("STMN1", "HMGB2"),
  c("MKI67", "BIRC5"),
  c("MZB1", "XBP1", "SDC1"),
  c("CD27", "TNFRSF13B", "GPR183")
),]
matrix = as.matrix(matrix)

rownames(tmp_BCR) = tmp_BCR$sequence_id
tmp_BCR = tmp_BCR[colnames(matrix),]
colnames(matrix) = tmp_BCR$cdr3_aa

tmp_obj$heatanno = ifelse(
  tmp_obj$Annotation %in% c("c12_Bgc_LZ-like", "c13_Bgc_DZ-like"),
  as.character(tmp_obj$Annotation),
  as.character(tmp_obj$Annotation_major_2)
)

ha = HeatmapAnnotation(Annotation = tmp_obj$heatanno,
                       col = list(
                         Annotation = c(
                           "c12_Bgc_LZ-like" = "#593202",
                           "c13_Bgc_DZ-like" = "#ECE4B7",
                           "Bm" = "#0072B5FF",
                           "ASC" = "#20854EFF"
                         )
                       ))

col_fun <- circlize::colorRamp2(c(0, 1, 2, 4),
                                c("#fffde7", "#ffe0b2", "#ff9800", "#e65100"))

pdf(paste0(fdir, "Figure3_", i_clone_id, "_heatmap.pdf"),
    width = 7,
    height = 5)

Heatmap(
  matrix,
  col = col_fun,
  cluster_rows = F,
  cluster_columns = T,
  column_names_gp = grid::gpar(fontsize = 6),
  row_names_gp = grid::gpar(fontsize = 6),
  width = ncol(matrix) * unit(0.1, "inch"),
  height = nrow(matrix) * unit(0.12, "inch"),
  name = "Exprssion",
  top_annotation = ha
)

dev.off()


#----- 03. CD45_obj -----#
CD45_obj = read_rds("obj/obj_CD45.rds")
CD45_meta = CD45_obj@meta.data

# Identify and label CD4 and CD8 T cells
CD45_meta <- CD45_meta %>%
  mutate(Annotation_CD45_major_2 = case_when(
    grepl("CD4", Annotation_CD45) ~ "CD4T",
    grepl("CD8", Annotation_CD45) ~ "CD8T",
    TRUE ~ Annotation_CD45_major
  ))

# Filter, summarize, and calculate percentages in one step
tmp <- CD45_meta %>%
  filter(Annotation_CD45_major_2 %in% c("Myeloid", "CD4T", "CD8T", "B")) %>%
  group_by(SampleID, Annotation_CD45_major_2, Annotation_CD45) %>%
  summarise(subtype_n = n(), .groups = 'drop') %>%
  group_by(SampleID) %>%
  mutate(Sample_n = sum(subtype_n)) %>%
  group_by(SampleID, Annotation_CD45_major_2) %>%
  mutate(major_n = sum(subtype_n)) %>%
  ungroup() %>%
  filter(Sample_n > 100) %>%
  mutate(
    sub_percent = round(subtype_n / Sample_n * 100, 2),
    major_percent = round(subtype_n / major_n * 100, 2)
  ) %>%
  select(Annotation_CD45, SampleID, major_percent) %>%
  rename(subtype = Annotation_CD45, freq = major_percent) %>%
  spread(key = subtype, value = freq, fill = 0)

# Summarize meta information and join with tmp
meta_brief <- CD45_meta %>%
  group_by(SampleID, Tissue, Cancer) %>%
  summarise()

tmp <- tmp %>%
  left_join(meta_brief, by = "SampleID")

# Filter for tumor tissue samples
df.T_freq <- tmp %>%
  filter(Tissue == "Tumor")

##--- Figure3G; pearson correlation between the abundance of Bgc and CD4 Tfh cells in tumors
plot_df = as.data.frame(df.T_freq)
plot_df$GCB = plot_df$`c12_Bgc_LZ-like`+plot_df$`c13_Bgc_DZ-like`

plot_df = plot_df[, c("Cancer",
                      "GCB",
                      "CD4.c16.Tfh.CXCR5")]

colnames(plot_df)[2:3] = c("x", "y")


ggplot(plot_df, aes(x = x, y = y)) +
  geom_point(aes(color = Cancer), size = 1, shape = 16, stroke = 0) +
  geom_smooth(
    method = 'lm',
    se = F,
    color = 'red',
    size = 0.3
  ) +
  ggpubr::stat_cor(method = "pearson",
                   output.type = "text",
                   size = 7 * 0.35) +
  scale_color_manual(values = Cancer_color_panel, name = "Cancer") +
  xlab("Percentage of GCB \n among TIBs") + ylab("Percentage of CD4 Tfh \n among CD4 T cells") +
  cowplot::theme_cowplot() +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    text = element_text(size = 8),
    plot.margin = unit(c(0, 0, 0, 0), "char"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    legend.position = "right"
  ) +
  xlim(-0.1, 60) +
  ylim(-0.1, 60) +
  guides(color = guide_legend(
    ncol = 3,
    override.aes = list(size = 2, alpha = 1)
  ))+
  theme(legend.text = element_text(size = 7),
        legend.spacing.y = unit(0, 'cm'),
        legend.key.height = unit(0.4, "cm"),
        legend.box.spacing = unit(0, 'cm'))

ggsave(paste0(fdir,"Figure3G.pdf"),width = 4, height = 2.2)

##--- FigureS4F; Heatmap showing the Pearson correlation between the frequencies of cell clusters in tumors
library(rstatix)
library(ggcorrplot)
library(ComplexHeatmap)
library(circlize)

m = df.T_freq[,2:82]
rownames(m) = df.T_freq$SampleID
corr = cor_mat(data = m)
p = cor_get_pval(corr)

corr_list = list(corr = corr,p = p)

corr_df = corr[,2:82] %>% as.data.frame()
rownames(corr_df) = corr$rowname

p_df = p[,2:82] %>% as.data.frame()
rownames(p_df) = p$rowname

tmp = corr_df[,c("c12_Bgc_LZ-like", "c13_Bgc_DZ-like")]
tmp = tmp[order(tmp$`c12_Bgc_LZ-like`,decreasing = T),]
tmp = tmp[1:10,]

summary(tmp)

col_fun <- colorRamp2(c(0, 0.5, 1),
                      c("#0F7B9F", "white", "#D83215"))

pdf(paste0(fdir,"FigureS4F.pdf"),width = 3,height = 3)
Heatmap(
  tmp,
  col = col_fun,
  cluster_rows = F,
  cluster_columns = F,
  column_names_gp = grid::gpar(fontsize = 6),
  row_names_gp = grid::gpar(fontsize = 6),
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", tmp[i, j]), x, y, gp = gpar(fontsize = 6))
  },
  width = ncol(tmp) * unit(0.4, "inch"),
  height = nrow(tmp) * unit(0.2, "inch"),
  name = paste0("Pearson correlation","\n","coefficient"),
)
dev.off()

p_df[rownames(tmp),colnames(tmp)] < 0.001

#----- 04. B obj DE -----# 
##--- Figure 3A; Top pathways enriched in intratumoral c12 and c13 Bgc cells compared with each other
library(Seurat)
library(readr)
library(dplyr)
library(ggplot2)

# calculate DE genes
tmp_obj = subset(obj,Annotation_major_2 == "Bgc" & Tissue == "Tumor")
Idents(tmp_obj) = tmp_obj$Annotation

all_markers = Seurat::FindMarkers(
  tmp_obj,
  ident.1 = "c12_Bgc_LZ-like",
  ident.2 = "c13_Bgc_DZ-like",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.25
)

all_markers$gene = rownames(all_markers)
write_rds(all_markers, "Bgc_DE.rds")

all_markers = read_rds("Bgc_DE.rds")

# perform GO
GO_a <- function(interest_gene, global_gene) {
  require(clusterProfiler)
  
  # Convert gene symbols to ENTREZ IDs for interest and global genes
  gene_entrez_id2GO <-
    clusterProfiler::bitr(interest_gene, "SYMBOL", "ENTREZID", "org.Hs.eg.db", drop = TRUE)$ENTREZID
  
  universe_gene_entrez_id2GO <-
    clusterProfiler::bitr(global_gene, "SYMBOL", "ENTREZID", "org.Hs.eg.db", drop = TRUE)$ENTREZID
  
  # Perform GO enrichment analysis and simplify results
  enr_res <- clusterProfiler::enrichGO(
    gene = gene_entrez_id2GO,
    universe = universe_gene_entrez_id2GO,
    ont = "BP",
    OrgDb = 'org.Hs.eg.db'
  )
  enr_res2 <- clusterProfiler::simplify(enr_res)
  
  # Generate plots
  p1 <- goplot(enr_res)
  p2 <- barplot(enr_res2, showCategory = 20, color = "p.adjust")
  
  # Return results and plots as a list
  list(
    go_res = enr_res,
    go_res_simplify = enr_res2,
    p1 = p1,
    p2 = p2
  )
}

GCB_c01 = all_markers %>% filter(p_val_adj < 0.05, avg_log2FC < -0.25)
GCB_c02 = all_markers %>% filter(p_val_adj < 0.05, avg_log2FC > 0.25)

GCB_C01_go = GO_a(rownames(GCB_c01),rownames(obj))
GCB_C02_go = GO_a(rownames(GCB_c02),rownames(obj))

# Process and combine GO results
process_GO <- function(GO_data, annotation) {
  GO_data$go_res_simplify@result %>%
    head(10) %>%
    mutate(
      annotation = annotation,
      log10q = -log10(abs(qvalue)))
}

GO_result <- bind_rows(
  process_GO(GCB_C02_go, "c12_Bgc_LZ-like"),
  process_GO(GCB_C01_go, "c13_Bgc_DZ-like")
)

GO_result$log10q = ifelse(GO_result$annotation == "c12_Bgc_LZ-like", GO_result$log10q * -1, GO_result$log10q)

# Plot the GO results
ggplot(GO_result) +
  geom_bar(aes(
    x = reorder(Description, log10q),
    y = log10q,
    fill = annotation
  ),
  stat = "identity",
  color = "white") +
  geom_text(aes(
    x = reorder(Description, log10q),
    y = 0,
    label = Description
  ),
  size = 7 * 0.35,
  angle = 0) +
  scale_fill_manual(
    name = "",
    values = c("c12_Bgc_LZ-like" = "#0f7b9f", "c13_Bgc_DZ-like" = "#d83215")
  ) +
  coord_flip() +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(size = 7),
    text = element_text(size = 7),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "", y = "") +
  ylim(-39, 60)

ggsave(paste0(fdir,"Figure3A.pdf"),width = 5,height = 4)  
