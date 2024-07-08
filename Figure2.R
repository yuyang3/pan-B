# Figure2.R

fdir = ""

## load packages
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## load data
obj = read_rds("obj/All_obj.rds")
meta = obj@meta.data

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
##--- Figure 2A; Proportions of TIBs in CD45+ cells across cancer types. One-way ANOVA
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

df.sample_count = CD45_meta %>% 
  group_by(SampleID,PatientID,Cancer,Tissue_short,Sample_type,Treatment_status) %>% 
  summarise(CD45_n=n())

for(i in sort(unique(CD45_meta$majortype))){
  
  tmp = CD45_meta %>% 
    filter(majortype == i) %>% 
    group_by(SampleID) %>% 
    summarise(n = n())
  
  tmp[is.na(tmp)] = 0
  
  df.sample_count = left_join(df.sample_count,tmp)
  df.sample_count[is.na(df.sample_count)] = 0
  
  df.sample_count$percent = 
    round(df.sample_count$n/df.sample_count$CD45_n * 100,2)
  
  i_names = paste0(i,"_n")    
  colnames(df.sample_count)[length(df.sample_count)-1] = i_names        
  
  i_names = paste0(i,"_percent")
  colnames(df.sample_count)[length(df.sample_count)] = i_names    
  
}

df.sample_count_flt = df.sample_count %>% filter(CD45_n > 100, Treatment_status == "treatment naïve")

# Define parameters
i <- "B"
tissue <- "T"
i_names <- paste0(i, "_percent")

# Select relevant columns and rename the last column to "percent"
tmp <- df.sample_count_flt[, c("SampleID", "Cancer", "Tissue_short", i_names)]
colnames(tmp)[length(tmp)] <- "percent"

# Filter data for specific tissue type and calculate order value based on median percent
tmp2 <- tmp %>% 
  filter(Tissue_short %in% tissue) %>%
  group_by(Cancer) %>% 
  mutate(order_value = median(percent))

# Perform ANOVA and extract p-value and F-value
aov_results <- aov(percent ~ Cancer, data = tmp2)
aov_summary <- summary(aov_results)
p_value <- aov_summary[[1]][["Pr(>F)"]][1]
F_value <- aov_summary[[1]][["F value"]][1]

# Create boxplot with jitter and customize appearance
ggplot(tmp2, aes(
  x = reorder(Cancer, order_value),
  y = percent,
  color = Cancer
)) +
  geom_boxplot(outlier.color = NA, lwd = 0.3) +
  geom_jitter(size = 0.2, width = 0.2) +
  scale_color_manual(values = Cancer_color_panel) +
  labs(x = "", y = expression(paste("% among CD45+ cells \nin tumor tissues")),
       title = paste0("P = ", signif(p_value, 2), ", F = ", round(F_value, 2)))+
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(size = 7,angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    text = element_text(size = 8),
    plot.title = element_text(size = 7, hjust = 0, face = 'plain'),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "char")
  ) 

# Save plot as PDF
ggsave(
  filename = paste0(fdir, "Figure2A.pdf"),
  device = "pdf",
  width = 3.7,
  height = 2.2
)

##--- Figure S3A; Proportions of B cells in CD45+ cells from ANTs across cancer types. One-way ANOVA test.
# Define parameters
i <- "B"
tissue <- "N"
i_names <- paste0(i, "_percent")

# Select relevant columns and rename the last column to "percent"
tmp <- df.sample_count_flt[, c("SampleID", "Cancer", "Tissue_short", i_names)]
colnames(tmp)[length(tmp)] <- "percent"

# Filter data for specific tissue type and calculate order value based on median percent
tmp2 <- tmp %>% 
  filter(Tissue_short %in% tissue) %>%
  group_by(Cancer) %>% 
  mutate(order_value = median(percent))

# Perform ANOVA and extract p-value and F-value
aov_results <- aov(percent ~ Cancer, data = tmp2)
aov_summary <- summary(aov_results)
p_value <- aov_summary[[1]][["Pr(>F)"]][1]
F_value <- aov_summary[[1]][["F value"]][1]

# Create boxplot with jitter and customize appearance
ggplot(tmp2, aes(
  x = reorder(Cancer, order_value),
  y = percent,
  color = Cancer
)) +
  geom_boxplot(outlier.color = NA, lwd = 0.3) +
  geom_jitter(size = 0.2, width = 0.2) +
  scale_color_manual(values = Cancer_color_panel) +
  labs(x = "", y = expression(paste("% among CD45+ cells \nin adjacent non-tumor tissues")),
       title = paste0("P = ", signif(p_value, 2), ", F = ", round(F_value, 2)))+
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(size = 7,angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    text = element_text(size = 8),
    plot.title = element_text(size = 7, hjust = 0, face = 'plain'),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "char")
  ) 

ggsave(
  filename = paste0(fdir,"FigureS3A.pdf"),
  device = "pdf",
  width = 3,
  height = 2.2
)


##--- Figure 2B; Pearson correlation of B cell proportions between tumors and ANTs in CRC, STAD and ESCA.
tmp = df.sample_count_flt %>% filter(Tissue_short %in% c("T","N")) %>%
  group_by(PatientID) %>% mutate(tissue_n = length(unique(Tissue_short))) 
tmp = filter(tmp,tissue_n==2)

tmp = subset(tmp,select = c(PatientID,Tissue_short,B_percent))
tmp = tmp %>% spread(key = Tissue_short,value = B_percent)

brief_meta = df.sample_count %>% group_by(PatientID, Cancer) %>% summarise()
tmp = tmp %>% left_join(brief_meta)

plot_df = tmp  %>% filter(Cancer %in% c("CRC","ESCA","STAD"))

ggplot(plot_df, aes(x = N, y = T)) +
  geom_point(aes(color = Cancer),size = 0.2) +
  geom_smooth(method = 'lm', se = F, color = 'red',size = 0.3) +
  ggpubr::stat_cor(method = "pearson",output.type = "text",size = 3) +
  scale_color_manual(values = Cancer_color_panel) +
  xlab("% among CD45+ cells in adjacent non-tumor tissues") +
  ylab("% among CD45+ cells in tumor tissues") +
  facet_wrap( ~ Cancer, scale = 'free') +
  cowplot::theme_cowplot() +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    text = element_text(size = 8),
    title = element_text(size = 8),
    plot.margin = unit(c(0, 0, 0, 0), "char"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    strip.background = element_rect(color = NA,fill = NA),
    strip.text = element_text(size = 8),
    legend.position = "None"
  )
ggsave(paste0(fdir,"Figure2B.pdf"),width = 4.3,height = 2)

##--- Figure S3B; Boxplot comparing the proportion of B cells in CD45+ cells between tumors 
plot_df = df.sample_count_flt %>% 
  filter(Tissue_short %in% c("T", "N"))

tmp = plot_df[, c("SampleID", "Cancer", "Tissue_short", "B_percent")]
colnames(tmp)[length(tmp)] = "percent"
tmp = filter(tmp, Cancer %in% c("NSCLC","BRCA","PACA"))

ggplot(tmp, aes(x = Tissue_short, y = percent, color = Tissue_short)) +
  geom_boxplot(outlier.color = NA,
               lwd = 0.3) +
  geom_jitter(size = 0.2,
              width = 0.2) +
  ggsignif::geom_signif(
    color = "black",
    comparisons = list(c("T", "N")),
    test = wilcox.test,
    textsize = 7*0.35
  ) +
  scale_color_manual(values = c("T" = "#5BC0EB","N" = "#9BC53D")) +
  xlab("") + ylab(expression(paste("% among CD45+ cells"))) +
  facet_wrap( ~ Cancer, scale = 'free',nrow = 1) +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(size = 8, hjust = 0.5, face = 'plain'),
    title = element_text(size = 8),
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    strip.background = element_rect(color = NA,fill = NA),
    strip.text = element_text(size = 8),
    legend.position = "Right"
  )

ggsave(
  filename = paste0(dir, "FigureS3B.pdf"),
  device = "pdf",
  width = 3.5,
  height = 2.2
)

##--- Figure 2C; Upset plot showing the infiltraion status of three major immune components within tumors.


##--- Figure 2D; TIB major lineage compositions across cancer types; only samples with B > 50 are shown; one-way ANOVA test
### major type abundance data frame
df.sample_count <- meta %>%
  group_by(SampleID, PatientID, Cancer, Tissue, Treatment_status) %>%
  summarise(B_n = n(), .groups = 'drop')

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

df.sample_count_flt = df.sample_count %>% filter(B_n > 50, Treatment_status == "treatment naïve")

### plot code
psum = list()

for (i in sort(unique(meta$Annotation_major_2))) {
  i_names = paste0(i, "_percent")
  
  tmp = df.sample_count_flt[, c("SampleID", "Cancer", "Tissue", eval(i_names))]
  colnames(tmp)[length(tmp)] = "percent"
  
  tissue = "Tumor"
  
  tmp2 = tmp %>% filter(Tissue %in% tissue)
  
  a = summary(aov(percent ~ Cancer, data = tmp2))
  p_value = a[[1]][["Pr(>F)"]][1]
  F_value = a[[1]][["F value"]][1]
  
  p = ggplot(tmp2, aes(
    x = Cancer,
    y = percent,
    color = Cancer
  )) +
    geom_boxplot(outlier.color = NA,
                 lwd = 0.3) +
    geom_jitter(size = 0.2,
                width = 0.2) +
    scale_color_manual(values = Cancer_color_panel) +
    xlab("") + ylab(expression(paste("% among TIBs"))) +
    ggtitle(paste0(i),
            paste0("P = ", signif(p_value, 2), ", F = ", round(F_value, 2))) +
    cowplot::theme_cowplot() +
    theme(
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(
        size = 10,
        hjust = 0.5,
        face = 'plain'
      ),
      title = element_text(size = 8),
      plot.subtitle = element_text(size = 7),
      axis.line = element_line(size = 0.3),
      axis.ticks = element_line(size = 0.3),
      legend.position = "None",
      plot.margin = unit(c(0, 0, 0, 0), "char")
    )
  
  psum[[i]] = p
}

psum[[1]] + psum[[2]] + psum[[3]] + psum[[4]] + plot_layout(ncol = 1)
ggsave(paste0(fdir,"Figure2D.pdf"),
       width = 4,height = 7.5)

##--- Figure S2F
psum = list()

for (i in sort(unique(meta$Annotation_major_2))) {
  i_names = paste0(i, "_percent")
  
  tmp = df.sample_count_flt[, c("SampleID", "Cancer", "Tissue", eval(i_names))]
  colnames(tmp)[length(tmp)] = "percent"
  
  tissue = "Adjacent non-tumor tissue"
  
  tmp2 = tmp %>% filter(Tissue %in% tissue)
  
  a = summary(aov(percent ~ Cancer, data = tmp2))
  p_value = a[[1]][["Pr(>F)"]][1]
  F_value = a[[1]][["F value"]][1]
  
  p = ggplot(tmp2, aes(
    x = Cancer,
    y = percent,
    color = Cancer
  )) +
    geom_boxplot(outlier.color = NA,
                 lwd = 0.3) +
    geom_jitter(size = 0.2,
                width = 0.2) +
    scale_color_manual(values = Cancer_color_panel) +
    xlab("") + ylab(expression(paste("% among TIBs"))) +
    ggtitle(paste0(i),
            paste0("P = ", signif(p_value, 2), ", F = ", round(F_value, 2))) +
    cowplot::theme_cowplot() +
    theme(
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(
        size = 10,
        hjust = 0.5,
        face = 'plain'
      ),
      title = element_text(size = 8),
      plot.subtitle = element_text(size = 7),
      axis.line = element_line(size = 0.3),
      axis.ticks = element_line(size = 0.3),
      legend.position = "None",
      plot.margin = unit(c(0, 0, 0, 0), "char")
    )
  
  psum[[i]] = p
}

psum[[1]] + psum[[2]] + psum[[3]] + psum[[4]] + plot_layout(ncol = 2)
ggsave(paste0(fdir,"FigureS3F.pdf"),
       width = 6,height = 4)

##--- Figure 2E; Cancer type stratification by hierarchical clustering of B cell subset proportions, with sample number annotated
tmp = df.sample_count_flt %>%
  filter(Tissue == "Tumor")

cluster_df = meta %>% filter(SampleID %in% tmp$SampleID) %>%
  group_by(Cancer) %>%
  mutate(sample_n = length(unique(SampleID))) %>%
  mutate(Cancer = paste0(Cancer, " (N = ", sample_n, ")")) %>% 
  # calculate percentage for each cancer
  group_by(Cancer, Annotation) %>%
  summarise(n = n()) %>% group_by(Cancer) %>% mutate(Cancer_n = sum(n)) %>% ungroup() %>%
  mutate(percentage = n / Cancer_n) %>%
  # spread data frame into data frame for further analysis
  subset(select = c(Annotation, percentage, Cancer)) %>%
  spread(key = Cancer, value = percentage, fill = 0) %>%
  as.data.frame()

rownames(cluster_df) = cluster_df$Annotation
cluster_df = subset(cluster_df, select = -c(Annotation))
seeds_df = t(cluster_df)

require(ggtree)

dist_mat <- dist(seeds_df, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
ggtree(hclust_avg) + geom_tiplab()

ggsave(
  paste0(fdir, "Figure2E_sub_tree.pdf"),
  height = 4,
  width = 2
)

plot_df = meta %>% filter(SampleID %in% tmp$SampleID) %>%
  group_by(Cancer) %>%
  mutate(sample_n = length(unique(SampleID))) %>%
  mutate(Cancer = paste0(Cancer, " (N = ", sample_n, ")")) %>% 
  group_by(Cancer, Annotation) %>%
  dplyr::summarise(count = n())

levels = hclust_avg[["labels"]][hclust_avg[["order"]]]
plot_df$Cancer <- factor(plot_df$Cancer, levels = levels)

ggplot(plot_df, aes(x = Cancer, y = count, fill = Annotation)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = Subset_color_panel,
                    name = "Subsets") +
  theme_cowplot() +
  xlab("") + ylab("Proportion") +
  theme(
    axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 7),
    text = element_text(size = 8),
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    legend.position = "bottom"
  ) + 
  guides(fill = guide_legend(ncol = 3))

ggsave(
  file.path(fdir,"Figure2E_sub_bar.pdf"),
  width = 4.2,
  height = 4.5
)
