# Figure2.R

fdir = ""

## load packages
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

## load data
obj = read_rds("All_obj.rds")
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

## Figure 2D; TIB major lineage compositions across cancer types; only samples with B > 50 are shown; one-way ANOVA test
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

library(patchwork)
psum[[1]] + psum[[2]] + psum[[3]] + psum[[4]] + plot_layout(ncol = 1)
ggsave(paste0(fdir,"Figure2D.pdf"),
       width = 4,height = 7.5)
