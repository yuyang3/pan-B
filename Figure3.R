fdir = ""

## load packages
library(readr)
library(dplyr)
library(ggplot2)

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

## Figure 3B; Pearson correlation between the abundance of c12 and c13 Bgc cells in tumor samples with B cells > 50.
### sub type abundance data frame
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

### plot code
plot_df = as.data.frame(df.sample_count_flt) %>% filter(Tissue == "Tumor")

plot_df = plot_df[, c("Cancer",
                      "c12_Bgc_LZ-like_percent",
                      "c13_Bgc_DZ-like_percent")]

colnames(plot_df)[2:3] = c("x", "y")
plot_df$Cancer = as.character(plot_df$Cancer)
max(plot_df$y)

ggplot(plot_df, aes(x = x, y = y)) +
  geom_point(aes(color = Cancer), size = 0.2) +
  geom_smooth(
    method = 'lm',
    se = F,
    color = 'red',
    size = 0.3
  ) +
  ggpubr::stat_cor(method = "pearson",
                   output.type = "text",
                   size = 7 * 0.35) +
  scale_color_manual(values = Cancer_color_panel) +
  xlab("c12_Bgc_LZ-like") + ylab("c13_Bgc_DZ-like") +
  cowplot::theme_cowplot() +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    title = element_text(size = 8),
    plot.margin = unit(c(0, 0, 0, 0), "char"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    strip.background = element_rect(color = NA, fill = NA),
    strip.text = element_text(size = 8),
    legend.position = "none"
  ) +
  xlim(-0.1, 42) +
  ylim(-0.1, 42)

ggsave(paste0(fdir, "Figure3B.pdf"),
       width = 2,
       height = 2)
