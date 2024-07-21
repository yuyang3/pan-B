fdir = ""
## I. load packages
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

## II. load data
obj = read_rds("obj/All_obj.rds")
meta = obj@meta.data
BCR = read_rds("BCR.rds")
CD45_obj = read_rds("obj/obj_CD45.rds")

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

##--- Figure 4A; Tissue preference of each B cell subset evaluated by the Ro/e index.

#--- 01. major type abundance data frame ---#
##--- Figure S5A; B cell major lineage compositions in the blood, ANTs and tumors
if(TRUE){
  df.sample_count = meta %>% 
    group_by(SampleID,PatientID,Cancer,Tissue,Treatment_status) %>%
    summarise(B_n=n())
  
  for(i in sort(unique(meta$Annotation_major_2))){
    tmp = meta %>% 
      filter(Annotation_major_2 == i) %>% 
      group_by(SampleID) %>% 
      summarise(n = n())
    
    tmp[is.na(tmp)] = 0
    
    df.sample_count = left_join(df.sample_count,tmp,by = "SampleID")
    df.sample_count[is.na(df.sample_count)] = 0
    
    df.sample_count$percent = 
      round(df.sample_count$n/df.sample_count$B_n * 100,2)
    
    i_names = paste0(i,"_n")    
    colnames(df.sample_count)[length(df.sample_count)-1] = i_names        
    
    i_names = paste0(i,"_percent")
    colnames(df.sample_count)[length(df.sample_count)] = i_names    
    
  }
  
  df.sample_count_flt = df.sample_count %>% filter(B_n > 50, Treatment_status == "treatment naïve")
}

# 
plot_df = df.sample_count_flt %>% filter(Tissue %in% c("Blood","Adjacent non-tumor","Tumor"))
plot_df = plot_df[, c(1, seq(8, ncol(plot_df), 2), 4)]

#
plot_df = plot_df %>% gather(key = "Annotation", value = "percent", -SampleID, -Tissue)
plot_df$Annotation = stringr::str_replace(plot_df$Annotation, "_percent", "")

#
plot_df$Tissue = factor(plot_df$Tissue, levels = c("Blood","Adjacent non-tumor","Tumor"))
plot_df$Annotation =  factor(plot_df$Annotation, levels = c("Bn",  "Bm", "Bgc", "ASC"))

ggplot(plot_df, aes(x = Annotation, y = percent, color = Annotation)) +
  geom_boxplot(outlier.color = NA,
               lwd = 0.3) +
  geom_jitter(
    size = 0.4,
    shape = 16,
    stroke = 0,
    width = 0.2
  ) +
  scale_color_manual(name = "",values = Major_color_panel) +
  xlab("") + ylab("% among B cells") +
  cowplot::theme_cowplot() +
  facet_grid(Tissue ~ ., switch = "y") +
  theme(
    axis.text.x = element_text(size = 6,angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 7),
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    strip.text = element_text(size = 5)
  ) +
  guides(fill = FALSE, color = FALSE)

ggsave(paste0(fdir, "FigureS4A.pdf"),height = 3,width = 1.5)

#--- 02. subset abundance data frame ---#
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

##--- Figure S5B; B cell subsets compositions in the blood, ANTs and tumors.
plot_df = df.sample_count_flt %>% filter(Tissue %in% c("Blood","Adjacent non-tumor","Tumor"))
plot_df = plot_df[, c(1, seq(8, ncol(plot_df), 2), 4)]

#
plot_df = plot_df %>% gather(key = "Annotation", value = "percent", -SampleID, -Tissue)
plot_df$Annotation = stringr::str_replace(plot_df$Annotation, "_percent", "")

#
plot_df$Tissue = factor(plot_df$Tissue, levels = c("Blood","Adjacent non-tumor","Tumor"))

ggplot(plot_df, aes(x = Annotation, y = percent, color = Annotation)) +
  geom_boxplot(outlier.color = NA,
               lwd = 0.3) +
  geom_jitter(
    size = 0.4,
    shape = 16,
    stroke = 0,
    width = 0.2
  ) +
  scale_color_manual(name = "",values = Subset_color_panel) +
  xlab("") + ylab("% among B cells") +
  cowplot::theme_cowplot() +
  facet_grid(Tissue ~ ., switch = "y") +
  theme(
    axis.text.x = element_text(size = 6,angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 7),
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    strip.text = element_text(size = 5)
  ) +
  guides(fill = FALSE, color = FALSE)

ggsave(paste0(fdir, "FigureS4B.pdf"),height = 3.6,width = 4)

##--- Figure S5H; Boxplots comparing the proportions of c16 and c17 PCs in total B cells between tumors and ANTs. 
i = "c16_PC_IGHG"
i = "c17_PC_IGHA"

if(TRUE){
  i_names = paste0(i, "_percent")
  
  tmp = df.sample_count_flt %>% filter(Tissue %in% c("Adjacent non-tumor","Tumor"))
  
  tmp = tmp[, c("SampleID", "Cancer", "Tissue", eval(i_names))]
  colnames(tmp)[length(tmp)] = "percent"
  
  tmp$Tissue = factor(tmp$Tissue, levels = c("Adjacent non-tumor","Tumor"))
  
  
  ggplot(tmp, aes(x = Tissue,
                  y = percent,
                  color = Tissue)) +
    geom_boxplot(outlier.color = NA,
                 lwd = 0.3) +
    geom_jitter(size = 0.3,
                shape = 16,
                stroke = 0,
                width = 0.2) +
    ggsignif::geom_signif(
      color = "black",
      comparisons = list(c("Tumor", "Adjacent non-tumor")),
      test = wilcox.test,
      step_increase = -1,
      textsize = 5*0.35
    ) +
    scale_color_manual(values = Tissue_color_panel) +
    xlab("") + ylab("% among B cells") +
    ggtitle(i)+
    cowplot::theme_cowplot() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 6),
      plot.title = element_text(size = 7, hjust = 0.5, face = 'plain'),
      text = element_text(size = 7),
      plot.margin = unit(c(0,0,0,0), "char"),
      axis.line = element_line(size = 0.3),
      axis.ticks = element_line(size = 0.3),
      strip.background = element_rect(color = NA,fill = NA),
      strip.text = element_text(size = 8),
      legend.position = "none"
    )
  
  ggsave(paste0(fdir,"FigureS5H_",i,"_NT.pdf"),width = 0.9,height = 1.5)
  
}

##--- Figure 4J; Boxplots comparing the abundance ration of c16 to c17 PCs between tumors and ANTs
tmp = df.sample_count_flt %>% filter(c16_PC_IGHG_n > 0,c17_PC_IGHA_n > 0)
tmp$Ig = log2(tmp$c16_PC_IGHG_n/tmp$c17_PC_IGHA_n)

tmp2 = tmp %>% filter(Tissue %in% c("Tumor", "Adjacent non-tumor")) %>% 
  group_by(Cancer) %>% dplyr::summarise(n = length(unique(Tissue))) %>% filter(n == 2)

tmp = tmp[,c("Cancer","Ig","Tissue")] %>% 
  filter(Tissue %in% c("Tumor", "Adjacent non-tumor"),
         Cancer %in% tmp2$Cancer)

tmp$Tissue = factor(tmp$Tissue, levels = c("Adjacent non-tumor","Tumor"))

ggplot(filter(tmp, Cancer %in% c("CRC", "ESCA", "NSCLC", "STAD")), 
              aes(x = Tissue, y = Ig)) +
  geom_boxplot(aes(color = Tissue), outlier.colour = NA, lwd = 0.3) +
  geom_jitter(aes(color = Tissue), size = 0.4,shape = 16, stroke = 0, width = 0.2) +
  ggsignif::geom_signif(
    comparisons = list(c("Tumor", "Adjacent non-tumor")),
    test = wilcox.test,
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
    step_increase = -0.1,
    textsize = 2.5,
    y_position = 5.75,
    size = 0.2
  )+
  scale_color_manual(values = Tissue_color_panel) +
  facet_wrap(~ Cancer, nrow = 1) +
  xlab("") +
  ylab("log2(c16_PC_IGHG/\nc17_PC_IGHA)") +
  theme_bw() +  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 7),
    panel.grid.major = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    plot.title = element_text(size = 8),
    legend.position = "none"
  ) 

ggsave(paste0(fdir,"Figure4J.pdf"),width = 2.4, height = 1.6)

##--- Figure S5C; Boxplot comparing the compositional diversity of B cells in the blood, ANTs and tumors. 
df =  meta %>% 
  # filter 
  filter(Tissue %in% c("Blood","Adjacent non-tumor","Tumor"),
         Treatment_status == "treatment naïve") %>%
  # group by
  group_by(Tissue,Annotation,SampleID) %>% summarise(cell_count = n()) %>% 
  group_by(SampleID) %>% mutate(cell_count_sum = sum(cell_count)) %>% 
  ungroup() %>% mutate(percent = cell_count/cell_count_sum) %>%
  # filter
  filter(cell_count_sum > 50) %>% 
  # calculate H 
  group_by(SampleID,Tissue) %>% summarise(H=philentropy::H(percent)) %>%
  group_by(Tissue) %>% mutate(diversity = H / max(H)) %>%
  group_by(Tissue) %>% mutate(order_value = median(diversity))

ggplot(df, aes(
  x = reorder(Tissue,-order_value),
  y = diversity,
  color = Tissue
)) +
  geom_boxplot(outlier.color = NA,
               lwd = 0.3) +
  ggforce::geom_sina(size = 0.4,shape = 16, stroke = 0) +
  ggsignif::geom_signif(
    color = "black",
    comparisons = list(c("Tumor", "Adjacent non-tumor"),
                       c("Tumor", "Blood"),
                       c("Adjacent non-tumor", "Blood")),
    test = wilcox.test,
    step_increase = 0.2,
    textsize = 5*0.35
  ) +
  scale_color_manual(values = Tissue_color_panel) +
  xlab("") + ylab("Diversity") +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 7),
    plot.margin = unit(c(0,0,0,0), "char"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    legend.position = "right"
  ) 

ggsave(
  paste0(fdir,"FigureS5C.pdf"),
  width = 2.1,
  height = 1.3
)


#----- 03. Ro/e -----#
if(TRUE){
  divMatrix <- function(m1, m2){
    ## Divide each element in turn in two same dimension matrixes
    ##
    ## Args:
    #' @m1: the first matrix
    #' @m2: the second matrix
    ##
    ## Returns:
    ## a matrix with the same dimension, row names and column names as m1. 
    ## result[i,j] = m1[i,j] / m2[i,j]
    dim_m1 <- dim(m1)
    dim_m2 <- dim(m2)
    if( sum(dim_m1 == dim_m2) == 2 ){
      div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
      row.names(div.result) <- row.names(m1)
      colnames(div.result) <- colnames(m1)
      for(i in 1:dim_m1[1]){
        for(j in 1:dim_m1[2]){
          div.result[i,j] <- m1[i,j] / m2[i,j]
        }
      }   
      return(div.result)
    }
    else{
      warning("The dimensions of m1 and m2 are different")
    }
  }
  
  ROIE <- function(crosstab){
    ## Calculate the Ro/e value from the given crosstab
    ##
    ## Args:
    #' @crosstab: the contingency table of given distribution
    ##
    ## Return:
    ## The Ro/e matrix 
    rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
    rowsum.matrix[,1] <- rowSums(crosstab)
    colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
    colsum.matrix[1,] <- colSums(crosstab)
    allsum <- sum(crosstab)
    roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
    row.names(roie) <- row.names(crosstab)
    colnames(roie) <- colnames(crosstab)
    return(roie)
  }
}

##--- Figure 4A; Tissue preference of each B cell subset evaluated by the Ro/e index
plot_df = meta %>% 
  filter(Tissue %in% c("Blood","Adjacent non-tumor","Tumor"))

plot_df$Annotation = as.character(plot_df$Annotation)
plot_df$Tissue = as.character(plot_df$Tissue)

summary <- table(plot_df[,c('Annotation','Tissue')])
roe <- as.data.frame(ROIE(summary))
roe <- roe[,c("Blood","Adjacent non-tumor","Tumor")]

require(ComplexHeatmap)
require(circlize)
col_fun <- colorRamp2(c(0, 1, 1.5, 2),
                      c("#fffde7", "#ffe0b2","#ff9800", "#e65100"))

pdf(paste0(fdir,"Figure4A.pdf"),width = 3,height = 3.5)
Heatmap(
  roe,
  col = col_fun,
  cluster_rows = T,
  cluster_columns = F,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D",
  column_names_gp = grid::gpar(fontsize = 6),
  row_names_gp = grid::gpar(fontsize = 6),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", roe[i, j]), x, y, gp = gpar(fontsize = 6))
  },
  width = ncol(roe) * unit(0.3, "inch"),
  height = nrow(roe) * unit(0.15, "inch"),
  name = "Ro/e"
)
dev.off()

##--- Figure S5I; Tissue preference of c16 and c17 PCs across cancer types evaluated by the Ro/e index.
result_df <- NULL

for (i_Cancer in sort(unique(meta$Cancer))) {
  
  meta_subset <- meta %>%
    filter(Cancer == i_Cancer) %>%
    filter(Tissue_short %in% c("N", "T")) %>% 
    mutate(sample_n = length(unique(SampleID))) %>%
    mutate(Cancer = paste0(Cancer, " (N = ", sample_n, ")"))
  
  
  
  if (sum(meta_subset$Tissue_short == "N") == 0 |
      sum(meta_subset$Tissue_short != "N") == 0) {
    next
  }
  
  meta_subset$Annotation = as.character(meta_subset$Annotation)
  meta_subset$Tissue_short = as.character(meta_subset$Tissue_short)
  
  summary <- table(meta_subset[, c('Annotation', 'Tissue_short')])
  roe <- as.data.frame(ROIE(summary))
  i_Cancer <- unique(meta_subset$Cancer)
  
  result_df[[i_Cancer]] <- data.frame(
    Cancer = i_Cancer,
    Annotation = rownames(roe),
    roe = roe[, "T"],
    stringsAsFactors = FALSE
  )
}

result_df <- bind_rows(result_df)

tmp = spread(result_df, Annotation, roe) %>% as.data.frame()
rownames(tmp) = tmp$Cancer
tmp = subset(tmp, select = -c(Cancer))
tmp = as.matrix(tmp) %>% t()
tmp = tmp[c("c16_PC_IGHG","c17_PC_IGHA"),]
tmp = t(tmp)

col_fun <- circlize::colorRamp2(c(0, 0.5, 1, 1.5, 2),
                                c("#0081a7", "#00afb9", "#fdfcdc", "#fed9b7", "#f07167"))
pdf(
  file.path(fdir,"FigureS5I.pdf"),
  width = 4,
  height = 4
)

ComplexHeatmap::Heatmap(
  tmp,
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_side = "right",
  column_names_side = "bottom",
  # column_names_rot = 45,
  name = "Ro/e",
  heatmap_legend_param = list(
    legend_direction = "horizontal",
    title_position = "topcenter",
    legend_width = unit(0.8, "inch"),
    title_gp = gpar(fontsize = 7),
    labels_gp = gpar(fontsize = 6)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", tmp[i, j]), x, y, gp = gpar(fontsize = 6))
  },
  width = ncol(tmp) * unit(0.2, "inch"),
  height = nrow(tmp) * unit(0.2, "inch"),
  column_names_gp = grid::gpar(fontsize = 7),
  row_names_gp = grid::gpar(fontsize = 7),
)
dev.off()


#----- 04. BCR associated analysis -----#
BCR$SHM_rate = (1-BCR$v_identity) * 100
BCR$CellID = BCR$sequence_id
tmp = BCR[,c("CellID","clone_id","clonal","c_call","SHM_rate")]
BCR_df = tmp %>% left_join(meta)

##--- Figure S5D; Clonal expansion levels of B cells in the blood, ANTs and tumors, with cells categorized by the clone size of their corresponding clones.
plot_df = BCR_df %>%
  filter(Tissue %in% c("Blood","Adjacent non-tumor","Tumor")) %>%
  group_by(clone_id) %>% mutate(clone_size = n())

plot_df$clone_size_factor = ifelse(plot_df$clone_size > 3, ">3", plot_df$clone_size)

plot_df = plot_df %>% 
  group_by(Tissue, clone_size_factor) %>% summarise(cell_count = n()) %>%
  group_by(Tissue) %>% mutate(cell_count_sum = sum(cell_count)) %>%
  ungroup() %>% mutate(percent = cell_count / cell_count_sum) %>%
  filter(clone_size_factor  != "1")

plot_df$clone_size_factor = factor(plot_df$clone_size_factor, levels = c(">3", "3", "2","1"))
plot_df$Tissue = factor(plot_df$Tissue, levels = c("Blood", "Adjacent non-tumor", "Tumor"))

ggplot(plot_df,
       aes(x = Tissue, y = percent, fill = clone_size_factor)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c(
      ">3" = "#e76f51",
      "3" = "#f4a261",
      "2" = "#e9c46a",
      "1" = "#2a9d8f"
    ),
    name = 'Clone size'
  ) +
  xlab("") + ylab("Proportion") +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 7),
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    legend.position = "right"
  )

ggsave(paste0(fdir,"FigureS5D.pdf"),width = 1.9,height = 2.2)

##--- Figure 4B; clonal expansion levels of B cell major lineages in tumors, with cells categorized by the clone size of their corresponding clones.
plot_df = BCR_df %>%
  group_by(clone_id) %>% mutate(clone_size = n()) %>% 
  filter(Tissue == "Tumor")

plot_df$clone_size_factor = ifelse(plot_df$clone_size > 3, ">3", plot_df$clone_size)

plot_df = plot_df %>% group_by(Annotation_major_2, clone_size_factor) %>% summarise(cell_count = n()) %>%
  group_by(Annotation_major_2) %>% mutate(cell_count_sum = sum(cell_count)) %>%
  ungroup() %>% mutate(percent = cell_count / cell_count_sum) %>%
  filter(clone_size_factor  != "1")

plot_df$clone_size_factor = factor(plot_df$clone_size_factor, levels = c(">3", "3", "2","1"))
plot_df$Annotation_major_2 =  factor(plot_df$Annotation_major_2, levels = c("Bn", "Bm", "Bgc", "ASC"))

ggplot(plot_df,
       aes(x = Annotation_major_2, y = percent, fill = clone_size_factor)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab("") + ylab("Proportion") +
  scale_fill_manual(
    values = c(
      ">3" = "#e76f51",
      "3" = "#f4a261",
      "2" = "#e9c46a",
      "1" = "#2a9d8f"
    ),
    name = 'Clone size'
  ) +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_text(size = 6,angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 7),
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    legend.position = "none"
  )

ggsave(paste0(fdir,"Figure4B.pdf"),width = 1.4,height = 2)

##--- Figure S5E; Clonal expansion levels of TIB subsets, with cells categorized by the clone size of their corresponding clones.
plot_df = BCR_df %>%
  group_by(clone_id) %>%
  mutate(clone_size = n()) %>%
  filter(Tissue == "Tumor")

plot_df$clone_size_factor = ifelse(plot_df$clone_size > 3, ">3", plot_df$clone_size)

plot_df = plot_df %>% group_by(Annotation, Annotation_major_2, clone_size_factor) %>% summarise(cell_count = n()) %>%
  group_by(Annotation, Annotation_major_2) %>% mutate(cell_count_sum = sum(cell_count)) %>%
  ungroup() %>% mutate(percent = cell_count / cell_count_sum) %>%
  filter(clone_size_factor  != "1")

plot_df$clone_size_factor = factor(plot_df$clone_size_factor, levels = c(">3", "3", "2", "1"))
plot_df$Annotation_major_2 = factor(plot_df$Annotation_major_2, levels = rev(c("Bn","Bm","Bgc","ASC")))
plot_df$Annotation = stringr::str_extract(plot_df$Annotation, "c[0-9]{2}")

ggplot(plot_df, aes(x = percent, y = Annotation, fill = clone_size_factor)) +
  geom_bar(stat = "identity") +
  cowplot::theme_cowplot() +
  xlab("Proportion") + ylab("") +
  scale_fill_manual(values = c(
    ">3" = "#e76f51",
    "3" = "#f4a261",
    "2" = "#e9c46a",
    "1" = "#2a9d8f"
  )) +
  theme(
    axis.text.x = element_text(
      size = 6
    ),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 7),
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    axis.line = element_line(
      linetype = 1,
      color = "black",
      size = 0.3
    ),
    axis.ticks = element_line(
      linetype = 1,
      color = "black",
      size = 0.3
    ),
    legend.position = "bottom"
  ) +
  labs(fill = 'Clone size') +
  facet_grid(
    facets = Annotation_major_2 ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  )+
  theme(panel.spacing.x = unit(0.2, "lines"), strip.text = element_text(size = 6))

ggsave(paste0(fdir,"FigureS5E.pdf"),width = 2,height = 3.5)

##--- Figure 4C; Boxplot comparing the percentage of clonally expanded cells (defined as cells from clonotypes with clone size > 1) in IgA and IgG isotype ASCs between tumors and ANTs
keep_sample <- BCR_df %>%
  group_by(SampleID) %>%
  summarise(n = n()) %>%
  filter(n > 2)

plot_df <- BCR_df %>%
  mutate(c_call_new = case_when(c_call %in% c("IGHA1","IGHA2") ~ "IGHA", 
                                c_call %in% c("IGHG1","IGHG2","IGHG3","IGHG4") ~ "IGHG")) %>% 
  filter(Annotation_major %in% c("ASC")) %>% 
  filter(Tissue %in% c("Adjacent non-tumor", "Tumor"),
         SampleID %in% keep_sample$SampleID) %>%
  group_by(Tissue, SampleID, c_call_new) %>%
  summarise(
    clonal_cell = sum(clonal),
    clonal_percent = sum(clonal) / n() * 100,
    n = n()
  ) %>%
  filter(c_call_new %in% c("IGHG","IGHA"))

plot_df$Tissue = factor(plot_df$Tissue,levels = c(c("Adjacent non-tumor", "Tumor")))

ggplot(plot_df,
       aes(x = Tissue,
           y = clonal_percent,
           color = Tissue)) +
  geom_boxplot(outlier.color = NA,
               lwd = 0.3) +
  geom_jitter(size = 0.2,
              width = 0.2) +
  ggsignif::geom_signif(
    color = "black",
    comparisons = list(c("Tumor", "Adjacent non-tumor")),
    test = wilcox.test,
    step_increase = -0.2,
    textsize = 0.35*6
  ) +
  scale_color_manual(values = Tissue_color_panel) +
  xlab("") + ylab("Percentage of clonally expanded cells") +
  facet_wrap(~c_call_new)+
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    title = element_text(size = 7),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    strip.background = element_rect(color = NA, fill = NA),
    strip.text = element_text(size = 7),
    legend.position = "none"
  )

ggsave(paste0(fdir,"Figure4C.pdf"),width = 1.7,height = 1.8)

##--- Figure 4D; scatter plot showing the median percentage of clonally expanded cells in each Bm subset from tumors and ANTs
keep_sample <- BCR_df %>%
  group_by(SampleID) %>%
  summarise(n = n()) %>%
  filter(n > 2)

plot_df <- BCR_df %>%
  filter(Annotation_major_2 %in% c("Bm")) %>%
  filter(Tissue %in% c("Adjacent non-tumor", "Tumor"),
         SampleID %in% keep_sample$SampleID) %>%
  group_by(Tissue, SampleID, Annotation) %>%
  summarise(
    clonal_cell = sum(clonal),
    clonal_percent = sum(clonal) / n() * 100,
    n = n()
  ) 

plot_df = plot_df %>% ungroup() %>%
  group_by(Tissue, Annotation) %>%
  summarise(clonal_percent = median(clonal_percent)) %>%
  spread(Tissue, clonal_percent)

ggplot(plot_df, aes(
  x = `Adjacent non-tumor`,
  y = Tumor,
  color = Annotation,
  label = Annotation
)) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey"
  ) +
  geom_point() +
  ggrepel::geom_text_repel(size = 0.35 * 5) +
  scale_color_manual(values = Subset_color_panel, name = "") +
  cowplot::theme_cowplot() +
  xlim(-0.01,30)+
  ylim(-0.01,30)+
  ylab("Percentage of clonally expanded \n cells in tumor") +
  xlab("Percentage of clonally expanded \n cells in adjacent non-tumor") +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 7),
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    legend.position = "none"
  )+
  coord_fixed()

ggsave(paste0(fdir,"Figure4D.pdf"),width = 2.3, height = 2.3)

##--- Figure 4E; Boxplot comparing the percentage of clonally expanded cells in c08 cells between tumors and ANTs
keep_sample <- BCR_df %>%
  group_by(SampleID) %>%
  summarise(n = n()) %>%
  filter(n > 2)

plot_df <- BCR_df %>%
  filter(Tissue %in% c("Adjacent non-tumor", "Tumor"),
         SampleID %in% keep_sample$SampleID) %>%
  group_by(Tissue, SampleID, Annotation) %>%
  summarise(
    clonal_cell = sum(clonal),
    clonal_percent = sum(clonal) / n() * 100,
    n = n()
  ) %>%
  filter(Annotation %in% c("c08_ABC_FCRL4"))

plot_df$Tissue = factor(plot_df$Tissue,levels = c(c("Adjacent non-tumor", "Tumor")))

ggplot(plot_df,
       aes(x = Tissue,
           y = clonal_percent,
           color = Tissue)) +
  geom_boxplot(outlier.color = NA,
               lwd = 0.3) +
  geom_jitter(size = 0.2,
              width = 0.2) +
  ggsignif::geom_signif(
    color = "black",
    comparisons = list(c("Tumor", "Adjacent non-tumor")),
    test = wilcox.test,
    step_increase = -0.2,
    textsize = 0.35*6
  ) +
  scale_color_manual(values = Tissue_color_panel) +
  xlab("") + ylab("Percentage of clonally expanded cells") +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    title = element_text(size = 7),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    strip.background = element_rect(color = NA, fill = NA),
    strip.text = element_text(size = 7),
    legend.position = "none"
  )

ggsave(paste0(fdir,"Figure4E.pdf"),width = 1.1,height = 1.6)

#--- ASC part ---#

##--- Figure 4I; Distribution of IgH isotypes among c15_cycling_ASC
tmp = BCR_df %>% 
  filter(Tissue == "Tumor",Annotation == "c15_cycling_ASC") %>% 
  group_by(c_call) %>% summarise(n = n()) %>%
  mutate(percent = n/sum(n))

##--- Figure 4K; Boxplots comparing the composition of ASC IgH isotypes between tumors and ANTs
keep_sample <- BCR_df %>%
  group_by(SampleID) %>%
  summarise(n = n()) %>%
  filter(n > 2)

plot_df <- BCR_df %>%
  filter(Annotation_major_2 == "ASC",
         Tissue %in% c("Adjacent non-tumor", "Tumor"),
         SampleID %in% keep_sample$SampleID,
         c_call %in% c("IGHE")== FALSE,
         !is.na(c_call))%>%
  group_by(SampleID,Tissue,c_call) %>%
  summarise(ign = n()) %>% 
  group_by(SampleID) %>%
  mutate(samplen = sum(ign)) %>%
  mutate(percent = ign/samplen*100)

plot_df$Tissue = factor(plot_df$Tissue,levels = c(c("Adjacent non-tumor", "Tumor")))
plot_df$c_call = factor(plot_df$c_call,levels = c("IGHM","IGHD","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4"))

ggplot(plot_df, aes(x = Tissue, y = percent)) +
  geom_boxplot(aes(color = Tissue), outlier.colour = NA, lwd = 0.3) +
  geom_jitter(aes(color = Tissue), size = 0.4,shape = 16, stroke = 0, width = 0.2) +
  ggsignif::geom_signif(
    comparisons = list(c("Tumor", "Adjacent non-tumor")),
    test = wilcox.test,
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
    step_increase = -0.1,
    textsize = 2,
    y_position = 90,
    size = 0.2
  )+
  scale_color_manual(values = Tissue_color_panel) +
  facet_wrap(~ c_call, nrow = 1) +
  xlab("") +
  ylab("Percent of isotype") +
  theme_bw() +  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 7),
    panel.grid.major = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    legend.position = "none"
  ) 

ggsave(paste0(fdir,"Figure4K.pdf"),width = 4.1,height = 1.7)

##--- Figure 5L; BCR repertoire overlaps across B cells with different IgH isotypes in tumors and ANTs.
i_tissue = "Tumor"
i_tissue = "Adjacent non-tumor"

if (TRUE) {
  df = BCR_df %>%
    filter(Tissue == i_tissue,
           is.na(c_call) == FALSE,
           c_call %in% c("IGHE") == FALSE)
  
  df$c_call = ifelse(df$c_call %in% c("IGHM", "IGHD"), "IGHMD", df$c_call)
  
  m = matrix(nrow = length(unique(df$c_call)), ncol = length(unique(df$c_call)))
  rownames(m) = c("IGHMD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHA2")
  colnames(m) = c("IGHMD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHA2")
  
  for (i in sort(unique(df$c_call))) {
    for (j in sort(unique(df$c_call))) {
      if (i == j) {
        # clone which contain i
        tmp1 = df %>% filter(c_call == i)
        m[i, j] = length(unique(tmp1$clone_id))
      } else{
        # clone which contain i and j
        tmp1 = df %>% filter(c_call %in% c(i, j)) %>% group_by(clone_id) %>%
          summarise(ccall_n = length(unique(c_call))) %>% filter(ccall_n == 2)
        m[i, j] = length(unique(tmp1$clone_id))
      }
    }
  }
  
  m[lower.tri(m)] <- 0
  
  # m
  n = m
  for (i in 1:nrow(m)) {
    for (j in 1:nrow(m)) {
      if (i == j) {
        n[i, j] = m[i, j]
      } else{
        n[i, j] = m[i, j] / sum(m[i, i]) * 100
      }
    }
  }
  
  m = n
  col_fun <- colorRamp2(c(0, 3, 6, 9),
                        c("#fffde7", "#ffe0b2", "#ff9800", "#e65100"))
  
  m = t(m)
  
  pdf(paste0(fdir, "Figure4L_", i_tissue, ".pdf"),
      width = 5,
      height = 5)
  Heatmap(
    m,
    col = col_fun,
    cluster_rows = F,
    cluster_columns = F,
    column_names_gp = grid::gpar(fontsize = 6),
    row_names_gp = grid::gpar(fontsize = 6),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", m[i, j]), x, y, gp = gpar(fontsize = 6))
    },
    width = ncol(m) * unit(0.3, "inch"),
    height = nrow(m) * unit(0.3, "inch"),
    name = "Number of overlaped clones"
  )
  dev.off()
}

#----- 05. ASC expresssion -----#
##--- Figure S5J; CCR10 expression in c16 and c17 PCs.
require(Seurat)

tmp_obj = subset(obj, Annotation %in% c("c16_PC_IGHG","c17_PC_IGHA"))
Idents(tmp_obj) = tmp_obj$Annotation

VlnPlot(
  object = tmp_obj,
  cols = Subset_color_panel,
  features = c("CCR10"),
  pt.size = 0
) +
  xlab("") +
  ylab("CCR10 expression") +
  ggtitle("") +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 7),
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    legend.position = "none"
  )

ggsave(paste0(fdir,"FigureS5J.pdf"),width = 1.8,height = 2.5)

##--- Figure S5L; Expression of Fc receptors in TME immune subtypes. Clusters in which less than 15% of cell expressing any of the eight Fc receptors are not shown.
tmp_obj = subset(CD45_obj, Tissue_short == "T" & Annotation_CD45_major != "B")

FcR = c(c("FCGR1A", "FCGR2A", "FCGR2C", "FCGR3A", "FCGR3B"),
        c("FCGR2B"),
        c("FCAMR", "FCAR"))

yy_Dotplot(
  seuratObj = tmp_obj,
  genes = FcR,
  coord_flip = TRUE,
  group.by = "Annotation_CD45",
  dot.scale = 2.5,
  cell_pct_cutoff = 15
) +
  scale_color_gradientn(colors  = brewer.pal(9, "RdBu")[9:1]) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor  = element_blank(),
        axis.line = element_line(size = 0.3),
        axis.ticks = element_line(size = 0.3),
        axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 6))

ggsave(paste0(fdir,"FigureS5L.pdf"),width = 5.3, height = 2.9)

#----- 06. Cell typist c08 cells -----#

## prepare
library(SeuratDisk)
obj$Annotation <- as.character(obj$Annotation)

tmp_obj = subset(obj,Annotation_major == "Bm")
tmp_obj = DietSeurat(tmp_obj,counts = F,scale.data = F)
SaveH5Seurat(tmp_obj, filename = "reference_memory.h5Seurat",overwrite = T)
Convert("reference_memory.h5Seurat", dest = "h5ad",overwrite = T)

tmp_obj = subset(obj,annotation == "c14_Bm_activated-cycling")
tmp_obj = DietSeurat(tmp_obj,counts = F,scale.data = F)
SaveH5Seurat(tmp_obj, filename = "pro_memory.h5Seurat",overwrite = T)
Convert("pro_memory.h5Seurat", dest = "h5ad",overwrite = T)

## run celltypist
# model, "celltypist_model_from_all_memory.pkl"
# predict, "predict_pro_memory.csv"


##---- Figure 4F;	Distribution of CellTypist-assigned subset labels in activated cycling Bm cells
df_celltypist = read_csv("predict_pro_memory.csv")
df_celltypist$predicted_labels %>% table()

##---- Figure 4G; Boxplot comparing the proportions of cycling cells in all FCRL4+ Bm cells between tumors and ANTs
tmp = data.frame(predicted_labels = 3:10,
                 Annotation_cycling = paste("Cycling_",sort(unique(meta$Annotation))[4:11],sep = ""))

df_celltypist = df_celltypist %>% left_join(tmp)

# sum(df_celltypist$CellID %in% meta$CellID)
tmp = meta %>% group_by(CellID, SampleID) %>% summarise()
df_celltypist = df_celltypist %>% subset(select = -c(SampleID)) %>% left_join(tmp)

tmp = df_celltypist %>% filter(Annotation_cycling == "Cycling_c08_ABC_FCRL4") %>%
  group_by(SampleID) %>% summarise(cycling_c08_n = n())

c08_meta = meta %>% filter(Annotation == "c08_ABC_FCRL4") %>%
  group_by(SampleID, Tissue) %>% summarise(c08_n = n()) %>% left_join(tmp)
c08_meta[is.na(c08_meta)] <- 0
c08_meta = c08_meta %>%
  mutate(cell_n = sum(c08_n, cycling_c08_n)) %>%
  mutate(percent = cycling_c08_n / cell_n * 100)

plot_df <- c08_meta %>%
  filter(Tissue %in% c("Adjacent non-tumor", "Tumor"),
         SampleID %in% df.sample_count_flt$SampleID,
         cell_n > 4)

plot_df$Tissue = factor(plot_df$Tissue,levels = c(c("Adjacent non-tumor", "Tumor")))

ggplot(plot_df,
       aes(x = Tissue,
           y = percent,
           color = Tissue)) +
  geom_boxplot(outlier.color = NA,
               lwd = 0.2) +
  geom_jitter(size = 0.3,
              shape = 16,
              stroke = 0,
              width = 0.2) +
  ggsignif::geom_signif(
    color = "black",
    comparisons = list(c("Tumor", "Adjacent non-tumor")),
    test = wilcox.test,
    step_increase = -0.2,
    textsize = 0.35*6,
    size = 0.2
  ) +
  scale_color_manual(values = Tissue_color_panel) +
  xlab("") + ylab("Percentage of clonally expanded cells") +
  cowplot::theme_cowplot() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    title = element_text(size = 7),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    strip.background = element_rect(color = NA, fill = NA),
    strip.text = element_text(size = 7),
    legend.position = "none"
  )

ggsave(paste0(fdir,"Figure4G.pdf"),width = 1.1,height = 1.6)
                     
