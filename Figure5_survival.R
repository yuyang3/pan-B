################## ------------------ Figures 5A-C and S7I ------------------ ##################
# 1. library
library(readr)
library(tidyverse)
library(Seurat)
library(foreach)
library(doParallel)
library(AUCell)
library(survival)
library(survminer)
library(meta)

# 2. params
dir_for_CD45_obj <- "data/obj_CD45_20240708.rds"
dir_for_B_obj <- "data/All_obj.rds"
dir_for_TCGA <- "data/bulk_analysis/04.TCGA_data/TCGA_data_tumor.rds"
dir_for_leukocyte_proportions <- "data/bulk_analysis/04.TCGA_data/TCGA_all_leuk_estimate.masked.20170107.tsv"
dir_for_cibersort <- "data/bulk_analysis/04.TCGA_data/TCGA.Kallisto.fullIDs.cibersort.relative.tsv"

dir_for_result <- "figures"
dir_for_data <- "data/bulk_analysis"

# 3. calculate DEGs for each B subset in the CD45 compartment
## load data
CD45_obj <- read_rds(dir_for_CD45_obj)
CD45_tumor <- subset(CD45_obj, subset = Tissue_short == "T")

## set new idents
Idents(CD45_tumor) <- CD45_tumor$Annotation_CD45
B_celltypes <- sort(unique(CD45_tumor$Annotation_CD45[CD45_tumor$Annotation_CD45_major == "B"]))

## calculate DEGs
dir.create(file.path(dir_for_data, "01.B_subset_signature_in_CD45"))
cl <- makeCluster(4)
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths(c("/home/cxy/R/x86_64-pc-linux-gnu-library/4.1", "/share/modules/R/R-4.1.2/lib/R/library")))
foreach(
    ident.1 = B_celltypes,
    .export = c("CD45_tumor", "dir_for_data"),
    .packages = c("Seurat", "readr", "tidyverse")
) %dopar% {
    markers <- Seurat::FindMarkers(CD45_tumor, ident.1 = ident.1, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox") %>%
        rownames_to_column("gene") %>%
        arrange(desc(avg_log2FC))
    write_tsv(markers, file.path(dir_for_data, "01.B_subset_signature_in_CD45", paste0(ident.1, "_wilcoxon_DE_genes.tsv")))
}
stopCluster(cl)

# 4. identify cell cycle-related genes
## calculate DEGs for cycling B subsets vs non-cycling B cells
### load data
B_obj <- read_rds(dir_for_B_obj)
B_obj_tumor <- subset(B_obj, subset = Tissue_short == "T" & Treatment_status == "treatment naïve")

### set new idents
B_obj_tumor$Annotation_tmp <- ifelse(B_obj_tumor$Annotation_major == "Bcycling", "cycling", "non-cycling")
B_obj_tumor$Annotation_tmp[B_obj_tumor$Annotation_tmp == "cycling"] <- as.character(B_obj_tumor$Annotation[B_obj_tumor$Annotation_tmp == "cycling"])
Idents(B_obj_tumor) <- B_obj_tumor$Annotation_tmp

### calculate DEGs
dir.create(file.path(dir_for_data, "02.cell_cycle_genes"))
cl <- makeCluster(3)
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths(c("/home/cxy/R/x86_64-pc-linux-gnu-library/4.1", "/share/modules/R/R-4.1.2/lib/R/library")))
foreach(
    ident.1 = setdiff(unique(B_obj_tumor$Annotation_tmp), "non-cycling"),
    .export = c("B_obj_tumor", "dir_for_data"),
    .packages = c("Seurat", "readr", "tidyverse")
) %dopar% {
    markers <- Seurat::FindMarkers(B_obj_tumor, ident.1 = ident.1, ident.2 = "non-cycling", only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox") %>%
        rownames_to_column("gene") %>%
        arrange(desc(avg_log2FC))
    write_tsv(markers, file.path(dir_for_data, "02.cell_cycle_genes", paste0(ident.1, "_wilcoxon_DE_genes.tsv")))
}
stopCluster(cl)

### select overlapping DEGs in three cycling B subsets
deg <- list()
deg_files <- list.files(file.path(dir_for_data, "02.cell_cycle_genes"))
deg_files <- deg_files[str_detect(deg_files, "_wilcoxon_DE_genes\\.tsv")]
for (i_file in deg_files) {
    markers_1 <- read_tsv(file.path(dir_for_data, "02.cell_cycle_genes", i_file))
    markers_1 <- markers_1 %>%
        dplyr::filter(
            p_val_adj < 0.05,
            avg_log2FC > 0.5
        ) %>%
        arrange(desc(avg_log2FC))
    markers_1 <- markers_1$gene
    deg[[i_file]] <- markers_1
}
cell_cycle_genes <- sort(Reduce(intersect, deg))

## combine with KEGG cell cycling genes
### KEGG_CELL_CYCLE: hsa04110
cell_cycle_genes_2 <- c("ABL1", "ANAPC1", "ANAPC10", "ANAPC11", "ANAPC13", "ANAPC2", "ANAPC4", "ANAPC5", "ANAPC7", "ATM", "ATR", "BUB1", "BUB1B", "BUB3", "CCNA1", "CCNA2", "CCNB1", "CCNB2", "CCNB3", "CCND1", "CCND2", "CCND3", "CCNE1", "CCNE2", "CCNH", "CDC14A", "CDC14B", "CDC16", "CDC20", "CDC23", "CDC25A", "CDC25B", "CDC25C", "CDC26", "CDC27", "CDC45", "CDC6", "CDC7", "CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDKN1A", "CDKN1B", "CDKN1C", "CDKN2A", "CDKN2B", "CDKN2C", "CDKN2D", "CHEK1", "CHEK2", "CREBBP", "CUL1", "DBF4", "E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "EP300", "ESPL1", "FZR1", "GADD45A", "GADD45B", "GADD45G", "GSK3B", "HDAC1", "HDAC2", "MAD1L1", "MAD2L1", "MAD2L2", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "MDM2", "MYC", "ORC1", "ORC2", "ORC3", "ORC4", "ORC5", "ORC6", "PCNA", "PKMYT1", "PLK1", "PRKDC", "PTTG1", "PTTG2", "RAD21", "RB1", "RBL1", "RBL2", "RBX1", "SFN", "SKP1", "SKP1P2", "SKP2", "SMAD2", "SMAD3", "SMAD4", "SMC1A", "SMC1B", "SMC3", "STAG1", "STAG2", "TFDP1", "TFDP2", "TGFB1", "TGFB2", "TGFB3", "TP53", "TTK", "WEE1", "WEE2", "YWHAB", "YWHAE", "YWHAG", "YWHAH", "YWHAQ", "YWHAZ", "ZBTB17")

cell_cycle_genes_combined <- sort(union(cell_cycle_genes, cell_cycle_genes_2))
write_rds(cell_cycle_genes_combined, file.path(dir_for_data, "02.cell_cycle_genes", "cell_cycle_genes_combined_with_KEGG.rds"), compress = "gz")

# 5. load TCGA data
TCGA_data <- read_rds(dir_for_TCGA)
expr <- TCGA_data$tpm
meta <- TCGA_data$meta
rm(TCGA_data)

## score B cell abundance by cibersort (B in CD45 pct (from cibersort) * CD45 in all pct (from methylation analysis))
## from https://gdc.cancer.gov/about-data/publications/panimmune
leukocyte_prop <- data.table::fread(dir_for_leukocyte_proportions)
cibersort <- data.table::fread(dir_for_cibersort)

leukocyte_prop$V2 <- lapply(leukocyte_prop$V2, function(x) substr(x, 1, 15))
leukocyte_prop <- leukocyte_prop %>% distinct(V2, .keep_all = TRUE)
sum(!meta$sample %in% leukocyte_prop$V2)
leukocyte_prop$V3[leukocyte_prop$V3 <= 0] <- 0

cibersort$SampleID <- lapply(cibersort$SampleID, function(x) substr(x, 1, 15))
cibersort$SampleID <- str_replace_all(cibersort$SampleID, "\\.", "-")
cibersort <- cibersort %>% distinct(SampleID, .keep_all = TRUE)
sum(!meta$sample %in% cibersort$SampleID)

meta$leukocyte_pct <- leukocyte_prop$V3[match(meta$sample, leukocyte_prop$V2)]
meta$B_pct_in_leukocyte <- rowSums(cibersort[match(meta$sample, cibersort$SampleID), c("B.cells.naive", "B.cells.memory", "Plasma.cells")])
meta$B_pct_in_all <- meta$B_pct_in_leukocyte * meta$leukocyte_pct

## build cells rankings for AUCell scoring in each cancer type
dir.create(file.path(dir_for_data, "04.TCGA_data", "cells_rankings_per_cancer"))
for (i_cancer in unique(meta$cancer.type.abbreviation)) {
    subset <- meta$cancer.type.abbreviation == i_cancer
    expr_subset <- expr[, subset]
    meta_subset <- meta[subset, ]
    cells_rankings <- AUCell::AUCell_buildRankings(as.matrix(expr_subset), nCores = 10, plotStats = FALSE)
    write_rds(cells_rankings, file.path(dir_for_data, "04.TCGA_data", "cells_rankings_per_cancer", paste0(i_cancer, ".rds")), compress = "gz")
}

## 6. survival function
run_survival <- function(
    category = NULL, # low, other, high
    survival_time,
    survival_status,
    correction_df = NULL,
    plot_forest = FALSE,
    plot_curve = FALSE,
    curve_line_size = 1,
    curve_censor_size = 4.5,
    return_meta_survival = FALSE) {
    require(survival)

    if (length(unique(category)) < 2) {
        stop("Not enough categories.\n")
    }
    meta_survival <- data.frame(
        category = category,
        survival_time = survival_time,
        survival_status = survival_status,
        stringsAsFactors = FALSE
    )
    # remove correction factor that only contains NA
    for (i in colnames(correction_df)) {
        if (sum(!is.na(correction_df[[i]])) == 0) {
            correction_df[[i]] <- NULL
        }
    }

    if (!is.null(correction_df)) {
        meta_survival <- cbind(meta_survival, correction_df)
    }
    meta_survival <- meta_survival[meta_survival$category != "other", ]
    meta_survival$category <-
        factor(meta_survival$category, levels = c("low", "high"))
    n_sample <- nrow(meta_survival)

    formula <-
        "survival::Surv(survival_time, survival_status) ~ category"
    if (!is.null(correction_df)) {
        for (i in colnames(correction_df)) {
            formula <- paste(formula, "+", i)
        }
    }
    formula <- as.formula(formula)

    res <- survival::coxph(formula, data = meta_survival)

    ret <- list()

    # save info for test variable
    result <-
        c(
            summary(res)$conf.int[1, ],
            summary(res)$coefficients[1, 3:5],
            n_sample
        )
    names(result) <-
        c(colnames(summary(res)$conf.int), "se(coef)", "z", "Pr(>|z|)", "n")

    ret$result <- result

    # save info for test variable and co-variables
    result_all <- cbind(
        summary(res)$conf.int,
        summary(res)$coefficients[, 3:5],
        rep(n_sample, n = nrow(summary(res)$conf.int))
    )
    colnames(result_all) <-
        c(colnames(summary(res)$conf.int), "se(coef)", "z", "Pr(>|z|)", "n")
    result_all <- result_all %>%
        as.data.frame() %>%
        rownames_to_column("Var")

    ret$result_all <- result_all

    if (plot_forest) {
        require(survminer)
        p_forest <- survminer::ggforest(res, data = meta_survival)
        ret$p_forest <- p_forest
    }

    if (plot_curve) {
        require(survminer)

        sfit <-
            survival::survfit(survival::Surv(survival_time, survival_status) ~ category,
                data = meta_survival
            )
        p_curve <- survminer::ggsurvplot(
            sfit,
            data = meta_survival,
            palette = c("#0173C1", "#EFBF02"),
            legend = "right",
            # legend.title = cancer,
            legend.labs = levels(meta_survival$category),
            censor.size = curve_censor_size,
            size = curve_line_size
        ) +
            labs(title = paste0(
                "p value = ",
                round(summary(res)$coefficients[1, 5], 4),
                "\n",
                # please check
                "HR (high vs. low) = ",
                round(summary(res)$coefficients[1, 2], 2) # please check
            )) +
            xlab("time (days)")
        ret$p_curve <- p_curve
    }

    if (return_meta_survival) {
        ret$meta_survival <- meta_survival
    }
    return(ret)
}

test_survival_correct_B_abundance <- function(genes,
                                              dir_for_cells_rankings,
                                              aucMaxRank_perc = 0.05,
                                              low_threshold = 0.5,
                                              high_threshold = 0.5,
                                              cancer_restricted = TRUE, # only plot cancers in our pan-cancer atlas
                                              year_cutoff = NULL) {
    hr_all_all <- list() # all cancers, all co-variables

    if (cancer_restricted) {
        cancers <- sort(unique(meta$cancer.type.abbreviation[meta$Cancer != "Other"]))
    } else {
        cancers <- sort(unique(meta$cancer.type.abbreviation))
    }

    ### 1. calculate per cancer
    for (i_cancer in cancers) {
        print(i_cancer)

        subset <- meta$cancer.type.abbreviation == i_cancer

        if (sum(subset) < 100) {
            next
        }
        meta_subset <- meta[subset, ]

        if (!is.null(year_cutoff)) {
            use_sample <- meta_subset$OS.time <= year_cutoff * 365
            use_sample[is.na(use_sample)] <- FALSE
            meta_subset <- meta_subset[use_sample, ]
        }

        ##################################################################################################################
        # category
        cells_rankings <- read_rds(file.path(dir_for_cells_rankings, paste0(i_cancer, ".rds")))
        cells_AUC <- AUCell::AUCell_calcAUC(list("1" = genes),
            cells_rankings,
            nCores = 10,
            aucMaxRank = ceiling(aucMaxRank_perc * nrow(cells_rankings))
        )
        abundance <- cells_AUC@assays@data$AUC

        if (!is.null(year_cutoff)) {
            abundance <- abundance[use_sample]
        }

        threshold <- quantile(abundance,
            probs = c(low_threshold, high_threshold),
            na.rm = TRUE
        )
        category <- rep("other", length(abundance))
        category[abundance <= threshold[1]] <- "low"
        category[abundance > threshold[2]] <- "high"

        ##################################################################################################################
        # correction df
        correction_df <- meta_subset[, c("age_at_initial_pathologic_diagnosis", "Gender", "Tumor_stage", "B_pct_in_all")]
        colnames(correction_df) <- c("Age", "Sex", "Tumor_stage", "B_abundance_in_all_cells")

        correction_df[[4]] <- ifelse(correction_df[[4]] > median(correction_df[[4]], na.rm = TRUE), "high", "low")
        correction_df[[4]] <- factor(correction_df[[4]], levels = c("low", "high"))
        ##################################################################################################################

        res_cancer <- run_survival(
            category = category,
            survival_time = meta_subset$OS.time,
            survival_status = meta_subset$OS,
            correction_df = correction_df,
            plot_forest = FALSE,
            plot_curve = FALSE,
            curve_censor_size = 2,
            curve_line_size = 0.5
        )

        # save hr_all_all
        result_all <- res_cancer$result_all
        result_all$cancer <- i_cancer
        hr_all_all[[i_cancer]] <- result_all
    }

    hr_all_all <- dplyr::bind_rows(hr_all_all)

    ### 2. derive pan cancer HR by meta analysis
    hr_plot_df_all <- list()
    for (i_var in unique(hr_all_all$Var)) { # for each co-variable
        hr_all_all_tmp <- hr_all_all %>% dplyr::filter(Var == i_var)
        # metagen(logHR, selogHR, sm="HR")
        tmp <- meta::metagen(log(hr_all_all_tmp$`exp(coef)`), hr_all_all_tmp$`se(coef)`, sm = "HR", fixed = FALSE)
        pan_cancer_hr <- data.frame(
            `exp(coef)` = exp(tmp$TE.random),
            `lower .95` = exp(tmp$lower.random),
            `upper .95` = exp(tmp$upper.random),
            `Pr(>|z|)` = tmp$pval.random,
            cancer = "Pan cancer",
            n = sum(hr_all_all_tmp$n),
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
        hr_plot_df_all_tmp <- hr_all_all_tmp %>%
            dplyr::select(`exp(coef)`, `lower .95`, `upper .95`, `Pr(>|z|)`, cancer, n) %>%
            dplyr::filter(is.finite(`exp(coef)`) & is.finite(`lower .95`) & is.finite(`upper .95`) & is.finite(`Pr(>|z|)`)) %>%
            rbind(pan_cancer_hr) %>%
            mutate(Var = i_var)
        hr_plot_df_all[[i_var]] <- hr_plot_df_all_tmp
    }
    hr_plot_df_all <- bind_rows(hr_plot_df_all)

    # significance
    # hr_plot_df_all$significance <- "ns"
    hr_plot_df_all$significance <- ""
    hr_plot_df_all$significance[hr_plot_df_all$`Pr(>|z|)` <= 0.05] <- "*"
    hr_plot_df_all$significance[hr_plot_df_all$`Pr(>|z|)` <= 0.01] <- "**"
    hr_plot_df_all$significance[hr_plot_df_all$`Pr(>|z|)` <= 0.001] <- "***"
    hr_plot_df_all$significance[hr_plot_df_all$`Pr(>|z|)` <= 0.0001] <- "****"

    ### 3. plot HR (only for first variable)
    i_var <- unique(hr_all_all$Var)[1]
    hr_plot_df <- hr_plot_df_all[hr_plot_df_all$Var == i_var, ]

    hr_plot_df$color <- "Neutral"
    hr_plot_df$color[hr_plot_df$`exp(coef)` < 1] <- "Better survival"
    hr_plot_df$color[hr_plot_df$`exp(coef)` > 1] <- "Worse survival"
    hr_plot_df$color <- factor(hr_plot_df$color)
    mycolor <- c("Neutral" = "grey", "Worse survival" = "#D83215", "Better survival" = "#0F7B9F")

    # hr_plot_df$xlab <- paste0(hr_plot_df$cancer," (n=",hr_plot_df$n,")")

    levels <- levels(with(hr_plot_df, reorder(cancer, `exp(coef)`, function(x) -x)))
    levels <- c(levels[levels != "Pan cancer"], "Pan cancer")
    hr_plot_df$cancer <- factor(hr_plot_df$cancer, levels = levels)

    p <- ggplot(data = hr_plot_df, aes(x = cancer, y = `exp(coef)`, ymin = `lower .95`, ymax = `upper .95`, color = color)) +
        geom_pointrange() +
        geom_hline(yintercept = 1, lty = 2) + # add a dotted line at x=1 after flip
        geom_text(aes(x = cancer, y = max(`upper .95`), label = significance), show.legend = FALSE) +
        scale_color_manual(values = mycolor, name = "Survival association") +
        xlab("") +
        ylab("Hazard ratio (95% CI)") +
        cowplot::theme_cowplot() +
        theme(
            text = element_text(size = 7, family = "ArialMT"),
            axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 7),
            plot.title = element_text(hjust = 0.5, size = 8),
            plot.margin = unit(c(1, 1, 1, 1), "char"),
            axis.line = element_line(linetype = 1, color = "black", size = 0.3),
            axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
        ) +
        # scale_x_discrete(labels = hr_plot_df$xlab) +
        scale_y_continuous(trans = "log2")
    return(list(p = p, hr_df = hr_plot_df, hr_df_all = hr_plot_df_all))
}

# 7. run survival
## params
top_n_DE_genes <- 20
year_cutoff <- 10
dir_for_signatures <- file.path(dir_for_data, "01.B_subset_signature_in_CD45")
remove_cell_cycle_genes <- TRUE
cell_cycle_genes <- read_rds(file.path(dir_for_data, "02.cell_cycle_genes", "cell_cycle_genes_combined_with_KEGG.rds"))

## run
result_df <- list()
for (i_file in list.files(dir_for_signatures)) {
    # get signature with cell cycle genes removed
    markers_1 <- read_tsv(file.path(dir_for_signatures, i_file))
    if (remove_cell_cycle_genes) {
        markers_1 <- markers_1 %>%
            dplyr::filter(
                p_val_adj < 0.05,
                gene %in% rownames(expr),
                !gene %in% cell_cycle_genes
            ) %>%
            arrange(desc(avg_log2FC)) %>%
            slice_head(n = top_n_DE_genes)
    } else {
        markers_1 <- markers_1 %>%
            dplyr::filter(
                p_val_adj < 0.05,
                gene %in% rownames(expr)
            ) %>%
            arrange(desc(avg_log2FC)) %>%
            slice_head(n = top_n_DE_genes)
    }
    markers_1 <- markers_1$gene

    # plot all cancers
    result_df_tmp <- test_survival_correct_B_abundance(
        genes = markers_1,
        dir_for_cells_rankings = file.path(dir_for_data, "04.TCGA_data", "cells_rankings_per_cancer"),
        aucMaxRank_perc = 0.05,
        low_threshold = 0.5,
        high_threshold = 0.5,
        cancer_restricted = TRUE,
        year_cutoff = year_cutoff
    )
    result_df_tmp <- result_df_tmp$hr_df
    result_df_tmp$celltype <- str_extract(i_file, "^.*(?=_wilcoxon)")

    result_df[[i_file]] <- result_df_tmp
}
result_df <- bind_rows(result_df)
write_rds(result_df, file.path(dir_for_data, "TCGA_survival_result.rds"), compress = "gz")

# 8. plot pan-cancer survival for each B subset
plot_df <- result_df
## arrange celltype
plot_df <- plot_df %>%
    dplyr::filter(cancer == "Pan cancer") %>%
    arrange(`exp(coef)`)
plot_df$celltype <- factor(plot_df$celltype, levels = plot_df$celltype)
## color
plot_df$color <- as.character(plot_df$color)
plot_df$color[plot_df$`Pr(>|z|)` >= 0.05 & plot_df$`exp(coef)` < 1] <- "Better survival (P > 0.05)"
plot_df$color[plot_df$`Pr(>|z|)` >= 0.05 & plot_df$`exp(coef)` > 1] <- "Worse survival (P > 0.05)"
plot_df$color <- factor(plot_df$color, levels = c(
    "Better survival", "Better survival (P > 0.05)",
    "Worse survival (P > 0.05)", "Worse survival"
))
mycolor <- c(
    "Better survival" = "#0F7B9F",
    "Better survival (P > 0.05)" = "#E0F3F8FF",
    "Worse survival (P > 0.05)" = "#FDDBC7FF",
    "Worse survival" = "#C3423F"
)
## plot
ggplot(data = plot_df, aes(x = celltype, y = `exp(coef)`, ymin = `lower .95`, ymax = `upper .95`, color = color)) +
    geom_pointrange() +
    geom_hline(yintercept = 1, colour = "grey40", linetype = "dashed", size = 0.2) + # add a dotted line at x=1 after flip
    geom_text(aes(x = celltype, y = max(`upper .95`), label = significance), show.legend = FALSE) +
    scale_color_manual(values = mycolor, name = "Survival association") +
    xlab("") +
    ylab("Hazard ratio (95% CI)") +
    cowplot::theme_cowplot() +
    theme(
        text = element_text(size = 7, family = "ArialMT"),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, size = 8, face = "plain"),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
    ) +
    ggtitle("Pan-cancer prognostic associations")
ggsave(
    filename = file.path(dir_for_result, "5A.survival_pan_cancer.pdf"),
    height = 2.6,
    width = 5
)

# 9. plot per-cancer survival for TAABs
plot_df <- result_df[result_df$celltype == "c08_ABC_FCRL4", ]
## arrange cancer
plot_df <- plot_df %>%
    arrange(`exp(coef)`)
plot_df$cancer <- factor(plot_df$cancer, levels = c(setdiff(plot_df$cancer, "Pan cancer"), "Pan cancer"))
## color
plot_df$color <- as.character(plot_df$color)
plot_df$color[plot_df$`Pr(>|z|)` >= 0.05 & plot_df$`exp(coef)` < 1] <- "Better survival (P > 0.05)"
plot_df$color[plot_df$`Pr(>|z|)` >= 0.05 & plot_df$`exp(coef)` > 1] <- "Worse survival (P > 0.05)"
plot_df$color <- factor(plot_df$color, levels = c(
    "Better survival", "Better survival (P > 0.05)",
    "Worse survival (P > 0.05)", "Worse survival"
))
mycolor <- c(
    "Better survival" = "#0F7B9F",
    "Better survival (P > 0.05)" = "#E0F3F8FF",
    "Worse survival (P > 0.05)" = "#FDDBC7FF",
    "Worse survival" = "#C3423F"
)
## plot
ggplot(data = plot_df, aes(x = cancer, y = `exp(coef)`, ymin = `lower .95`, ymax = `upper .95`, color = color)) +
    geom_pointrange() +
    geom_hline(yintercept = 1, colour = "grey40", linetype = "dashed", size = 0.2) + # add a dotted line at x=1 after flip
    geom_text(aes(x = cancer, y = max(`upper .95`), label = significance), show.legend = FALSE) +
    scale_color_manual(values = mycolor, name = "Survival association") +
    xlab("") +
    ylab("Hazard ratio (95% CI)") +
    cowplot::theme_cowplot() +
    theme(
        text = element_text(size = 7, family = "ArialMT"),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, size = 8),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
    ) +
    scale_y_continuous(trans = "log2")
ggsave(
    filename = file.path(dir_for_result, "5B.survival_TAAB.pdf"),
    height = 2.1,
    width = 5
)

# 10. plot per-cancer survival for TAABs
plot_df <- result_df[result_df$celltype == "c06_Bm_stress-response", ]
## arrange cancer
plot_df <- plot_df %>%
    arrange(`exp(coef)`)
plot_df$cancer <- factor(plot_df$cancer, levels = c(setdiff(plot_df$cancer, "Pan cancer"), "Pan cancer"))
## color
plot_df$color <- as.character(plot_df$color)
plot_df$color[plot_df$`Pr(>|z|)` >= 0.05 & plot_df$`exp(coef)` < 1] <- "Better survival (P > 0.05)"
plot_df$color[plot_df$`Pr(>|z|)` >= 0.05 & plot_df$`exp(coef)` > 1] <- "Worse survival (P > 0.05)"
plot_df$color <- factor(plot_df$color, levels = c(
    "Better survival", "Better survival (P > 0.05)",
    "Worse survival (P > 0.05)", "Worse survival"
))
mycolor <- c(
    "Better survival" = "#0F7B9F",
    "Better survival (P > 0.05)" = "#E0F3F8FF",
    "Worse survival (P > 0.05)" = "#FDDBC7FF",
    "Worse survival" = "#C3423F"
)
## plot
ggplot(data = plot_df, aes(x = cancer, y = `exp(coef)`, ymin = `lower .95`, ymax = `upper .95`, color = color)) +
    geom_pointrange() +
    geom_hline(yintercept = 1, colour = "grey40", linetype = "dashed", size = 0.2) + # add a dotted line at x=1 after flip
    geom_text(aes(x = cancer, y = max(`upper .95`), label = significance), show.legend = FALSE) +
    scale_color_manual(values = mycolor, name = "Survival association") +
    xlab("") +
    ylab("Hazard ratio (95% CI)") +
    cowplot::theme_cowplot() +
    theme(
        text = element_text(size = 7, family = "ArialMT"),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, size = 8),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
    ) +
    scale_y_continuous(trans = "log2")
ggsave(
    filename = file.path(dir_for_result, "5C.survival_stress_response_Bm.pdf"),
    height = 2.1,
    width = 5
)

# 11. plot survival for TAAB vs total B abundance
## run survival
markers_1 <- read_tsv(file.path(dir_for_signatures, "c08_ABC_FCRL4_wilcoxon_DE_genes.tsv"))
if (remove_cell_cycle_genes) {
    markers_1 <- markers_1 %>%
        dplyr::filter(
            p_val_adj < 0.05,
            gene %in% rownames(expr),
            !gene %in% cell_cycle_genes
        ) %>%
        arrange(desc(avg_log2FC)) %>%
        slice_head(n = top_n_DE_genes)
} else {
    markers_1 <- markers_1 %>%
        dplyr::filter(
            p_val_adj < 0.05,
            gene %in% rownames(expr)
        ) %>%
        arrange(desc(avg_log2FC)) %>%
        slice_head(n = top_n_DE_genes)
}
markers_1 <- markers_1$gene
result_df_tmp <- test_survival_correct_B_abundance(
    genes = markers_1,
    dir_for_cells_rankings = file.path(dir_for_data, "04.TCGA_data", "cells_rankings_per_cancer"),
    aucMaxRank_perc = 0.05,
    low_threshold = 0.5,
    high_threshold = 0.5,
    cancer_restricted = TRUE,
    year_cutoff = year_cutoff
)
## plot
plot_df <- result_df_tmp$hr_df_all
plot_df <- plot_df %>%
    dplyr::filter(Var %in% c("categoryhigh", "B_abundance_in_all_cellshigh"))
plot_df$Var[plot_df$Var == "categoryhigh"] <- "TAAB"
plot_df$Var[plot_df$Var == "B_abundance_in_all_cellshigh"] <- "B cell"

plot_df$sig <- NA
plot_df$sig[plot_df$`Pr(>|z|)` < 0.05 & plot_df$`exp(coef)` < 1] <- "blue"
plot_df$sig[plot_df$`Pr(>|z|)` < 0.05 & plot_df$`exp(coef)` > 1] <- "red"

plot_df$logHR <- log(plot_df$`exp(coef)`)

plot_df_tmp <- plot_df %>%
    dplyr::filter(Var == "TAAB") %>%
    arrange(`exp(coef)`)

plot_df$cancer <- factor(plot_df$cancer, levels = c(setdiff(plot_df_tmp$cancer, "Pan cancer"), "Pan cancer"))

ggplot(plot_df, aes(x = cancer, y = Var)) +
    geom_tile(aes(fill = logHR, width = 0.9, height = 0.9), color = plot_df$sig, size = 0.35) +
    coord_equal() +
    scale_fill_gradient2(high = "#D83215", low = "#0F7B9F", name = "Log(Hazard ratio)") +
    cowplot::theme_cowplot() +
    theme(
        panel.grid = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        text = element_text(size = 7, family = "ArialMT"),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        plot.title = element_text(hjust = 0.5, size = 8, family = "ArialMT", face = "plain"),
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
    ) +
    xlab("") +
    ylab("")

ggsave(file.path(dir_for_result, "S7I.survival_TAAB_vs_total_B_abundance.pdf"),
    width = 4.5,
    height = 2
)
