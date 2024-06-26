## load packages
library(readr)
library(dplyr)
library(ggplot2)
library(scales)

## load data
obj = read_rds("All_obj.rds")
meta = obj@meta.data

## load path
fdir = ""

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
  ylab("Number of cells") + 
  cowplot::theme_cowplot() +
  scale_y_continuous(
    trans = "sqrt",
    breaks = c(1000, 10000, 50000, 100000, 200000),
    labels = comma
  ) +
  theme(
    axis.text.x = element_text(
      size = 7,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 7),
    text = element_text(size = 8),
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
    )
  )

ggsave(
  filename = paste0(
    fdir,"Figure1A_Cell_count.pdf"
  ),
  device = "pdf",
  width = 3.5,
  height = 2.5
)
