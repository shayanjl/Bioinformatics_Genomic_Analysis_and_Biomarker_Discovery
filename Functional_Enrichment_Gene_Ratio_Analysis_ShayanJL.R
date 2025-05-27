# 1. Set working directory
# setwd("D:/XlocationX")

# 2. Install required packages if not already installed
required_pkgs <- c("readxl", "dplyr", "ggplot2", "forcats", "scales")
new_pkgs <- setdiff(required_pkgs, installed.packages()[,"Package"])
if(length(new_pkgs)) install.packages(new_pkgs)

# 3. Load libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(forcats)
library(scales)

# 4. Read the cleaned Excel file for dataset 200
df <- read_excel("200_cleaned_overlap.xlsx")

# 5. Compute gene ratio and sort by Gene Ratio
df <- df %>%
  mutate(Gene_Ratio = Overlap_Count / Overlap_Total) %>%
  arrange(Gene_Ratio)

# 6. Create categories of adjusted p-value based on its quartiles
breaks <- quantile(df$`Adjusted P-value`, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
labels <- c(
  paste0("≤", signif(breaks[2], 2)),
  paste0(signif(breaks[2], 2), "–", signif(breaks[3], 2)),
  paste0(signif(breaks[3], 2), "–", signif(breaks[4], 2)),
  paste0(">", signif(breaks[4], 2))
)
df <- df %>%
  mutate(p_adj_cat = cut(
    `Adjusted P-value`,
    breaks = breaks,
    labels = labels,
    include.lowest = TRUE
  ))

# 7. Define a 4-color palette (red → magenta → blue → darkblue)
pal <- c("red", "magenta", "blue", "darkblue")

# 8. Plot with discrete legend for p.adjust categories
p <- ggplot(df, aes(
  x    = Gene_Ratio,
  y    = fct_reorder(Term, Gene_Ratio),
  fill = p_adj_cat
)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = pal,
    name   = "p.adjust",
    drop   = FALSE
  ) +
  labs(
    x     = "Gene Ratio",
    y     = NULL,
    title = "Functional Enrichment by Gene Ratio"
  ) +
  theme_minimal(base_size = 8, base_family = "Times New Roman") +
  theme(
    plot.title       = element_text(size = 10, face = "bold"),
    axis.text.y      = element_text(size = 4),
    axis.text.x      = element_text(size = 6),
    axis.title.x     = element_text(size = 8),
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 6),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )

# 9. Draw and save the plot
print(p)
ggsave(
  "new-functional_enrichment_200_by_gene_ratio_categorized_padj.jpg",
  plot = p,
  width = 6, height = 10, dpi = 300
)
