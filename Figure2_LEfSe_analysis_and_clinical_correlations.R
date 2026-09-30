# ==============================================================================
# Figure 2: Oropharyngeal dysbiosis and clinical correlations in ILD subgroups
#
# This script generates:
#   Figure 2A: LEfSe differential species (ILD-only vs Infection-related ILD)
#   Figure 2B: correlations between LEfSe-identified species and clinical indices
#   Figure 2C: correlations between top15 species and clinical indices
#
# Output: output/Figure2/
# ==============================================================================

# ------------------------------ 0. Packages ----------------------------------
library(tidyverse)
library(microeco)
library(ggplot2)
library(reshape2)
library(dplyr)
library(Hmisc)
library(readr)

# ------------------------------ 1. Paths -------------------------------------
data_dir   <- "data"
output_dir <- "output/Figure2"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# Figure 2A: LEfSe differential species
# ==============================================================================

# ----------------------------- Read data -------------------------------------
counts_data <- read.table(
  file.path(data_dir, "species_counts_filtered.tsv"),
  header = TRUE, sep = "\t", check.names = FALSE,
  stringsAsFactors = FALSE
)

metadata <- read.csv(
  file.path(data_dir, "metadata.csv"),
  header = TRUE, stringsAsFactors = FALSE
)

rownames(counts_data) <- counts_data$name
otu_table <- counts_data[, -c(1, 2, 3)]

metadata$SampleID <- as.character(metadata$SampleID)

common_samples <- intersect(colnames(otu_table), metadata$SampleID)
otu_table <- otu_table[, common_samples, drop = FALSE]
metadata  <- metadata[match(common_samples, metadata$SampleID), ]

target_groups <- c("ILD-only", "Infection-related ILD")
metadata  <- metadata[metadata$SubGroup %in% target_groups, ]
otu_table <- otu_table[, metadata$SampleID, drop = FALSE]

# ----------------------------- CPM normalization -----------------------------
otu_table_normalized <- sweep(otu_table, 2, colSums(otu_table), "/") * 1e6
otu_table_normalized <- as.data.frame(otu_table_normalized)

sample_table <- metadata
rownames(sample_table) <- sample_table$SampleID

tax_table <- data.frame(
  Species = rownames(otu_table_normalized),
  row.names = rownames(otu_table_normalized)
)

dataset <- microtable$new(
  otu_table = otu_table_normalized,
  sample_table = sample_table,
  tax_table = tax_table
)

# ----------------------------- LEfSe ----------------------------------------
lefse <- trans_diff$new(
  dataset = dataset,
  method = "lefse",
  group = "SubGroup",
  taxa_level = "Species",
  alpha = 0.05,
  lda_score = 3,
  p_adjust_method = "BH"
)

res <- lefse$res_diff
lda_col <- grep("lda", colnames(res), ignore.case = TRUE, value = TRUE)[1]
res_filtered <- res[res[[lda_col]] >= 3, ]

lefse_out <- file.path(output_dir, "LEfSe_ILD-only_vs_Infection-related_ILD_LDA3.csv")
write.csv(res_filtered, lefse_out, row.names = FALSE)

# ----------------------------- Bar plot -------------------------------------
plot_data <- res_filtered %>%
  arrange(desc(LDA)) %>%
  mutate(
    Taxa  = factor(Taxa, levels = rev(Taxa)),
    Group = factor(Group, levels = c("ILD-only", "Infection-related ILD"))
  )

unified_palette <- c(
  "ILD-only"              = "#0077BB",
  "Infection-related ILD" = "#CC3311"
)

p_final <- ggplot(plot_data, aes(x = LDA, y = Taxa, fill = Group)) +
  geom_bar(stat = "identity", width = 0.72, color = NA) +
  scale_fill_manual(values = unified_palette) +
  labs(x = "LDA Score", y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text.y = element_text(size = 9, face = "italic", color = "black"),
    axis.text.x = element_text(size = 9.5, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 9.5),
    legend.background = element_blank()
  )

height_use <- max(4, min(12, nrow(plot_data) * 0.22 + 1.5))
file_base <- file.path(output_dir, "Figure2A_LEfSe_LDA_barplot")

ggsave(paste0(file_base, ".pdf"), p_final,
       device = cairo_pdf, width = 8.5, height = height_use)
ggsave(paste0(file_base, ".png"), p_final,
       width = 8.5, height = height_use, dpi = 300)
ggsave(paste0(file_base, ".tiff"), p_final,
       width = 8.5, height = height_use, dpi = 600, compression = "lzw")

# ==============================================================================
# Figure 2B: correlations between LEfSe-identified species and clinical indices
# ==============================================================================

# ----------------------------- Read data -------------------------------------
otu_file   <- file.path(data_dir, "species_relative_abundance.tsv")
meta_file  <- file.path(data_dir, "metadata.csv")

otu   <- read.table(otu_file, header = TRUE, sep = "\t",
                    row.names = 1, check.names = FALSE)
meta  <- read.csv(meta_file, stringsAsFactors = FALSE)
lefse <- read.csv(lefse_out)

diff_species <- lefse$Taxa
otu <- otu[diff_species, , drop = FALSE]

common <- intersect(colnames(otu), meta$SampleID)
otu  <- otu[, common]
meta <- meta[match(common, meta$SampleID), ]

clin <- c("WBC_count", "Neu_count", "Lym_count",
          "Eos_count", "CRP", "OxygenationIndex")

rho_mat <- matrix(NA, nrow = nrow(otu), ncol = length(clin))
p_mat   <- matrix(NA, nrow = nrow(otu), ncol = length(clin))

rownames(rho_mat) <- rownames(otu)
colnames(rho_mat) <- clin
rownames(p_mat)   <- rownames(otu)
colnames(p_mat)   <- clin

for (i in rownames(otu)) {
  for (j in clin) {
    test <- cor.test(as.numeric(otu[i, ]), meta[[j]], method = "spearman")
    rho_mat[i, j] <- test$estimate
    p_mat[i, j]   <- test$p.value
  }
}

p_adj <- apply(p_mat, 2, p.adjust, method = "BH")

write.csv(rho_mat, file.path(output_dir, "Figure2B_LEfSe_species_clinical_rho.csv"))
write.csv(p_mat,   file.path(output_dir, "Figure2B_LEfSe_species_clinical_p.csv"))
write.csv(p_adj,   file.path(output_dir, "Figure2B_LEfSe_species_clinical_padj.csv"))

# ----------------------------- Heatmap --------------------------------------
rho  <- rho_mat
padj <- p_adj

colnames(rho)[colnames(rho) == "OxygenationIndex"]   <- "Oxygenation_Index"
colnames(padj)[colnames(padj) == "OxygenationIndex"] <- "Oxygenation_Index"

rho_long  <- melt(as.matrix(rho),  varnames = c("Species", "Clinical"),
                  value.name = "rho")
padj_long <- melt(as.matrix(padj), varnames = c("Species", "Clinical"),
                  value.name = "FDR")

df <- left_join(rho_long, padj_long, by = c("Species", "Clinical"))

df <- df %>%
  mutate(sig = case_when(
    FDR < 0.001 ~ "***",
    FDR < 0.01  ~ "**",
    FDR < 0.05  ~ "*",
    TRUE        ~ ""
  ))

target_order <- lefse %>%
  arrange(desc(LDA)) %>%
  pull(Taxa) %>%
  rev()

target_order <- intersect(target_order, rownames(rho))
df$Species <- factor(df$Species, levels = target_order)

color_low  <- "#3B7A9E"
color_mid  <- "#FFFFFF"
color_high <- "#B84A39"

p <- ggplot(df, aes(x = Clinical, y = Species)) +
  geom_tile(aes(fill = rho), color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = color_low, mid = color_mid, high = color_high,
    midpoint = 0, limits = c(-1, 1),
    name = "Correlation coefficient (rho)"
  ) +
  geom_text(aes(label = sig), color = "black", size = 4.5, vjust = 0.75) +
  scale_x_discrete(
    expand = c(0, 0),
    labels = c("Oxygenation_Index" = expression(paste("PaO"[2], "/", "FiO"[2])))
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_void(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               color = "black", size = 11),
    axis.text.y = element_text(face = "italic", color = "black", size = 10),
    legend.title = element_text(size = 11, face = "plain"),
    legend.text = element_text(size = 10),
    legend.position = "right"
  ) +
  labs(x = NULL, y = NULL)

fig_width  <- 7.5
fig_height <- max(6.0, nrow(rho) * 0.25 + 2.0)

ggsave(file.path(output_dir, "Figure2B_LEfSe_clinical_correlation_heatmap.pdf"),
       p, device = cairo_pdf, width = fig_width, height = fig_height, bg = "white")
ggsave(file.path(output_dir, "Figure2B_LEfSe_clinical_correlation_heatmap.png"),
       p, width = fig_width, height = fig_height, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "Figure2B_LEfSe_clinical_correlation_heatmap.tiff"),
       p, width = fig_width, height = fig_height, dpi = 600,
       compression = "lzw", bg = "white")

# ==============================================================================
# Figure 2C: correlations between top15 species and clinical indices
# ==============================================================================

# ----------------------------- Read data -------------------------------------
abundance_file   <- file.path(data_dir, "species_counts_filtered.tsv")
top15_stats_file <- file.path(data_dir, "SCI_Table_Species_Top15_Abundance.tsv")
metadata_file    <- file.path(data_dir, "metadata.csv")

counts_data <- read_tsv(abundance_file, col_names = TRUE, show_col_types = FALSE)
sample_cols <- setdiff(colnames(counts_data),
                       c("name", "taxonomy_id", "taxonomy_lvl"))

abundance <- counts_data %>%
  select(name, all_of(sample_cols)) %>%
  mutate(across(all_of(sample_cols), function(x) {
    total <- sum(x, na.rm = TRUE)
    if (total == 0) return(0) else return(x / total)
  }))

top15_table <- read_tsv(top15_stats_file, show_col_types = FALSE)
top_species <- top15_table[[1]]
top_species <- top_species[top_species != "Others" & !is.na(top_species)]

metadata <- read.csv(metadata_file, header = TRUE, stringsAsFactors = FALSE)

abundance_long <- abundance %>%
  pivot_longer(cols = all_of(sample_cols),
               names_to = "sample", values_to = "abundance") %>%
  mutate(sample = as.character(sample))

metadata <- metadata %>% mutate(SampleID = as.character(SampleID))

data_combined <- abundance_long %>%
  left_join(metadata, by = c("sample" = "SampleID")) %>%
  filter(!is.na(Batch))

existing_species <- intersect(top_species, unique(data_combined$name))
top_species <- existing_species

species_wide <- data_combined %>%
  filter(name %in% top_species) %>%
  select(sample, name, abundance) %>%
  pivot_wider(id_cols = sample, names_from = name,
              values_from = abundance, values_fill = 0)

clinical_cols <- c("WBC_count", "Neu_count", "Lym_count",
                   "Eos_count", "CRP", "OxygenationIndex")

clinical_data <- data_combined %>%
  select(sample, all_of(clinical_cols)) %>%
  distinct(sample, .keep_all = TRUE)

cor_data <- species_wide %>%
  left_join(clinical_data, by = "sample") %>%
  as.data.frame()
rownames(cor_data) <- cor_data$sample

species_mat  <- as.matrix(cor_data[, top_species, drop = FALSE])
clinical_mat <- as.matrix(cor_data[, clinical_cols, drop = FALSE])
all_mat      <- cbind(species_mat, clinical_mat)

cor_result <- Hmisc::rcorr(all_mat, type = "spearman")
R <- cor_result$r
P <- cor_result$P

R_species_clinical <- R[top_species, clinical_cols, drop = FALSE]
P_species_clinical <- P[top_species, clinical_cols, drop = FALSE]

P_adj_matrix <- apply(P_species_clinical, 2,
                      function(x) p.adjust(x, method = "fdr"))
rownames(P_adj_matrix) <- rownames(P_species_clinical)

write.csv(R_species_clinical,
          file.path(output_dir, "Figure2C_top15_clinical_rho.csv"),
          row.names = TRUE)
write.csv(P_species_clinical,
          file.path(output_dir, "Figure2C_top15_clinical_p.csv"),
          row.names = TRUE)
write.csv(P_adj_matrix,
          file.path(output_dir, "Figure2C_top15_clinical_padj.csv"),
          row.names = TRUE)

# ----------------------------- Heatmap --------------------------------------
rho  <- R_species_clinical
padj <- P_adj_matrix

colnames(rho)[colnames(rho) == "OxygenationIndex"]   <- "Oxygenation_Index"
colnames(padj)[colnames(padj) == "OxygenationIndex"] <- "Oxygenation_Index"

rho_long  <- melt(as.matrix(rho),  varnames = c("Species", "Clinical"),
                  value.name = "rho")
padj_long <- melt(as.matrix(padj), varnames = c("Species", "Clinical"),
                  value.name = "FDR")

df <- left_join(rho_long, padj_long, by = c("Species", "Clinical"))

df <- df %>%
  mutate(sig = case_when(
    FDR < 0.001 ~ "***",
    FDR < 0.01  ~ "**",
    FDR < 0.05  ~ "*",
    TRUE        ~ ""
  ))

species_order <- top15_table[[1]]
species_order <- species_order[species_order != "Others" & !is.na(species_order)]
species_order <- intersect(species_order, rownames(rho))
species_order <- rev(species_order)

df$Species  <- factor(df$Species, levels = species_order)
df$Clinical <- factor(df$Clinical, levels = colnames(rho))

color_low  <- "#3B7A9E"
color_mid  <- "#FFFFFF"
color_high <- "#B84A39"

p <- ggplot(df, aes(x = Clinical, y = Species)) +
  geom_tile(aes(fill = rho), color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = color_low, mid = color_mid, high = color_high,
    midpoint = 0, limits = c(-1, 1),
    name = "Correlation coefficient (rho)"
  ) +
  geom_text(aes(label = sig), color = "black", size = 4.5, vjust = 0.75) +
  scale_x_discrete(
    expand = c(0, 0),
    labels = c("Oxygenation_Index" = expression(paste("PaO"[2], "/", "FiO"[2])))
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_void(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               color = "black", size = 11),
    axis.text.y = element_text(face = "italic", color = "black", size = 10),
    legend.title = element_text(size = 11, face = "plain"),
    legend.text = element_text(size = 10),
    legend.position = "right"
  ) +
  labs(x = NULL, y = NULL)

fig_width  <- 7.5
fig_height <- max(6.0, length(species_order) * 0.25 + 2.0)

ggsave(file.path(output_dir, "Figure2C_top15_clinical_correlation_heatmap.pdf"),
       p, device = cairo_pdf, width = fig_width, height = fig_height, bg = "white")
ggsave(file.path(output_dir, "Figure2C_top15_clinical_correlation_heatmap.png"),
       p, width = fig_width, height = fig_height, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "Figure2C_top15_clinical_correlation_heatmap.tiff"),
       p, width = fig_width, height = fig_height, dpi = 600,
       compression = "lzw", bg = "white")

# ==============================================================================
# Session info
# ==============================================================================
sessionInfo()