# ==============================================================================
# Figure 5: Virulence factor (VF) profiles across ILD subgroups
#
# This script generates:
#   Figure 5A: Total VF gene abundance comparison
#   Figure 5B-5J: Nine VF functional categories
#     (Adherence, Antimicrobial activity/Competitive advantage, Biofilm,
#      Effector delivery system, Immune modulation, Motility,
#      Nutritional/Metabolic factor, Regulation, Stress survival)
#
# Output: output/Figure5/
# ==============================================================================

# ------------------------------ 0. Packages ----------------------------------
library(tidyverse)
library(data.table)
library(ggpubr)
library(stringr)

# ------------------------------ 1. Paths -------------------------------------
data_dir   <- "data"
output_dir <- "output/Figure5"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

group_levels <- c("ILD-only", "ILD with infection", "Infection-related ILD")
group_colors <- c(
  "ILD-only"              = "#0077BB",
  "ILD with infection"    = "#EE7733",
  "Infection-related ILD" = "#CC3311"
)

# ------------------------------ 2. Read data ---------------------------------
vfdb_file <- file.path(data_dir, "protein_vfdb_setA_identity70_len50.tsv")
tpm_file  <- file.path(data_dir, "gene.TPM")
meta_file <- file.path(data_dir, "metadata.csv")

vfdb <- fread(vfdb_file, header = FALSE, data.table = FALSE)
colnames(vfdb)[1]  <- "gene_id"
colnames(vfdb)[13] <- "description"

tpm <- fread(tpm_file, header = TRUE, data.table = FALSE)
colnames(tpm)[1] <- "gene_id"

meta <- fread(meta_file, data.table = FALSE) %>%
  mutate(SampleID = as.character(SampleID)) %>%
  filter(SubGroup %in% group_levels) %>%
  mutate(SubGroup = factor(SubGroup, levels = group_levels))

# ==============================================================================
# Figure 5A: Total VF gene abundance
# ==============================================================================
vf_orfs <- unique(vfdb$gene_id)

tpm_vf <- tpm[tpm$gene_id %in% vf_orfs, ]
if (nrow(tpm_vf) == 0) stop("No VF ORFs matched TPM matrix.")

tpm_vf_mat <- tpm_vf[, -1, drop = FALSE]
vf_total_tpm <- colSums(tpm_vf_mat, na.rm = TRUE)

vf_total_df <- data.frame(
  SampleID     = names(vf_total_tpm),
  VF_total_TPM = as.numeric(vf_total_tpm),
  stringsAsFactors = FALSE
)

write_tsv(vf_total_df, file.path(output_dir, "Figure5A_VF_total_TPM_by_sample.tsv"))

plot_df_A <- vf_total_df %>%
  left_join(meta %>% select(SampleID, SubGroup), by = "SampleID") %>%
  filter(!is.na(SubGroup)) %>%
  mutate(SubGroup = factor(SubGroup, levels = group_levels))

plot_df_A$log_VF_TPM <- log10(plot_df_A$VF_total_TPM + 1)

# Overall test
test_all <- kruskal.test(VF_total_TPM ~ SubGroup, data = plot_df_A)
writeLines(
  paste0("Kruskal-Wallis chi-squared = ", round(test_all$statistic, 3),
         ", df = ", test_all$parameter,
         ", p = ", signif(test_all$p.value, 4)),
  file.path(output_dir, "Figure5A_VF_total_group_test.txt")
)

# Pairwise comparisons (raw p values shown on the plot, consistent with Figure 5A)
res_pair <- compare_means(VF_total_TPM ~ SubGroup, data = plot_df_A,
                          method = "wilcox.test", p.adjust.method = "BH")
write.csv(res_pair,
          file.path(output_dir, "Figure5A_VF_total_pairwise_comparison.csv"),
          row.names = FALSE)

comparisons <- combn(group_levels, 2, simplify = FALSE)

p_A <- ggplot(plot_df_A, aes(x = SubGroup, y = log_VF_TPM, fill = SubGroup)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA,
               color = "#2B2B2B", linewidth = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.45, color = "#4D4D4D") +
  ggtitle("Total VF gene abundance") +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text.x     = element_text(size = 10.5, face = "bold", color = "black"),
    axis.text.y     = element_text(size = 10, color = "black"),
    axis.title      = element_text(size = 12, face = "bold"),
    axis.ticks      = element_line(color = "black"),
    plot.title      = element_text(size = 11, face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  scale_fill_manual(values = group_colors) +
  labs(x = NULL, y = expression("Relative abundance (log"[10]("(TPM+1)"))) +
  stat_compare_means(comparisons = comparisons,
                     method       = "wilcox.test",
                     label        = "p.signif",
                     hide.ns      = FALSE,
                     size         = 4.2,
                     tip.length   = 0.015)

ggsave(file.path(output_dir, "Figure5A_total_VF_abundance.pdf"),
       p_A, width = 6, height = 5.5)
ggsave(file.path(output_dir, "Figure5A_total_VF_abundance.png"),
       p_A, width = 6, height = 5.5, dpi = 600)
ggsave(file.path(output_dir, "Figure5A_total_VF_abundance.tiff"),
       p_A, width = 6, height = 5.5, dpi = 600, compression = "lzw")

# ==============================================================================
# Figure 5B-5J: Nine VF functional categories
# ==============================================================================

# Valid VF functional categories
valid_vfg <- c(
  "Adherence",
  "Antimicrobial activity/Competitive advantage",
  "Biofilm",
  "Effector delivery system",
  "Exoenzyme",
  "Exotoxin",
  "Immune modulation",
  "Invasion",
  "Motility",
  "Nutritional/Metabolic factor",
  "Post-translational modification",
  "Regulation",
  "Stress survival",
  "Others"
)

# Clean VFDB annotations: extract VFCID and VFG_category
vfdb_clean <- vfdb %>%
  mutate(
    VFCID        = str_extract(description, "VFC\\d{4}"),
    VFG_category = str_extract(description, "(?<=- ).*(?= \\(VFC)")
  ) %>%
  filter(
    !is.na(VFCID),
    !is.na(VFG_category),
    VFG_category %in% valid_vfg
  ) %>%
  distinct(gene_id, VFCID, VFG_category)

# Long TPM table
tpm_long <- tpm %>%
  pivot_longer(cols = -gene_id, names_to = "SampleID", values_to = "TPM") %>%
  mutate(SampleID = as.character(SampleID))

vf_tpm <- inner_join(vfdb_clean, tpm_long, by = "gene_id")

# Aggregate TPM by category and sample
vfg_cat_tpm <- vf_tpm %>%
  group_by(VFG_category, SampleID) %>%
  summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

write_tsv(vfg_cat_tpm, file.path(output_dir, "Figure5B_J_VFG_category_TPM.tsv"))

# Merge group information
vfg_cat_tpm <- vfg_cat_tpm %>%
  left_join(meta %>% select(SampleID, SubGroup), by = "SampleID") %>%
  filter(!is.na(SubGroup)) %>%
  mutate(SubGroup = factor(SubGroup, levels = group_levels))

# Overall Kruskal-Wallis + BH correction
diff_stats <- vfg_cat_tpm %>%
  group_by(VFG_category) %>%
  summarise(
    p_value = kruskal.test(TPM ~ SubGroup)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
    sig = case_when(
      q_value < 0.001 ~ "***",
      q_value < 0.01  ~ "**",
      q_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )
write_tsv(diff_stats, file.path(output_dir, "Figure5B_J_VFG_category_diff_stats.tsv"))

# Pairwise comparisons (Wilcoxon + BH, per category)
pairwise_list <- list()
for (cat in unique(vfg_cat_tpm$VFG_category)) {
  data_cat <- vfg_cat_tpm %>% filter(VFG_category == cat)
  if (n_distinct(data_cat$SubGroup) < 2) next
  res <- compare_means(TPM ~ SubGroup, data = data_cat,
                       method = "wilcox.test", p.adjust.method = "BH") %>%
    mutate(VFG_category = cat)
  pairwise_list[[cat]] <- res
}
pairwise_df <- bind_rows(pairwise_list) %>%
  mutate(p.adj.signif = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*",
    TRUE          ~ "ns"
  ))
write.csv(pairwise_df,
          file.path(output_dir, "Figure5B_J_VFG_category_pairwise_comparison.csv"),
          row.names = FALSE)

# Nine individual category plots
selected_categories <- c(
  "Adherence",
  "Antimicrobial activity/Competitive advantage",
  "Biofilm",
  "Effector delivery system",
  "Immune modulation",
  "Motility",
  "Nutritional/Metabolic factor",
  "Regulation",
  "Stress survival"
)

for (cat in selected_categories) {
  df_cat <- vfg_cat_tpm %>% filter(VFG_category == cat)
  if (nrow(df_cat) == 0 || n_distinct(df_cat$SubGroup) < 2) next
  
  df_cat$log_TPM <- log10(df_cat$TPM + 1)
  
  pw_cat <- pairwise_df %>%
    filter(VFG_category == cat) %>%
    select(group1, group2, p.adj.signif)
  
  max_y <- max(df_cat$log_TPM, na.rm = TRUE)
  min_y <- min(df_cat$log_TPM, na.rm = TRUE)
  step  <- (max_y - min_y) * 0.10
  pw_cat$y.position <- max_y + step * seq_len(nrow(pw_cat))
  
  p_single <- ggplot(df_cat, aes(x = SubGroup, y = log_TPM, fill = SubGroup)) +
    geom_boxplot(alpha = 0.85, outlier.shape = NA,
                 color = "#2B2B2B", linewidth = 0.5, width = 0.6) +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.45, color = "#4D4D4D") +
    ggtitle(cat) +
    theme_bw() +
    theme(
      panel.grid      = element_blank(),
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.text.x     = element_text(size = 10.5, face = "bold", color = "black"),
      axis.text.y     = element_text(size = 10, color = "black"),
      axis.title      = element_text(size = 12, face = "bold"),
      axis.ticks      = element_line(color = "black"),
      plot.title      = element_text(size = 11, face = "bold", hjust = 0.5),
      legend.position = "none"
    ) +
    scale_fill_manual(values = group_colors) +
    labs(x = NULL, y = expression("Relative abundance (log"[10]("(TPM+1)"))) +
    stat_pvalue_manual(pw_cat, label = "p.adj.signif",
                       tip.length = 0.015, hide.ns = FALSE, size = 4.2)
  
  clean_name <- gsub("/", "_", cat)
  clean_name <- gsub(" ", "_", clean_name)
  
  ggsave(file.path(output_dir, paste0("Figure5_", clean_name, ".pdf")),
         p_single, width = 6, height = 5.5)
  ggsave(file.path(output_dir, paste0("Figure5_", clean_name, ".png")),
         p_single, width = 6, height = 5.5, dpi = 600)
  ggsave(file.path(output_dir, paste0("Figure5_", clean_name, ".tiff")),
         p_single, width = 6, height = 5.5, dpi = 600, compression = "lzw")
}

# ==============================================================================
# Session info
# ==============================================================================
sessionInfo()