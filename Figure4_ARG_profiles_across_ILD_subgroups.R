# ==============================================================================
# Figure 4: Antimicrobial resistance gene (ARG) profiles across ILD subgroups
#
# This script generates:
#   Figure 4A: Total ARG abundance comparison (raw p values)
#   Figure 4B: Top4 drug class abundance (BH-corrected pairwise comparisons)
#   Figure 4C: Top4 resistance mechanism abundance (BH-corrected pairwise comparisons)
#   Figure 4D: Top4 ARG family abundance (BH-corrected pairwise comparisons)
#
# Output: output/Figure4/
# ==============================================================================

# ------------------------------ 0. Packages ----------------------------------
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

# ------------------------------ 1. Paths -------------------------------------
data_dir   <- "data"
output_dir <- "output/Figure4"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ------------------------------ 2. Read data ---------------------------------
gene <- fread(file.path(data_dir, "gene.TPM"),
              header = TRUE, sep = "\t", check.names = FALSE, data.table = FALSE)
rgi  <- fread(file.path(data_dir, "NRprotein_rgi_filtered.txt"),
              header = TRUE, sep = "\t", quote = "", check.names = FALSE, data.table = FALSE)
meta <- read.csv(file.path(data_dir, "metadata.csv"), header = TRUE)

colnames(gene)[1] <- "ORF_ID"
rgi$ORF_clean <- sub(" .*", "", rgi$ORF_ID)

# Unified group levels and palette
group_levels <- c("ILD-only", "ILD with infection", "Infection-related ILD")
group_colors <- c(
  "ILD-only"              = "#0077BB",
  "ILD with infection"    = "#EE7733",
  "Infection-related ILD" = "#CC3311"
)

meta$SampleID <- as.character(meta$SampleID)
meta_sub <- meta %>%
  filter(SubGroup %in% group_levels) %>%
  mutate(SubGroup = factor(SubGroup, levels = group_levels))

# ==============================================================================
# Figure 4A: Total ARG abundance
# ==============================================================================
arg_orfs <- unique(rgi$ORF_clean)
tpm_arg <- gene[gene$ORF_ID %in% arg_orfs, ]
tpm_arg_mat <- tpm_arg[, -1, drop = FALSE]
arg_total_tpm <- colSums(tpm_arg_mat, na.rm = TRUE)

arg_total_df <- data.frame(
  SampleID      = names(arg_total_tpm),
  ARG_total_TPM = as.numeric(arg_total_tpm),
  stringsAsFactors = FALSE
)

plot_df_A <- arg_total_df %>%
  left_join(meta_sub %>% select(SampleID, SubGroup), by = "SampleID") %>%
  filter(!is.na(SubGroup)) %>%
  mutate(SubGroup = factor(SubGroup, levels = group_levels))

# Pairwise Wilcoxon
comb <- combn(group_levels, 2, simplify = FALSE)
res_list <- lapply(comb, function(g) {
  df_sub <- plot_df_A %>% filter(SubGroup %in% g)
  test <- wilcox.test(ARG_total_TPM ~ SubGroup, data = df_sub)
  data.frame(Group1 = g[1], Group2 = g[2],
             W_statistic = test$statistic, p_value = test$p.value)
})
res_df <- do.call(rbind, res_list)
write.table(res_df,
            file.path(output_dir, "Figure4A_total_ARG_pairwise_test.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Descriptive statistics
summary_df <- plot_df_A %>%
  group_by(SubGroup) %>%
  summarise(
    n      = n(),
    mean   = mean(ARG_total_TPM),
    sd     = sd(ARG_total_TPM),
    se     = sd / sqrt(n),
    median = median(ARG_total_TPM),
    Q1     = quantile(ARG_total_TPM, 0.25),
    Q3     = quantile(ARG_total_TPM, 0.75),
    min    = min(ARG_total_TPM),
    max    = max(ARG_total_TPM),
    .groups = "drop"
  )
write.csv(summary_df,
          file.path(output_dir, "Figure4A_total_ARG_descriptive_stats.csv"),
          row.names = FALSE)

plot_df_A$log_TPM <- log10(plot_df_A$ARG_total_TPM + 1)

p_A <- ggplot(plot_df_A, aes(x = SubGroup, y = log_TPM, fill = SubGroup)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85,
               color = "#2B2B2B", linewidth = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.45, color = "#4D4D4D") +
  theme_bw() +
  theme(
    panel.grid      = element_blank(),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.text.x     = element_text(size = 10.5, face = "bold", color = "black"),
    axis.text.y     = element_text(size = 10, color = "black"),
    axis.title      = element_text(size = 12, face = "bold"),
    legend.position = "none",
    axis.ticks      = element_line(color = "black")
  ) +
  scale_fill_manual(values = group_colors) +
  labs(x = NULL, y = expression("Relative abundance (log"[10]("(TPM+1)"))) +
  stat_compare_means(comparisons = comb,
                     method       = "wilcox.test",
                     label        = "p.signif",
                     tip.length   = 0.015,
                     vjust        = 0.2,
                     hide.ns      = FALSE,
                     size         = 4.2)

ggsave(file.path(output_dir, "Figure4A_total_ARG_abundance.pdf"),
       p_A, width = 6, height = 5)
ggsave(file.path(output_dir, "Figure4A_total_ARG_abundance.png"),
       p_A, width = 6, height = 5, dpi = 600)
ggsave(file.path(output_dir, "Figure4A_total_ARG_abundance.tiff"),
       p_A, width = 6, height = 5, dpi = 600, compression = "lzw")

# ==============================================================================
# Helper function for Figure 4B-4D
# ==============================================================================
make_facet_plot <- function(plot_df, category_col, top_n,
                            out_prefix, output_dir,
                            group_colors, group_levels) {
  top_cats <- plot_df %>%
    group_by(.data[[category_col]]) %>%
    summarise(med_TPM = median(Total_TPM, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(med_TPM)) %>%
    slice_head(n = top_n) %>%
    pull(.data[[category_col]])
  
  plot_top <- plot_df %>%
    filter(.data[[category_col]] %in% top_cats) %>%
    mutate(!!category_col := factor(.data[[category_col]], levels = top_cats))
  
  # Pairwise Wilcoxon + BH correction within each category
  pairwise_list <- list()
  for (cat in top_cats) {
    df_cat <- plot_df %>% filter(.data[[category_col]] == cat)
    groups_cat <- unique(as.character(df_cat$SubGroup))
    if (length(groups_cat) >= 2) {
      comb_cat <- combn(groups_cat, 2, simplify = FALSE)
      for (pair in comb_cat) {
        df_pair <- df_cat %>% filter(SubGroup %in% pair)
        test <- wilcox.test(Total_TPM ~ SubGroup, data = df_pair)
        pairwise_list <- append(pairwise_list, list(
          data.frame(
            Category    = cat,
            Group1      = pair[1],
            Group2      = pair[2],
            W_statistic = ifelse(is.null(test$statistic), NA, test$statistic),
            p_value     = test$p.value,
            stringsAsFactors = FALSE
          )
        ))
      }
    }
  }
  pairwise_df <- do.call(rbind, pairwise_list)
  
  if (nrow(pairwise_df) > 0) {
    pairwise_df <- pairwise_df %>%
      group_by(Category) %>%
      mutate(p_adjusted = p.adjust(p_value, method = "BH")) %>%
      ungroup() %>%
      mutate(p_adj_signif = case_when(
        p_adjusted < 0.001 ~ "***",
        p_adjusted < 0.01  ~ "**",
        p_adjusted < 0.05  ~ "*",
        TRUE               ~ "ns"
      ))
  }
  
  write.table(pairwise_df,
              file.path(output_dir, paste0(out_prefix, "_pairwise_BH.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  plot_top$log_TPM <- log10(plot_top$Total_TPM + 1)
  
  facet_stat <- pairwise_df %>%
    mutate(group1 = Group1, group2 = Group2,
           p.adj = p_adjusted, p.adj.signif = p_adj_signif) %>%
    select(Category, group1, group2, p.adj, p.adj.signif)
  
  facet_stat$y.position <- NA
  for (cat in top_cats) {
    idx <- which(facet_stat$Category == cat)
    if (length(idx) == 0) next
    y_max <- max(plot_top$log_TPM[plot_top[[category_col]] == cat], na.rm = TRUE)
    y_min <- min(plot_top$log_TPM[plot_top[[category_col]] == cat], na.rm = TRUE)
    step  <- (y_max - y_min) * 0.10
    facet_stat$y.position[idx] <- y_max + step * seq_along(idx)
  }
  facet_stat$Category <- factor(facet_stat$Category, levels = top_cats)
  
  p <- ggplot(plot_top, aes(x = SubGroup, y = log_TPM, fill = SubGroup)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85,
                 color = "#2B2B2B", linewidth = 0.5, width = 0.6) +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.45, color = "#4D4D4D") +
    facet_wrap(vars(.data[[category_col]]), scales = "free_y", ncol = 2) +
    theme_bw() +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 10.5,
                                      face = "bold", color = "black"),
      axis.text.y      = element_text(size = 10, color = "black"),
      axis.title       = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "#F2F4F4", color = "black", linewidth = 0.5),
      strip.text       = element_text(face = "bold", size = 9, color = "black"),
      legend.position  = "none"
    ) +
    scale_fill_manual(values = group_colors) +
    labs(x = NULL, y = expression("Relative abundance (log"[10]("(TPM+1)"))) +
    stat_pvalue_manual(facet_stat, label = "p.adj.signif",
                       tip.length = 0.015, hide.ns = FALSE, size = 4.2)
  
  ggsave(file.path(output_dir, paste0(out_prefix, ".pdf")),
         p, width = 10, height = 8)
  ggsave(file.path(output_dir, paste0(out_prefix, ".png")),
         p, width = 10, height = 8, dpi = 300)
  ggsave(file.path(output_dir, paste0(out_prefix, ".tiff")),
         p, width = 10, height = 8, dpi = 600, compression = "lzw")
  
  return(p)
}

# ==============================================================================
# Figure 4B: Top4 drug class abundance
# ==============================================================================
drug_map <- rgi %>%
  select(ORF_clean, `Drug Class`) %>%
  filter(!is.na(`Drug Class`) & `Drug Class` != "" & `Drug Class` != "n/a") %>%
  separate_rows(`Drug Class`, sep = ";\\s*") %>%
  rename(ORF_ID = ORF_clean, DRUG_CLASS = `Drug Class`) %>%
  distinct() %>%
  filter(!is.na(DRUG_CLASS))

tpm_arg_B <- gene[gene$ORF_ID %in% unique(drug_map$ORF_ID), ]
tpm_long_B <- tpm_arg_B %>%
  pivot_longer(cols = -ORF_ID, names_to = "SampleID", values_to = "TPM") %>%
  left_join(drug_map, by = "ORF_ID")

drug_abundance <- tpm_long_B %>%
  group_by(SampleID, DRUG_CLASS) %>%
  summarise(Total_TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

plot_df_B <- drug_abundance %>%
  left_join(meta_sub %>% select(SampleID, SubGroup), by = "SampleID") %>%
  filter(!is.na(SubGroup)) %>%
  mutate(SubGroup = factor(SubGroup, levels = group_levels))

p_B <- make_facet_plot(
  plot_df      = plot_df_B,
  category_col = "DRUG_CLASS",
  top_n        = 4,
  out_prefix   = "Figure4B_drug_class_top4",
  output_dir   = output_dir,
  group_colors = group_colors,
  group_levels = group_levels
)

# ==============================================================================
# Figure 4C: Top4 resistance mechanism abundance
# ==============================================================================
mech_map <- rgi %>%
  select(ORF_clean, `Resistance Mechanism`) %>%
  rename(ORF_ID = ORF_clean, MECHANISM = `Resistance Mechanism`) %>%
  distinct() %>%
  filter(!is.na(ORF_ID) & !is.na(MECHANISM))

tpm_arg_C <- gene[gene$ORF_ID %in% unique(mech_map$ORF_ID), ]
tpm_long_C <- tpm_arg_C %>%
  pivot_longer(cols = -ORF_ID, names_to = "SampleID", values_to = "TPM") %>%
  left_join(mech_map, by = "ORF_ID")

mech_abundance <- tpm_long_C %>%
  group_by(SampleID, MECHANISM) %>%
  summarise(Total_TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

plot_df_C <- mech_abundance %>%
  left_join(meta_sub %>% select(SampleID, SubGroup), by = "SampleID") %>%
  filter(!is.na(SubGroup)) %>%
  mutate(SubGroup = factor(SubGroup, levels = group_levels))

p_C <- make_facet_plot(
  plot_df      = plot_df_C,
  category_col = "MECHANISM",
  top_n        = 4,
  out_prefix   = "Figure4C_resistance_mechanism_top4",
  output_dir   = output_dir,
  group_colors = group_colors,
  group_levels = group_levels
)

# ==============================================================================
# Figure 4D: Top4 ARG family abundance
# ==============================================================================
family_map <- rgi %>%
  select(ORF_clean, `AMR Gene Family`) %>%
  rename(ORF_ID = ORF_clean, ARG_FAMILY = `AMR Gene Family`) %>%
  distinct() %>%
  filter(!is.na(ORF_ID) & !is.na(ARG_FAMILY))

tpm_arg_D <- gene[gene$ORF_ID %in% unique(family_map$ORF_ID), ]
tpm_long_D <- tpm_arg_D %>%
  pivot_longer(cols = -ORF_ID, names_to = "SampleID", values_to = "TPM") %>%
  left_join(family_map, by = "ORF_ID")

family_abundance <- tpm_long_D %>%
  group_by(SampleID, ARG_FAMILY) %>%
  summarise(Total_TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

plot_df_D <- family_abundance %>%
  left_join(meta_sub %>% select(SampleID, SubGroup), by = "SampleID") %>%
  filter(!is.na(SubGroup)) %>%
  mutate(SubGroup = factor(SubGroup, levels = group_levels))

p_D <- make_facet_plot(
  plot_df      = plot_df_D,
  category_col = "ARG_FAMILY",
  top_n        = 4,
  out_prefix   = "Figure4D_ARG_family_top4",
  output_dir   = output_dir,
  group_colors = group_colors,
  group_levels = group_levels
)

# ==============================================================================
# Session info
# ==============================================================================
sessionInfo()