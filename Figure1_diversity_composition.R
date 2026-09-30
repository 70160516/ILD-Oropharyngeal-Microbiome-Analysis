# ==============================================================================
# Figure 1: Oropharyngeal microbial diversity and community structure
# across ILD subgroups
# ==============================================================================

# ------------------------------ 0. Packages ----------------------------------
library(tidyverse)
library(ggpubr)
library(cowplot)
library(rstatix)
library(vegan)
library(scales)

# ------------------------------ 1. Paths -------------------------------------
data_dir   <- "data"
output_dir <- "output/Figure1"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ------------------------------ 2. Group levels and palette ------------------
group_levels <- c(
  "ILD-only",
  "ILD with infection",
  "Infection-related ILD"
)

global_group_palette <- c(
  "ILD-only"              = "#0077BB",
  "ILD with infection"    = "#EE7733",
  "Infection-related ILD" = "#CC3311"
)

# Helper: save plot in PDF, PNG, TIFF
save_plot <- function(plot, filename, width, height) {
  ggsave(paste0(filename, ".pdf"), plot, device = cairo_pdf, width = width, height = height)
  ggsave(paste0(filename, ".png"), plot, width = width, height = height, dpi = 600, bg = "white")
  ggsave(paste0(filename, ".tiff"), plot, width = width, height = height, dpi = 600,
         compression = "lzw", bg = "white")
}

# ==============================================================================
# Figure 1A-C: Alpha diversity
# ==============================================================================
alpha_df <- read_csv(file.path(data_dir, "alpha_diversity_indices.csv"),
                     show_col_types = FALSE)
alpha_df$SubGroup <- factor(alpha_df$SubGroup, levels = group_levels)

metrics <- c("Shannon", "Simpson", "Observed")
adjusted_results <- list()

for (m in metrics) {
  test_res <- alpha_df %>%
    wilcox_test(as.formula(paste(m, "~ SubGroup"))) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj",
                     cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", "ns"))
  test_res$Metric <- m
  
  max_val <- max(alpha_df[[m]], na.rm = TRUE)
  val_range <- max_val - min(alpha_df[[m]], na.rm = TRUE)
  test_res <- test_res %>%
    mutate(y.position = c(max_val + 0.05 * val_range,
                          max_val + 0.18 * val_range,
                          max_val + 0.31 * val_range))
  adjusted_results[[m]] <- test_res
}

final_table <- bind_rows(adjusted_results)
write_csv(final_table, file.path(output_dir, "alpha_diversity_statistics_BH.csv"))

# Group descriptive statistics
group_stats <- alpha_df %>%
  pivot_longer(cols = all_of(metrics), names_to = "Metric", values_to = "Value") %>%
  group_by(SubGroup, Metric) %>%
  summarise(Mean = mean(Value, na.rm = TRUE),
            SD = sd(Value, na.rm = TRUE),
            Min = min(Value, na.rm = TRUE),
            Max = max(Value, na.rm = TRUE),
            N = n(), .groups = "drop") %>%
  mutate(Mean_SD = sprintf("%.3f ± %.3f", Mean, SD))

write_csv(group_stats, file.path(output_dir, "alpha_diversity_group_stats.csv"))

# Sample-level data
samples_output <- alpha_df %>%
  select(SampleID, SubGroup, Shannon, Simpson, Observed)
write_csv(samples_output, file.path(output_dir, "alpha_diversity_samples_with_groups.csv"))

# Plot function for alpha diversity
plot_alpha_custom <- function(df, stat_df, metric) {
  max_val <- max(df[[metric]], na.rm = TRUE)
  min_val <- min(df[[metric]], na.rm = TRUE)
  val_range <- max_val - min_val
  
  df_plot <- df %>%
    mutate(SubGroup_Clean = factor(SubGroup, levels = group_levels))
  
  stat_plot <- stat_df %>%
    filter(Metric == metric)
  
  ggplot(df_plot, aes(x = SubGroup_Clean, y = .data[[metric]], fill = SubGroup_Clean)) +
    geom_boxplot(alpha = 0.85, outlier.shape = NA, color = "#2B2B2B", linewidth = 0.5) +
    geom_jitter(position = position_jitter(0.14), size = 1.2, alpha = 0.45,
                color = "#4D4D4D") +
    stat_pvalue_manual(stat_plot, label = "p.adj.signif",
                       xmin = "group1", xmax = "group2",
                       y.position = "y.position", tip_length = 0.015,
                       size = 4.2, vjust = 0.2) +
    labs(x = NULL, y = metric) +
    ylim(min_val - 0.05 * val_range, max_val + 0.45 * val_range) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color = "black", linewidth = 0.6),
          axis.text.x = element_text(size = 10.5, color = "black", face = "bold"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.y = element_text(size = 12, face = "bold"),
          legend.position = "none") +
    scale_fill_manual(values = global_group_palette)
}

p1 <- plot_alpha_custom(alpha_df, final_table, "Shannon")
p2 <- plot_alpha_custom(alpha_df, final_table, "Simpson")
p3 <- plot_alpha_custom(alpha_df, final_table, "Observed")

combined_alpha <- plot_grid(p1, p2, p3, ncol = 3, align = "h", axis = "tb")
save_plot(combined_alpha, file.path(output_dir, "Figure1A_C_alpha_diversity"), 11.5, 4.5)

# ==============================================================================
# Figure 1D: Beta diversity PCoA
# ==============================================================================
rel_df <- read_csv(file.path(data_dir, "relative_abundance_matrix.csv"),
                   show_col_types = FALSE) %>%
  mutate(SampleID = as.character(SampleID))

disp_df <- read_csv(file.path(data_dir, "beta_dispersion_distances.csv"),
                    show_col_types = FALSE) %>%
  mutate(SampleID = as.character(SampleID))

global_p <- read_csv(file.path(data_dir, "beta_permanova_global.csv"),
                     show_col_types = FALSE)

write_csv(global_p, file.path(output_dir, "PERMANOVA_global_results.csv"))

disp_df <- disp_df %>%
  mutate(SubGroup_clean = factor(SubGroup, levels = group_levels))

rel_mat <- rel_df %>% column_to_rownames("SampleID") %>% as.matrix()
bc_dist <- vegdist(rel_mat, method = "bray")
pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE)

pcoa_points <- as.data.frame(pcoa_res$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_eig <- pcoa_res$eig
var_exp <- (pcoa_eig / sum(pcoa_eig[pcoa_eig > 0])) * 100

pcoa_df <- pcoa_points %>%
  rownames_to_column("SampleID") %>%
  left_join(disp_df %>% select(SampleID, SubGroup_clean), by = "SampleID") %>%
  mutate(SubGroup_clean = factor(SubGroup_clean, levels = group_levels))

write_csv(pcoa_df, file.path(output_dir, "PCoA_coordinates_with_groups.csv"))

raw_r2 <- global_p$R2[1]
raw_p  <- global_p$P_value[1]
dynamic_stat_label <- sprintf("PERMANOVA\nR² = %.4f\nP = %.3f", raw_r2, raw_p)

p_pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2,
                              color = SubGroup_clean, fill = SubGroup_clean)) +
  stat_ellipse(geom = "polygon", alpha = 0.05, level = 0.95,
               linewidth = 0.5, linetype = "dashed") +
  geom_point(size = 3.8, alpha = 0.85) +
  labs(x = sprintf("PCoA1 (%.2f%%)", var_exp[1]),
       y = sprintf("PCoA2 (%.2f%%)", var_exp[2]),
       color = "Clinical Group", fill = "Clinical Group") +
  annotate("text", x = -Inf, y = Inf, label = dynamic_stat_label,
           hjust = 0, vjust = 1.4, size = 4.2, color = "#2B2B2B") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 0.6),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9),
        plot.margin = unit(c(0.6, 0.6, 0.6, 0.6), "cm")) +
  scale_color_manual(values = global_group_palette) +
  scale_fill_manual(values = global_group_palette)

save_plot(p_pcoa, file.path(output_dir, "Figure1D_beta_diversity_PCoA"), 7.5, 6.0)

# ==============================================================================
# Figure 1E and 1F: Taxonomic composition bar plots
# ==============================================================================
plot_taxa_barplot <- function(counts_file, meta_file, taxon_level = c("Genus", "Species"),
                              top_n = 15, output_prefix) {
  taxon_level <- match.arg(taxon_level)
  counts <- read_tsv(counts_file, show_col_types = FALSE)
  meta <- read_csv(meta_file, show_col_types = FALSE) %>%
    mutate(SampleID = as.character(SampleID), SubGroup = str_trim(SubGroup)) %>%
    filter(SubGroup %in% group_levels)
  
  sample_ids <- intersect(colnames(counts), meta$SampleID)
  if (length(sample_ids) == 0) stop("Sample IDs do not match between counts and metadata.")
  
  counts_long <- counts %>%
    select(name, all_of(sample_ids)) %>%
    pivot_longer(-name, names_to = "SampleID", values_to = "count") %>%
    mutate(SampleID = as.character(SampleID)) %>%
    left_join(meta %>% select(SampleID, SubGroup), by = "SampleID") %>%
    group_by(SampleID) %>%
    mutate(RelAbund = count / sum(count, na.rm = TRUE)) %>%
    ungroup()
  
  group_mean <- counts_long %>%
    group_by(SubGroup, name) %>%
    summarise(MeanRelAbund = mean(RelAbund, na.rm = TRUE), .groups = "drop")
  
  top_taxa_info <- counts_long %>%
    group_by(name) %>%
    summarise(TotalMean = mean(RelAbund, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(TotalMean)) %>%
    slice_head(n = top_n)
  
  top_taxa <- top_taxa_info %>% pull(name)
  
  rel_agg <- group_mean %>%
    mutate(Taxon = if_else(name %in% top_taxa, name, "Others")) %>%
    group_by(SubGroup, Taxon) %>%
    summarise(MeanRelAbund = sum(MeanRelAbund), .groups = "drop")
  
  global_mean <- rel_agg %>%
    group_by(Taxon) %>%
    summarise(Global_Mean = mean(MeanRelAbund), .groups = "drop")
  
  group_mean_plot <- rel_agg %>%
    left_join(global_mean, by = "Taxon") %>%
    mutate(Global_Mean = if_else(Taxon == "Others", -Inf, Global_Mean)) %>%
    arrange(desc(Global_Mean)) %>%
    mutate(Taxon = factor(Taxon, levels = unique(Taxon))) %>%
    select(-Global_Mean) %>%
    mutate(SubGroup_Clean = factor(SubGroup, levels = group_levels))
  
  classic_order <- levels(group_mean_plot$Taxon)
  
  sci_palette <- c(
    "#CAFFBF", "#FFC6FF", "#BDB2FF", "#F67280", "#C06C84",
    "#355C7D", "#6C5B7B", "#94363E", "#BC5353", "#D87A67",
    "#E4A185", "#A7DBCE", "#84C3BE", "#5A9EAA", "#3A6B88",
    "Others" = "#E5E5E5"
  )
  names(sci_palette)[1:top_n] <- classic_order[classic_order != "Others"]
  
  p <- ggplot(group_mean_plot, aes(x = SubGroup_Clean, y = MeanRelAbund, fill = Taxon)) +
    geom_bar(stat = "identity", position = "fill", width = 0.55,
             color = "#4D4D4D", linewidth = 0.3) +
    scale_fill_manual(values = sci_palette, breaks = classic_order) +
    scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
    labs(x = NULL, y = "Relative Abundance",
         fill = paste(taxon_level, "Composition")) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(color = "black", face = "bold", size = 11, lineheight = 0.9),
          axis.text.y = element_text(color = "black", size = 11),
          axis.title = element_text(color = "black", face = "bold", size = 13),
          axis.line = element_line(color = "black", linewidth = 0.6),
          legend.title = element_text(color = "black", face = "bold", size = 14),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.5, "cm")) +
    guides(fill = guide_legend(byrow = TRUE,
                               label.theme = element_text(face = "italic", size = 12),
                               override.aes = list(color = NA)))
  
  save_plot(p, file.path(output_dir, output_prefix), 7.5, 6.5)
  
  # Export tables
  full_group_mean <- group_mean %>%
    pivot_wider(names_from = SubGroup, values_from = MeanRelAbund, values_fill = 0) %>%
    arrange(desc(`ILD-only` + `ILD with infection` + `Infection-related ILD`))
  
  write_tsv(full_group_mean,
            file.path(output_dir, paste0(output_prefix, "_group_mean_relabund_full.tsv")))
  write_tsv(group_mean_plot,
            file.path(output_dir, paste0(output_prefix, "_plot_data_mean_top15.tsv")))
  
  return(p)
}

p_genus <- plot_taxa_barplot(
  counts_file = file.path(data_dir, "genus_counts_filtered.tsv"),
  meta_file   = file.path(data_dir, "metadata.csv"),
  taxon_level = "Genus",
  top_n = 15,
  output_prefix = "Figure1E_genus_barplot_top15"
)

p_species <- plot_taxa_barplot(
  counts_file = file.path(data_dir, "species_counts_filtered.tsv"),
  meta_file   = file.path(data_dir, "metadata.csv"),
  taxon_level = "Species",
  top_n = 15,
  output_prefix = "Figure1F_species_barplot_top15"
)

# ==============================================================================
# Session info
# ==============================================================================
sessionInfo()