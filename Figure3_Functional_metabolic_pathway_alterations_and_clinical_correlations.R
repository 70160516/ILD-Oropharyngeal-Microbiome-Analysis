# ==============================================================================
# Figure 3: Functional metabolic pathway alterations and clinical correlations
#
# This script generates:
#   Figure 3A: Bar plot of differential pathways (MaAsLin2, q < 0.2)
#   Figure 3B: Heatmap of differential pathways across samples
#   Figure 3C: Correlation heatmap between differential pathways and clinical indices
#
# Output: output/Figure3/
# ==============================================================================

# ------------------------------ 0. Packages ----------------------------------
library(tidyverse)
library(Maaslin2)
library(pheatmap)
library(psych)
library(RColorBrewer)
library(stringr)

# ------------------------------ 1. Paths -------------------------------------
data_dir   <- "data"
output_dir <- "output/Figure3"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ------------------------------ 2. Read metadata and filter groups -----------
metadata <- read.csv(file.path(data_dir, "metadata.csv"),
                     check.names = FALSE, stringsAsFactors = FALSE)

metadata_sub <- metadata %>%
  filter(SubGroup %in% c("ILD-only", "ILD with infection")) %>%
  mutate(SampleID = as.character(SampleID))

# ------------------------------ 3. Read pathway abundance and clean ----------
raw_abundance <- read.table(
  file.path(data_dir, "merged_pathabundance_clean_relab_unstratified.tsv"),
  header = TRUE, sep = "\t", quote = "", comment.char = "",
  check.names = FALSE, stringsAsFactors = FALSE
)
colnames(raw_abundance)[1] <- "Pathway"

abundance_clean <- raw_abundance %>%
  filter(!str_detect(Pathway, "UNMAPPED|UNINTEGRATED"))

# Standardise sample column names: extract first 4 digits
raw_cols <- colnames(abundance_clean)[-1]
clean_cols <- sub("^([0-9]{4})_.*", "\\1", raw_cols)
colnames(abundance_clean)[-1] <- clean_cols

# ------------------------------ 4. Sample alignment --------------------------
common_samples <- intersect(metadata_sub$SampleID, clean_cols)
if (length(common_samples) == 0) stop("No common samples found.")

metadata_final <- metadata_sub %>%
  filter(SampleID %in% common_samples) %>%
  arrange(match(SampleID, common_samples)) %>%
  column_to_rownames("SampleID")

abundance_final <- abundance_clean %>%
  select(Pathway, all_of(rownames(metadata_final))) %>%
  column_to_rownames("Pathway")

# ------------------------------ 5. Filter pathways ---------------------------
keep_abundance  <- apply(abundance_final, 1, function(x) max(x) >= 0.001)
keep_prevalence <- apply(abundance_final, 1, function(x) sum(x > 0) / length(x) >= 0.10)
abundance_filtered <- abundance_final[keep_abundance & keep_prevalence, ]

# ------------------------------ 6. Run MaAsLin2 ------------------------------
input_data <- as.data.frame(t(abundance_filtered))

maaslin_dir <- file.path(output_dir, "MaAsLin2")
if (!dir.exists(maaslin_dir)) dir.create(maaslin_dir, recursive = TRUE)

fit_data <- Maaslin2(
  input_data       = input_data,
  input_metadata   = metadata_final,
  output           = maaslin_dir,
  fixed_effects    = c("SubGroup"),
  reference        = c("SubGroup,ILD-only"),
  max_significance = 0.2,
  min_abundance    = 0,
  min_prevalence   = 0,
  plot_heatmap     = FALSE,
  plot_scatter     = FALSE
)

# ------------------------------ 7. Restore pathway names ---------------------
name_mapping <- data.frame(
  original_name = rownames(abundance_filtered),
  feature       = make.names(rownames(abundance_filtered)),
  stringsAsFactors = FALSE
)

sig_file <- file.path(maaslin_dir, "significant_results.tsv")
all_file <- file.path(maaslin_dir, "all_results.tsv")

if (file.exists(sig_file) && file.info(sig_file)$size > 0) {
  sig_data <- read.table(sig_file, header = TRUE, sep = "\t", quote = "",
                         check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(sig_data) > 0) {
    sig_data_restored <- sig_data %>%
      mutate(feature = as.character(feature)) %>%
      left_join(name_mapping, by = "feature") %>%
      mutate(feature = ifelse(!is.na(original_name), original_name, feature)) %>%
      select(-original_name)
    write.table(sig_data_restored, sig_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

if (file.exists(all_file) && file.info(all_file)$size > 0) {
  all_data <- read.table(all_file, header = TRUE, sep = "\t", quote = "",
                         check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(all_data) > 0) {
    all_data_restored <- all_data %>%
      mutate(feature = as.character(feature)) %>%
      left_join(name_mapping, by = "feature") %>%
      mutate(feature = ifelse(!is.na(original_name), original_name, feature)) %>%
      select(-original_name)
    write.table(all_data_restored, all_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

# ==============================================================================
# Figure 3A: Bar plot of differential pathways
# ==============================================================================
df_sig <- read.table(sig_file, header = TRUE, sep = "\t", quote = "",
                     check.names = FALSE, stringsAsFactors = FALSE)

if (nrow(df_sig) > 0) {
  df_plot_bar <- df_sig %>%
    arrange(coef) %>%
    mutate(
      feature_clean = factor(feature, levels = unique(feature)),
      Direction = case_when(
        coef > 0 ~ "Enriched in ILD with infection",
        coef < 0 ~ "Enriched in ILD-only"
      )
    )
  
  num_features <- nrow(df_plot_bar)
  ideal_height <- max(4, min(12, num_features * 0.16 + 1.5))
  
  p_bar <- ggplot(df_plot_bar, aes(x = coef, y = feature_clean, fill = Direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_col(width = 0.72, color = NA) +
    scale_fill_manual(
      values = c(
        "Enriched in ILD-only"           = "#0077BB",
        "Enriched in ILD with infection" = "#EE7733"
      )
    ) +
    scale_y_discrete(position = "right") +
    theme_bw() +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.text.y      = element_text(size = 6, color = "black", hjust = 0),
      axis.text.x      = element_text(size = 9.5, color = "black"),
      axis.title.x     = element_text(size = 11, face = "bold"),
      legend.position  = "top",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 9.5),
      legend.background = element_blank()
    ) +
    labs(x = "Effect Size (MaAsLin2 Coefficient)", y = NULL)
  
  ggsave(file.path(output_dir, "Figure3A_differential_pathways_barplot.pdf"),
         p_bar, device = cairo_pdf, width = 9.5, height = ideal_height)
  ggsave(file.path(output_dir, "Figure3A_differential_pathways_barplot.png"),
         p_bar, width = 9.5, height = ideal_height, dpi = 300)
  ggsave(file.path(output_dir, "Figure3A_differential_pathways_barplot.tiff"),
         p_bar, width = 9.5, height = ideal_height, dpi = 600, compression = "lzw")
}

# ==============================================================================
# Figure 3B: Heatmap of differential pathways across samples
# ==============================================================================
sig_pathways <- df_sig$feature

heatmap_mat <- abundance_filtered %>%
  filter(rownames(.) %in% sig_pathways) %>%
  select(all_of(rownames(metadata_final)))

anno_col <- metadata_final %>% select(SubGroup)
anno_colors <- list(SubGroup = c("ILD-only" = "#0077BB",
                                 "ILD with infection" = "#EE7733"))

draw_heatmap <- function() {
  pheatmap(as.matrix(heatmap_mat),
           scale            = "row",
           annotation_col   = anno_col,
           annotation_colors = anno_colors,
           cluster_rows     = TRUE,
           cluster_cols     = TRUE,
           treeheight_row   = 0,
           treeheight_col   = 0,
           show_colnames    = FALSE,
           show_rownames    = TRUE,
           fontsize_row     = 8,
           cellheight       = 16,
           cellwidth        = NA,
           main             = NA,
           color            = colorRampPalette(c("#3B7A9E", "#FFFFFF", "#B84A39"))(100))
}

pdf(file.path(output_dir, "Figure3B_differential_pathways_heatmap.pdf"),
    width = 16, height = 5.5)
draw_heatmap()
dev.off()

png(file.path(output_dir, "Figure3B_differential_pathways_heatmap.png"),
    width = 16, height = 5.5, units = "in", res = 300)
draw_heatmap()
dev.off()

# Save clustering order of pathways (used later for correlation heatmap)
ph_temp <- pheatmap(as.matrix(heatmap_mat),
                    scale            = "row",
                    annotation_col   = anno_col,
                    annotation_colors = anno_colors,
                    cluster_rows     = TRUE,
                    cluster_cols     = TRUE,
                    treeheight_row   = 0,
                    treeheight_col   = 0,
                    show_colnames    = FALSE,
                    show_rownames    = TRUE,
                    fontsize_row     = 8,
                    cellheight       = 16,
                    cellwidth        = NA,
                    main             = NA,
                    color            = colorRampPalette(c("#3B7A9E", "#FFFFFF", "#B84A39"))(100),
                    silent           = TRUE)

clustered_pathways <- rownames(heatmap_mat)[ph_temp$tree_row$order]
writeLines(clustered_pathways,
           file.path(output_dir, "pathway_order_from_clustering.txt"))

# ==============================================================================
# Figure 3C: Correlation between differential pathways and clinical indices
# ==============================================================================
clinical_vars <- c("WBC_count", "Neu_count", "Lym_count",
                   "Eos_count", "CRP", "OxygenationIndex")

metadata_clin <- metadata_sub %>%
  filter(SampleID %in% common_samples) %>%
  select(SampleID, all_of(clinical_vars)) %>%
  column_to_rownames("SampleID")

abundance_corr <- abundance_filtered %>%
  filter(rownames(.) %in% sig_pathways) %>%
  select(all_of(rownames(metadata_clin))) %>%
  t() %>%
  as.data.frame()

metadata_clin <- metadata_clin[rownames(abundance_corr), ]

corr_res <- corr.test(as.matrix(abundance_corr),
                      as.matrix(metadata_clin),
                      method = "spearman", adjust = "none")

r_matrix <- corr_res$r
p_matrix <- corr_res$p
r_matrix[is.na(r_matrix)] <- 0
p_matrix[is.na(p_matrix)] <- 1

q_matrix <- p_matrix
for (j in 1:ncol(p_matrix)) {
  q_matrix[, j] <- p.adjust(p_matrix[, j], method = "BH")
}

write.csv(r_matrix, file.path(output_dir, "Figure3C_pathway_clinical_rho.csv"))
write.csv(p_matrix, file.path(output_dir, "Figure3C_pathway_clinical_p.csv"))
write.csv(q_matrix, file.path(output_dir, "Figure3C_pathway_clinical_padj.csv"))

text_matrix <- matrix("", nrow = nrow(q_matrix), ncol = ncol(q_matrix))
for (i in 1:nrow(q_matrix)) {
  for (j in 1:ncol(q_matrix)) {
    if (q_matrix[i, j] < 0.01) {
      text_matrix[i, j] <- "**"
    } else if (q_matrix[i, j] < 0.05) {
      text_matrix[i, j] <- "*"
    }
  }
}
rownames(text_matrix) <- rownames(r_matrix)
colnames(text_matrix) <- colnames(r_matrix)

colnames(r_matrix) <- gsub("_", " ", colnames(r_matrix))
colnames(r_matrix) <- ifelse(grepl("Oxygenation", colnames(r_matrix), ignore.case = TRUE),
                             "PaO2/FiO2", colnames(r_matrix))
colnames(text_matrix) <- colnames(r_matrix)

order_file <- file.path(output_dir, "pathway_order_from_clustering.txt")
if (file.exists(order_file)) {
  ordered_pathways <- readLines(order_file)
  ordered_pathways <- ordered_pathways[ordered_pathways %in% rownames(r_matrix)]
  if (length(ordered_pathways) == nrow(r_matrix)) {
    r_matrix <- r_matrix[ordered_pathways, , drop = FALSE]
    text_matrix <- text_matrix[ordered_pathways, , drop = FALSE]
  }
}

color_palette <- colorRampPalette(c("#3B7A9E", "#FFFFFF", "#B84A39"))(100)

plot_args <- list(
  mat               = r_matrix,
  color             = color_palette,
  breaks            = seq(-0.5, 0.5, length.out = 101),
  display_numbers   = text_matrix,
  fontsize_number   = 15,
  number_color      = "black",
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  angle_col         = 45,
  cellwidth         = NA,
  cellheight        = NA,
  fontsize_row      = 8,
  fontsize_col      = 10.5,
  border_color      = "white"
)

pdf(file.path(output_dir, "Figure3C_pathway_clinical_correlation_heatmap.pdf"),
    width = 12, height = 8)
do.call(pheatmap, plot_args)
dev.off()

png(file.path(output_dir, "Figure3C_pathway_clinical_correlation_heatmap.png"),
    width = 12, height = 8, units = "in", res = 300)
do.call(pheatmap, plot_args)
dev.off()

# ==============================================================================
# Figure 3D: Genus-level taxonomic contributions to eight core differential
# pathways (ILD-only vs ILD with infection)
#
# Output: output/Figure3/
# ==============================================================================

# ------------------------------ 0. Packages ----------------------------------
library(tidyverse)
library(stringr)

# ------------------------------ 1. Paths -------------------------------------
data_dir   <- "data"
output_dir <- "output/Figure3"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ------------------------------ 2. Parameters --------------------------------
target_groups <- c("ILD-only", "ILD with infection")

target_pathways <- c(
  "ARGSYN-PWY", "GLUTORN-PWY", "HISTSYN-PWY",
  "PWY-6700", "PWY-6703", "RIBOSYN2-PWY",
  "PWY-5973", "PWY-7663"
)

# ------------------------------ 3. Read and clean data -----------------------
meta <- read_csv(file.path(data_dir, "metadata.csv"), show_col_types = FALSE) %>%
  filter(SubGroup %in% target_groups) %>%
  mutate(SampleID = as.character(SampleID)) %>%
  select(SampleID, SubGroup)

abund_raw <- read_tsv(
  file.path(data_dir, "merged_pathabundance_clean_relab.tsv"),
  show_col_types = FALSE
)
colnames(abund_raw)[1] <- "Pathway"

clean_colnames <- colnames(abund_raw) %>%
  map_chr(~ ifelse(.x == "Pathway", "Pathway", substr(.x, 1, 4)))
colnames(abund_raw) <- clean_colnames

# ------------------------------ 4. Extract target pathways ------------------
plot_data_raw <- abund_raw %>%
  filter(str_detect(Pathway, paste0("^(", paste(target_pathways, collapse = "|"), ")"))) %>%
  pivot_longer(cols = -Pathway, names_to = "SampleID", values_to = "Abundance") %>%
  inner_join(meta, by = "SampleID")

# ------------------------------ 5. Parse to genus level ---------------------
plot_data_raw <- plot_data_raw %>%
  separate(Pathway, into = c("Pathway_Name", "Species_Raw"),
           sep = "\\|", fill = "right") %>%
  filter(!is.na(Species_Raw)) %>%
  mutate(Genus = case_when(
    Species_Raw == "unclassified" ~ "unclassified",
    str_detect(Species_Raw, "g__") ~ str_extract(Species_Raw, "(?<=g__)[^.]+"),
    TRUE ~ "unclassified"
  )) %>%
  group_by(SampleID, SubGroup, Pathway_Name, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# ------------------------------ 6. Merge low-abundance genera ---------------
group_mean_data <- plot_data_raw %>%
  group_by(SubGroup, Pathway_Name, Genus) %>%
  summarise(Mean_Abundance = mean(Abundance), .groups = "drop")

main_genera_df <- group_mean_data %>%
  group_by(Pathway_Name, Genus) %>%
  summarise(grand_mean = mean(Mean_Abundance), .groups = "drop")

main_genera <- main_genera_df %>%
  filter(grand_mean >= 0.0001 | Genus == "unclassified") %>%
  pull(Genus) %>% unique()

plot_data <- group_mean_data %>%
  mutate(Genus_Final = ifelse(Genus %in% main_genera, Genus, "Others")) %>%
  group_by(SubGroup, Pathway_Name, Genus_Final) %>%
  summarise(Mean_Abundance = sum(Mean_Abundance), .groups = "drop")

# ------------------------------ 7. Export percentage grid -------------------
percentage_data <- plot_data %>%
  group_by(SubGroup, Pathway_Name) %>%
  mutate(Percentage = (Mean_Abundance / sum(Mean_Abundance)) * 100) %>%
  ungroup()

percentage_grid <- percentage_data %>%
  mutate(Pathway_Short = str_extract(Pathway_Name, "^[^:]+"),
         Group_Genus   = paste0(SubGroup, "_", Genus_Final)) %>%
  select(Pathway_Short, Group_Genus, Percentage) %>%
  pivot_wider(names_from = Group_Genus, values_from = Percentage, values_fill = 0)

write_csv(percentage_grid,
          file.path(output_dir, "Figure3D_core_8_pathways_percentage.csv"))

# ------------------------------ 8. Fix order and wrap names -----------------
genus_rank <- plot_data %>%
  group_by(Genus_Final) %>%
  summarise(total_abundance = sum(Mean_Abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  pull(Genus_Final)

special_cats <- c("unclassified", "Others")
clean_rank <- genus_rank[!genus_rank %in% special_cats]
final_genus_order <- c(special_cats, rev(clean_rank))

plot_data <- plot_data %>%
  mutate(Genus_Final = factor(Genus_Final, levels = final_genus_order)) %>%
  mutate(Pathway_Wrapped = stringr::str_wrap(Pathway_Name, width = 38))

wrapped_levels <- map_chr(target_pathways, function(p) {
  full_name <- unique(plot_data$Pathway_Wrapped[
    str_detect(plot_data$Pathway_Name, paste0("^", p, "\\b"))
  ])
  return(full_name[1])
})

plot_data <- plot_data %>%
  filter(!is.na(Pathway_Wrapped)) %>%
  mutate(Pathway_Wrapped = factor(Pathway_Wrapped, levels = na.omit(wrapped_levels)))

# ------------------------------ 9. Plot -------------------------------------
genus_count <- length(unique(plot_data$Genus_Final))

p_absolute <- ggplot(plot_data, aes(x = SubGroup, y = Mean_Abundance, fill = Genus_Final)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.52) +
  facet_wrap(~ Pathway_Wrapped, ncol = 4, scales = "free_y") +
  theme_bw() +
  labs(x = NULL, y = "Mean Relative Abundance", fill = "Microbial Genus")

p_percent <- ggplot(plot_data, aes(x = SubGroup, y = Mean_Abundance, fill = Genus_Final)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), width = 0.52) +
  facet_wrap(~ Pathway_Wrapped, ncol = 4, scales = "fixed") +
  theme_bw() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Taxonomic Contribution Percentage", fill = "Microbial Genus")

theme_unified <- theme(
  panel.grid       = element_blank(),
  panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
  axis.text.x      = element_text(size = 8, face = "bold", color = "black"),
  axis.text.y      = element_text(size = 7, color = "black"),
  axis.ticks       = element_line(color = "black"),
  strip.text       = element_text(size = 6.5, face = "bold", color = "black", lineheight = 1.1),
  strip.background = element_rect(fill = "#F2F4F4", color = "black", linewidth = 0.5),
  panel.spacing.x  = unit(1.0, "lines"),
  panel.spacing.y  = unit(1.2, "lines"),
  legend.position  = "right",
  legend.title     = element_text(size = 8, face = "bold"),
  legend.text      = element_text(size = 7, face = "italic")
)

p_absolute <- p_absolute + theme_unified
p_percent  <- p_percent  + theme_unified

sci_colors <- c(
  "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666",
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#A6CEE3",
  "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6"
)

if (genus_count <= length(sci_colors)) {
  p_absolute <- p_absolute + scale_fill_manual(
    values = sci_colors[1:genus_count],
    breaks = rev(levels(plot_data$Genus_Final))
  )
  p_percent <- p_percent + scale_fill_manual(
    values = sci_colors[1:genus_count],
    breaks = rev(levels(plot_data$Genus_Final))
  )
} else {
  p_absolute <- p_absolute + scale_fill_viridis_d(option = "turbo",
                                                  breaks = rev(levels(plot_data$Genus_Final)))
  p_percent  <- p_percent  + scale_fill_viridis_d(option = "turbo",
                                                  breaks = rev(levels(plot_data$Genus_Final)))
}

# ------------------------------ 10. Export ----------------------------------
ggsave(file.path(output_dir, "Figure3D_core_8_pathways_absolute.pdf"),
       p_absolute, width = 14.5, height = 7.5, limitsize = FALSE)
ggsave(file.path(output_dir, "Figure3D_core_8_pathways_absolute.png"),
       p_absolute, width = 14.5, height = 7.5, dpi = 300, limitsize = FALSE)

ggsave(file.path(output_dir, "Figure3D_core_8_pathways_normalized.pdf"),
       p_percent, width = 14.5, height = 7.5, limitsize = FALSE)
ggsave(file.path(output_dir, "Figure3D_core_8_pathways_normalized.png"),
       p_percent, width = 14.5, height = 7.5, dpi = 300, limitsize = FALSE)

# ==============================================================================
# Session info
# ==============================================================================
sessionInfo()