# 02_highsplit_cluster_analysis.R
# Per-cluster comparison of Phase 1 vs BANKSY on the 14 Phase-1 clusters that BANKSY split into ≥10 subclusters. 
# Produces:
#   - all_highsplit_clusters_comparison.csv (per-cluster metrics, both methods)
#   - purity_scatter.png (purity difference vs cell-type count)
#   - pattern_anova_summary.txt (ANOVA test of purity diff across patterns)

source("/scratch/prj/mmg_grp_single_cell_lab/projects/michael/pipeline_output/CRC_Comparative_Analysis/scripts/00_helpers.R")
library(ggplot2)
library(dplyr)

ensure_dirs()

# The 14 high-split Phase 1 clusters, identified by inspecting BANKSY's split counts (clusters where ≥10 BANKSY subclusters had >100 spots within them).
# Pattern annotations follow the architectural classification described in Chapter 4 based on visual inspection of the subcluster spatial layouts.

HIGH_SPLIT_CLUSTERS <- tibble::tribble(
  ~Sample, ~Cluster, ~Pattern,
  "P1",    "3",      "Organised",
  "P1",    "7",      "Organised",
  "P1",    "8",      "Intermixed",
  "P1",    "12",     "Organised",
  "P2",    "0",      "Intermixed",
  "P2",    "2",      "Intermixed",
  "P2",    "4",      "Intermixed",
  "P2",    "6",      "Homogeneous",
  "P2",    "12",     "Intermixed",
  "P5",    "2",      "Intermixed",
  "P5",    "3",      "Homogeneous",
  "P5",    "4",      "Homogeneous",
  "P5",    "6",      "Homogeneous",
  "P5",    "14",     "Homogeneous"
)

# Analysis function — loads sample objects, subclusters, target cluster, and computes per-method metrics. 

analyze_one_cluster <- function(sample_name, cluster_id, pattern,
                                phase1_obj, banksy_obj) {
  message(sprintf("  --- %s cluster %s (%s) ---",
                  sample_name, cluster_id, pattern))
  
  res <- subcluster_and_compose(phase1_obj, banksy_obj, cluster_id,
                                resolution = 0.8, min_banksy_spots = 100,
                                verbose = TRUE)
  
  # Cell-type counts: columns of composition tables with any non-zero entry
  n_ct_standard <- sum(colSums(res$standard_composition > 0) > 0)
  n_ct_banksy   <- sum(colSums(res$banksy_composition  > 0) > 0)
  n_ct_union    <- length(unique(c(colnames(res$standard_composition),
                                   colnames(res$banksy_composition))))
  
  # Chi-square statistics, with safe failure on degenerate tables
  chi_std <- tryCatch(chisq.test(res$standard_composition)$statistic,
                      error = function(e) NA_real_, warning = function(w) {
                        suppressWarnings(chisq.test(res$standard_composition)$statistic)
                      })
  chi_bks <- tryCatch(chisq.test(res$banksy_composition)$statistic,
                      error = function(e) NA_real_, warning = function(w) {
                        suppressWarnings(chisq.test(res$banksy_composition)$statistic)
                      })
  
  std_purity  <- calculate_purity(res$standard_prop)
  bks_purity  <- calculate_purity(res$banksy_prop)
  std_entropy <- calculate_entropy(res$standard_prop)
  bks_entropy <- calculate_entropy(res$banksy_prop)
  
  data.frame(
    Sample                = sample_name,
    Cluster               = cluster_id,
    Pattern               = pattern,
    N_Spots               = length(res$cluster_spots),
    Standard_Subclusters  = res$n_standard,
    BANKSY_Subclusters    = res$n_banksy_significant,
    CellTypes_Standard    = n_ct_standard,
    CellTypes_BANKSY      = n_ct_banksy,
    CellTypes_Union       = n_ct_union,
    Standard_Purity       = round(std_purity,  1),
    BANKSY_Purity         = round(bks_purity,  1),
    Purity_Diff           = round(std_purity - bks_purity, 1),
    Standard_Entropy      = round(std_entropy, 2),
    BANKSY_Entropy        = round(bks_entropy, 2),
    Standard_ChiSq        = round(chi_std, 1),
    BANKSY_ChiSq          = round(chi_bks, 1),
    stringsAsFactors      = FALSE
  )
}

# Batch loop: load sample objects once, run all clusters for that sample.
results <- list()

for (sample_name in SAMPLES) {
  message("\n=== Loading ", sample_name, " ===")
  phase1 <- readRDS(phase1_path(sample_name))
  banksy <- readRDS(banksy_path(sample_name))
  
  sample_clusters <- HIGH_SPLIT_CLUSTERS[HIGH_SPLIT_CLUSTERS$Sample == sample_name, ]
  
  for (i in seq_len(nrow(sample_clusters))) {
    row <- sample_clusters[i, ]
    key <- paste(sample_name, row$Cluster, sep = "_")
    results[[key]] <- analyze_one_cluster(
      sample_name = sample_name,
      cluster_id  = row$Cluster,
      pattern     = row$Pattern,
      phase1_obj  = phase1,
      banksy_obj  = banksy
    )
  }
  
  rm(phase1, banksy); gc()
}

all_results <- do.call(rbind, results)
rownames(all_results) <- NULL

print(all_results)

write.csv(all_results,
          file.path(OUTPUT_DIR, "all_highsplit_clusters_comparison.csv"),
          row.names = FALSE)

# Purity scatter plot
p_scatter <- ggplot(all_results,
                    aes(x = CellTypes_Union, y = Purity_Diff,
                        color = Sample)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_text(aes(label = Cluster),
            hjust = -0.8, vjust = 0.5, size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("P1" = "#E41A1C",
                                "P2" = "#377EB8",
                                "P5" = "#4DAF4A")) +
  labs(
    title    = "Standard vs BANKSY Purity: Patient-Specific Patterns",
    x        = "Number of Cell Types (union of methods)",
    y        = "Purity Difference (Standard − BANKSY)\nPositive = Standard better",
    caption  = "14 high-split clusters (≥10 BANKSY subclusters) across 3 CRC samples"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(OUTPUT_DIR, "purity_scatter.png"),
       p_scatter, width = 10, height = 7, dpi = 120)

# Cross-pattern ANOVA on Purity_Diff
message("Pattern ANOVA")
pattern_means <- all_results %>%
  group_by(Pattern) %>%
  summarise(n      = n(),
            mean   = mean(Purity_Diff),
            sd     = sd(Purity_Diff),
            .groups = "drop")
print(pattern_means)

aov_fit <- aov(Purity_Diff ~ Pattern, data = all_results)
aov_summary <- summary(aov_fit)
print(aov_summary)

# Pearson correlation: does cell-type count predict purity difference?
cor_test <- cor.test(all_results$CellTypes_Union, all_results$Purity_Diff,
                     method = "pearson")
message(sprintf("\nPearson correlation: r = %.2f, p = %.3f",
                cor_test$estimate, cor_test$p.value))

# Write summary text
sink(file.path(OUTPUT_DIR, "pattern_anova_summary.txt"))
cat("Cross-pattern ANOVA on purity difference (Standard − BANKSY)\n")
cat("Per-pattern means (Purity_Diff):\n")
print(as.data.frame(pattern_means))
cat("\n\nANOVA summary:\n")
print(aov_summary)
cat("\nPearson correlation (CellTypes_Union vs Purity_Diff):\n")
cat(sprintf("  r = %.4f, p = %.4f, n = %d\n",
            cor_test$estimate, cor_test$p.value, nrow(all_results)))
sink()

message("Outputs in:", OUTPUT_DIR)
