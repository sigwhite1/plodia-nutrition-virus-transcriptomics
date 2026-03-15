# R/12_immune_pathway_dynamics.R
# Purpose: Pathway-level statistical analysis of immune gene dynamics across
#          diet x time, using gene-blocked linear models.
#
# ANALYSES:
#   1) Pathway mean logFC timeseries per diet (± SE across genes)
#   2) Diet × Time interaction test per pathway (gene-blocked model)
#   3) Early (≤24h) vs Late (≥96h) phase comparison + Diet × Phase test
#   4) Diet main effect controlling for time + gene identity
#   5) Diet variance differences (Levene test) per pathway
#
# Reads:
#   - intermediate/config.rds
#   - intermediate/immune_plot_df.rds   (from R_11; one row per gene x diet x time)
#     NOTE: immune_plot_df now uses 'category' (GO:BP term) instead of 'Set'
#     (pathway name). All grouping in this script uses 'category'.
#
# Writes:
#   - results/tables/immune_interaction_tests.csv
#   - results/tables/immune_phase_tests.csv
#   - results/tables/immune_diet_main_tests.csv
#   - results/tables/immune_variance_tests.csv
#   - results/figures/immune_pathway_mean_timeseries.png
#   - results/figures/immune_pathway_early_vs_late.png
#   - results/figures/immune_pathway_variance.png

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(purrr)
  library(readr)
  library(car)     # leveneTest
})

# -----------------------------
# 0) Load config + inputs
# -----------------------------
config            <- readRDS(file.path("intermediate", "config.rds"))
paths             <- config$paths
time_levels_hours <- config$time_levels_hours   # numeric: 0.5 1 2 ... 168

dir.create(paths$figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$tables_dir,  showWarnings = FALSE, recursive = TRUE)

immune_plot_df <- readRDS(file.path(paths$intermediate_dir, "immune_plot_df.rds"))

# Validate required columns
req_cols <- c("category", "FBgn", "ProbeName", "diet", "time_h", "logFC", "adj.P.Val")
missing_cols <- setdiff(req_cols, colnames(immune_plot_df))
if (length(missing_cols) > 0) {
  stop("immune_plot_df missing required columns: ",
       paste(missing_cols, collapse = ", "),
       "\nRe-run R_11_immune_gene_evaluation.R first.")
}

# R_11 applies display labels to diet ("Low resource"/"High resource").
# Recode back to raw values ("low"/"high") for modelling and consistent joins.
immune_plot_df <- immune_plot_df %>%
  dplyr::mutate(
    diet = dplyr::case_when(
      as.character(diet) %in% c("low",  "Low resource")  ~ "low",
      as.character(diet) %in% c("high", "High resource") ~ "high",
      TRUE ~ as.character(diet)
    )
  )

# Diet color palette and labels
diet_colors <- c("low" = "#00BFC4", "high" = "#F8766D")
diet_labels <- c("low" = "Low resource", "high" = "High resource")

# -----------------------------
# 1) Collapse probe-level -> gene-level
#    (mean logFC across probes per FBgn x diet x time_h)
# -----------------------------
gene_level_df <- immune_plot_df %>%
  dplyr::group_by(category, FBgn, diet, time_h) %>%
  dplyr::summarise(
    logFC_gene = mean(logFC,     na.rm = TRUE),
    adjP_gene  = min(adj.P.Val,  na.rm = TRUE),
    n_probes   = dplyr::n_distinct(ProbeName),
    .groups    = "drop"
  ) %>%
  dplyr::mutate(
    time_h     = factor(time_h, levels = time_levels_hours, ordered = TRUE),
    time_h_num = as.numeric(as.character(time_h)),  # for use in lm()
    diet       = factor(diet,   levels = c("low", "high")),
    category   = factor(category)
  )

message("gene_level_df: ", nrow(gene_level_df), " rows | ",
        dplyr::n_distinct(gene_level_df$FBgn), " genes | ",
        dplyr::n_distinct(gene_level_df$category), " categories")

# -----------------------------
# 2) Pathway-level mean timeseries
# -----------------------------
pathway_mean <- gene_level_df %>%
  dplyr::group_by(category, diet, time_h) %>%
  dplyr::summarise(
    mean_logFC = mean(logFC_gene, na.rm = TRUE),
    sd_logFC   = sd(logFC_gene,   na.rm = TRUE),
    n_genes    = dplyr::n_distinct(FBgn),
    se_logFC   = sd_logFC / sqrt(pmax(n_genes, 1)),
    .groups    = "drop"
  )

p_pathway <- ggplot2::ggplot(
  pathway_mean,
  ggplot2::aes(x = time_h, y = mean_logFC,
               color = diet, group = diet)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey50") +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean_logFC - se_logFC,
                 ymax = mean_logFC + se_logFC),
    width = 0.15
  ) +
  ggplot2::facet_wrap(~ category, ncol = 3, scales = "free_y") +
  ggplot2::scale_color_manual(values = diet_colors, labels = diet_labels) +
  ggplot2::labs(
    title    = "GO immune category mean response across diet \u00d7 time",
    subtitle = "Mean gene-level logFC within each GO immune category (\u00b1 SE across genes)",
    x        = "Time after infection (hours)",
    y        = "Mean gene logFC",
    color    = "Diet"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    panel.grid       = ggplot2::element_blank(),
    legend.position  = "top",
    strip.background = ggplot2::element_rect(fill = "grey85", color = "grey30"),
    strip.text       = ggplot2::element_text(face = "bold"),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1)
  )

print(p_pathway)
out_ts <- file.path(paths$figures_dir, "immune_pathway_mean_timeseries.png")
ggplot2::ggsave(out_ts, p_pathway, width = 12, height = 8, dpi = 300)
message("Saved: ", out_ts)

# -----------------------------
# 3) Diet × Time interaction test per pathway (gene-blocked model)
#    logFC_gene ~ Diet * Time + FBgn
# -----------------------------
# FBgn is included as a blocking factor to account for the fact that some
# genes consistently have higher or lower logFC regardless of treatment.
# This is equivalent to a repeated-measures design.
# Before fitting, check each pathway has sufficient variation to fit the model
pathway_coverage <- gene_level_df %>%
  dplyr::group_by(category) %>%
  dplyr::summarise(
    n_diets    = dplyr::n_distinct(diet),
    n_times    = dplyr::n_distinct(time_h),
    n_genes    = dplyr::n_distinct(FBgn),
    n_obs      = dplyr::n(),
    .groups    = "drop"
  )
message("Pathway coverage check:")
print(pathway_coverage)

# Flag sets that cannot support the full interaction model
# (need >= 2 diet levels, >= 2 time levels, >= 2 genes)
sets_ok <- pathway_coverage %>%
  dplyr::filter(n_diets >= 2, n_times >= 2, n_genes >= 2) %>%
  dplyr::pull(category)

sets_skipped <- setdiff(pathway_coverage$category, sets_ok)
if (length(sets_skipped) > 0) {
  message("Skipping interaction test for sets with insufficient variation: ",
          paste(sets_skipped, collapse = ", "))
}

interaction_tests <- gene_level_df %>%
  dplyr::filter(category %in% sets_ok) %>%
  dplyr::group_by(category) %>%
  dplyr::group_modify(~ {
    df  <- .x
    # Drop any factor levels not present in this subset
    df  <- dplyr::mutate(df,
                         diet   = droplevels(diet),
                         time_h = droplevels(time_h),
                         FBgn   = droplevels(factor(FBgn)))
    p_int <- tryCatch({
      fit  <- lm(logFC_gene ~ diet * time_h_num + FBgn, data = df)
      a    <- anova(fit)
      term <- grep("diet.*time_h_num|time_h_num.*diet", rownames(a), value = TRUE)[1]
      if (!is.na(term) && term %in% rownames(a)) a[term, "Pr(>F)"] else NA_real_
    }, error = function(e) {
      message("  Model failed for category ", unique(df$category), ": ", conditionMessage(e))
      NA_real_
    })
    tibble::tibble(
      n_genes       = dplyr::n_distinct(df$FBgn),
      n_obs         = nrow(df),
      p_interaction = p_int
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    fdr_interaction = p.adjust(p_interaction, method = "BH")
  ) %>%
  dplyr::arrange(fdr_interaction)

message("\nDiet x Time interaction tests (gene-blocked):")
print(interaction_tests)
readr::write_csv(interaction_tests,
                 file.path(paths$tables_dir, "immune_interaction_tests.csv"))

# -----------------------------
# 4) Early vs Late phase analysis
#    Early: time_h <= 24h  |  Late: time_h >= 96h
# -----------------------------
gene_phase_df <- gene_level_df %>%
  dplyr::mutate(
    Phase = dplyr::case_when(
      as.numeric(as.character(time_h)) <= 24 ~ "Early",
      as.numeric(as.character(time_h)) >= 96 ~ "Late",
      TRUE                                   ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Phase)) %>%
  dplyr::mutate(Phase = factor(Phase, levels = c("Early", "Late")))

pathway_phase_mean <- gene_phase_df %>%
  dplyr::group_by(category, diet, Phase) %>%
  dplyr::summarise(
    mean_logFC = mean(logFC_gene, na.rm = TRUE),
    sd_logFC   = sd(logFC_gene,   na.rm = TRUE),
    n_genes    = dplyr::n_distinct(FBgn),
    se_logFC   = sd_logFC / sqrt(pmax(n_genes, 1)),
    .groups    = "drop"
  )

p_phase <- ggplot2::ggplot(
  pathway_phase_mean,
  ggplot2::aes(x = Phase, y = mean_logFC,
               color = diet, group = diet)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey50") +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean_logFC - se_logFC,
                 ymax = mean_logFC + se_logFC),
    width = 0.15
  ) +
  ggplot2::facet_wrap(~ category, ncol = 3, scales = "free_y") +
  ggplot2::scale_color_manual(values = diet_colors, labels = diet_labels) +
  ggplot2::labs(
    title    = "Early vs Late immune response (GO immune categories)",
    subtitle = "Early \u226424h; Late \u226596h (mean gene-level logFC \u00b1 SE across genes)",
    x        = NULL,
    y        = "Mean gene logFC",
    color    = "Diet"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    panel.grid       = ggplot2::element_blank(),
    legend.position  = "top",
    strip.background = ggplot2::element_rect(fill = "grey85", color = "grey30"),
    strip.text       = ggplot2::element_text(face = "bold")
  )

print(p_phase)
out_phase <- file.path(paths$figures_dir, "immune_pathway_early_vs_late.png")
ggplot2::ggsave(out_phase, p_phase, width = 12, height = 8, dpi = 300)
message("Saved: ", out_phase)

# Diet × Phase interaction test (gene-blocked)
phase_tests <- gene_phase_df %>%
  dplyr::group_by(category) %>%
  dplyr::group_modify(~ {
    df  <- dplyr::mutate(.x,
                         diet  = droplevels(diet),
                         Phase = droplevels(Phase),
                         FBgn  = droplevels(factor(FBgn)))
    p_int <- tryCatch({
      fit <- lm(logFC_gene ~ diet * Phase + FBgn, data = df)
      a   <- anova(fit)
      if ("diet:Phase" %in% rownames(a)) a["diet:Phase", "Pr(>F)"] else NA_real_
    }, error = function(e) NA_real_)
    tibble::tibble(
      n_genes      = dplyr::n_distinct(df$FBgn),
      n_obs        = nrow(df),
      p_diet_phase = p_int
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(fdr_diet_phase = p.adjust(p_diet_phase, method = "BH")) %>%
  dplyr::arrange(fdr_diet_phase)

message("\nDiet x Phase interaction tests:")
print(phase_tests)
readr::write_csv(phase_tests,
                 file.path(paths$tables_dir, "immune_phase_tests.csv"))

# -----------------------------
# 5) Diet main effect (controlling for time + gene identity)
#    logFC_gene ~ Diet + Time + FBgn
# -----------------------------
diet_main_tests <- gene_level_df %>%
  dplyr::group_by(category) %>%
  dplyr::group_modify(~ {
    df  <- dplyr::mutate(.x,
                         diet   = droplevels(diet),
                         time_h = droplevels(time_h),
                         FBgn   = droplevels(factor(FBgn)))
    p_diet <- tryCatch({
      fit <- lm(logFC_gene ~ diet + time_h_num + FBgn, data = df)
      a   <- anova(fit)
      if ("diet" %in% rownames(a)) a["diet", "Pr(>F)"] else NA_real_
    }, error = function(e) NA_real_)
    tibble::tibble(
      n_genes     = dplyr::n_distinct(df$FBgn),
      n_obs       = nrow(df),
      p_diet_main = p_diet
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(fdr_diet_main = p.adjust(p_diet_main, method = "BH")) %>%
  dplyr::arrange(fdr_diet_main)

message("\nDiet main effect tests:")
print(diet_main_tests)
readr::write_csv(diet_main_tests,
                 file.path(paths$tables_dir, "immune_diet_main_tests.csv"))

# -----------------------------
# 6) Diet variance differences (Levene test within each pathway)
# -----------------------------
variance_tests <- gene_level_df %>%
  dplyr::group_by(category) %>%
  dplyr::group_modify(~ {
    df <- .x
    lv <- tryCatch(
      car::leveneTest(logFC_gene ~ diet, data = df),
      error = function(e) NULL
    )
    p_var <- if (!is.null(lv)) lv[1, "Pr(>F)"] else NA_real_
    
    var_tbl <- df %>%
      dplyr::group_by(diet) %>%
      dplyr::summarise(var_logFC = var(logFC_gene, na.rm = TRUE), .groups = "drop")
    
    diets_present <- as.character(var_tbl$diet)
    var_low  <- var_tbl$var_logFC[diets_present == "low"]
    var_high <- var_tbl$var_logFC[diets_present == "high"]
    
    tibble::tibble(
      p_variance = p_var,
      var_low    = if (length(var_low)  == 1) var_low  else NA_real_,
      var_high   = if (length(var_high) == 1) var_high else NA_real_
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(fdr_variance = p.adjust(p_variance, method = "BH")) %>%
  dplyr::arrange(fdr_variance)

message("\nVariance tests (Levene):")
print(variance_tests)
readr::write_csv(variance_tests,
                 file.path(paths$tables_dir, "immune_variance_tests.csv"))

# Variance bar plot
p_var_plot <- gene_level_df %>%
  dplyr::group_by(category, diet) %>%
  dplyr::summarise(var_logFC = var(logFC_gene, na.rm = TRUE), .groups = "drop") %>%
  ggplot2::ggplot(ggplot2::aes(x = diet, y = var_logFC, fill = diet)) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(~ category, scales = "free_y") +
  ggplot2::scale_fill_manual(values = diet_colors, labels = diet_labels) +
  ggplot2::labs(
    title = "Within-category variance of gene-level logFC (GO immune categories)",
    y     = "Variance across genes",
    x     = NULL,
    fill  = "Diet"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 legend.position = "top")

print(p_var_plot)
out_var <- file.path(paths$figures_dir, "immune_pathway_variance.png")
ggplot2::ggsave(out_var, p_var_plot, width = 10, height = 6, dpi = 300)
message("Saved: ", out_var)
