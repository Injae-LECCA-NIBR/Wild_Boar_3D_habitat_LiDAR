## ===========================================================================
## Behavior-specific 3D habitat associations of wild boar (Sus scrofa)
## in an urban-edge forest, Bukhansan National Park, Seoul, Republic of Korea
##
## Data: thermal-drone behavioral observations georeferenced by ray-terrain
##   intersection, paired with LiDAR-derived 3D habitat covariates.
## Behaviors: traveling, resting_day, resting_night.
##
## Pipeline: (1) covariate screening (Spearman + iterative VIF); (2) track-level
##   Kruskal-Wallis + Dunn post-hoc; (3) per-covariate univariate logistic odds
##   ratios (OR per +1 SD) for the three pairwise behavioral contrasts;
##   (4) forest-plot and descriptive boxplot figures.
##
## Unit of analysis: the TRACK (a continuous thermal-observation sequence).
##   Points within a track are not independent (over-sampling of movement bouts),
##   so covariates are aggregated to track means and tracks are analyzed as the
##   independent unit. A track-ID random-intercept GLMM is not identifiable here
##   because most tracks are behavior-pure (each belongs to a single state),
##   confounding the random intercept with the response; Section 5 demonstrates
##   this empirically.
##
## Terminology (manuscript convention): results describe habitat USE / ASSOCIATION,
##   not selection or preference. OR > 1 = covariate takes higher values at the
##   focal-behavior tracks than at the reference-behavior tracks.
## ===========================================================================


# ---- 0. Packages -----------------------------------------------------------
# Packages are NOT installed automatically. If any are missing, the script
# stops and prints the install.packages() call to run once.
required <- c("tidyverse", "rlang", "purrr", "broom", "car",
              "rstatix", "ggpubr", "corrplot", "patchwork", "lme4")
missing <- required[!required %in% rownames(installed.packages())]
if (length(missing))
  stop("Missing packages: ", paste(missing, collapse = ", "),
       "\n  install.packages(c(",
       paste0('"', missing, '"', collapse = ", "), "))")
invisible(lapply(required, library, character.only = TRUE))


# ---- 1. Load data ----------------------------------------------------------
# Set the working directory to the repository root (or open the .Rproj) so the
# relative path resolves to repo/data/WB_environment.csv
df_raw <- read.csv("G:/Bukhansan_2025/##writing/Revision/wildboar_environment_dataset.csv", stringsAsFactors = FALSE)

# Output directory for tables and figures (created if absent).
out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)


cat("Columns:\n");  print(names(df_raw))
cat("\nRaw behavior counts:\n"); print(table(df_raw$behavior))
stopifnot("trackID" %in% names(df_raw))


# ---- 2. Preprocessing ------------------------------------------------------
target_behaviors <- c("traveling", "resting_day", "resting_night")

df <- df_raw %>%
  dplyr::mutate(
    behavior = trimws(as.character(behavior)),
    aspect_clean = dplyr::na_if(aspect, -1),
    northness    = cos(aspect_clean * pi / 180),
    eastness     = sin(aspect_clean * pi / 180),
    log_Deck  = log1p(as.numeric(deck)),
    log_Trail = log1p(as.numeric(trail))
  ) %>%
  dplyr::filter(behavior %in% target_behaviors) %>%
  dplyr::mutate(behavior = factor(behavior, levels = target_behaviors))

cat("\nBehavior counts after filter:\n"); print(table(df$behavior))


# ---- 3. Covariate screening: types, NoData, Spearman, one-pass VIF ---------
all_candidate_vars <- c("DEM", "slope", "TRI", "TWI", "TPI",
                        "northness", "eastness",
                        "CHM", "canopy_cover", "gap_fraction", "LAI",
                        "density0", "log_Deck", "log_Trail")
all_candidate_vars <- intersect(all_candidate_vars, names(df))

for (v in all_candidate_vars) {
  if (!is.numeric(df[[v]])) suppressWarnings(df[[v]] <- as.numeric(df[[v]]))
  n_bad <- sum(df[[v]] < -9000, na.rm = TRUE)
  if (n_bad > 0) { df[[v]][df[[v]] < -9000] <- NA }
}

# Spearman correlation among ALL 14 candidates (screening basis)
df_cor  <- df %>% dplyr::select(dplyr::all_of(all_candidate_vars)) %>% na.omit()
cor_mat <- cor(df_cor, method = "spearman")
corrplot::corrplot(cor_mat, method = "color", type = "upper", diag = FALSE,
                   addCoef.col = "black", number.cex = 0.55,
                   tl.col = "black", tl.srt = 45)
title("Spearman correlation (14 candidate covariates)", line = 1)

# Drop two flagged variables before VIF:
## canopy_cover + gap_fraction = 1 (structural redundancy) -> drop gap_fraction
## slope vs TRI high collinearity (r = 0.98)               -> retain TRI, drop slope
candidate_vars <- setdiff(all_candidate_vars, c("gap_fraction", "slope"))

# One-pass iterative VIF reduction (drop highest VIF until all < threshold).
# VIF depends only on predictor structure; the dummy response is just a vehicle.

reduce_vif <- function(data, vars, threshold = 5) {
  vars <- intersect(vars, names(data))
  repeat {
    d <- data %>% dplyr::select(dplyr::all_of(vars)) %>% na.omit() %>%
      scale() %>% as.data.frame()
    if (ncol(d) < 2) break
    set.seed(1); d$.dummy <- rnorm(nrow(d))
    v <- car::vif(lm(.dummy ~ ., data = d))
    
    cat("\nCurrent VIF values:\n")
    print(round(sort(v, decreasing = TRUE), 2))
    
    if (max(v) < threshold) {
      cat(sprintf("\nAll VIF < %g. Stop.\n", threshold))
      break
    }
    drop_v <- names(which.max(v))
    cat(sprintf("Drop '%s' (VIF = %.2f)\n", drop_v, max(v)))
    vars <- setdiff(vars, drop_v)
  }
  vars
}

cat("\nOne-pass VIF reduction (threshold = 5):\n")
final_vars <- reduce_vif(df, candidate_vars, threshold = 5)
cat("\nFinal covariates:\n"); print(final_vars)

vif_check <- function(data, vars) {
  d <- data %>% dplyr::select(dplyr::all_of(vars)) %>% na.omit() %>%
    scale() %>% as.data.frame()
  set.seed(1); d$.dummy <- rnorm(nrow(d))
  sort(car::vif(lm(.dummy ~ ., data = d)), decreasing = TRUE)
}
cat("\nFinal VIF values:\n"); print(round(vif_check(df, final_vars), 2))

# ---- 4. Labels (full lookup; works for whichever covariates survive VIF) ---
var_labels <- c(
  DEM = "Elevation", TRI = "TRI", TWI = "TWI", TPI = "TPI",
  northness = "Northness", eastness = "Eastness",
  CHM = "Canopy Height", canopy_cover = "Canopy Cover", LAI = "Leaf Area Index",
  density0 = "Understory Density",
  log_Deck = "Distance to Deck Roads", log_Trail = "Distance to Hiking Trails"
)
label_order <- unname(var_labels[final_vars])


# ---- 5. Track structure + behavior-purity diagnostic -----------------------
# (a) points and tracks per behavior
sample_structure <- df %>%
  dplyr::group_by(behavior) %>%
  dplyr::summarise(points = dplyr::n(),
                   tracks = dplyr::n_distinct(trackID), .groups = "drop")
cat("\nSample structure:\n"); print(sample_structure)

# (b) behavior purity of tracks
# A track is behavior-pure when it contains a single behavioral state
# (n_behaviors == 1). These pure tracks are what make a track-ID random
# intercept non-identifiable.
purity <- df %>%
  dplyr::group_by(trackID) %>%
  dplyr::summarise(n_behaviors = dplyr::n_distinct(behavior), .groups = "drop") %>%
  dplyr::count(n_behaviors, name = "n_tracks")
cat("\nTracks by number of behaviors contained (1 = behavior-pure):\n")
print(purity)

# (c) Why a track-ID random-intercept GLMM is not used.
# A mixed model with (1 | trackID) would, in principle, account for repeated
# sampling within tracks. It is non-identifiable here: most tracks are
# behavior-pure (see 5b), so trackID is confounded with the response and
# separates it almost perfectly. This is a property of the random-intercept
# structure itself, not of any covariate, so an intercept-only fixed part is
# enough to expose it -- and is in fact cleaner, since fixed covariates would
# otherwise absorb part of the between-track separation. Fitting the model
# confirms this empirically: the random-intercept SD diverges on the logit scale.
# NOTE: glmer() below is expected to emit a convergence/singularity warning;
# that warning is the intended illustration of non-identifiability, not an error.

glmm_dat <- df %>%
  dplyr::filter(behavior %in% c("resting_day", "traveling")) %>%
  dplyr::mutate(y = as.integer(behavior == "resting_day"))

glmm_fit <- lme4::glmer(
  y ~ 1 + (1 | trackID),
  data = glmm_dat, family = binomial,
  control = lme4::glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 2e5)))

# Random-intercept SD on the logit scale.
re_sd <- as.data.frame(lme4::VarCorr(glmm_fit))$sdcor[1]
cat(sprintf("\nTrack-ID random-intercept SD (logit scale) = %.1f\n", re_sd))
cat("Identifiable logistic intercepts are O(1); an SD of this magnitude means the\n",
    "random intercept spans many orders of magnitude in odds, i.e. trackID\n",
    "separates the response almost perfectly (non-identifiable).\n")


# ---- 6. Track-level aggregation (track = independent unit) ------------------
# A behavior-pure track yields one row; a (rare) mixed track contributes one
# row per behavioral state it contains.
track_df <- df %>%
  dplyr::group_by(trackID, behavior) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(final_vars),
                                 ~ mean(.x, na.rm = TRUE)),
                   n_points = dplyr::n(), .groups = "drop")
cat(sprintf("\nTrack-level rows: %d (from %d unique tracks)\n",
            nrow(track_df), dplyr::n_distinct(df$trackID)))


# ---- 7. Table 1 — descriptive (track-level) + Kruskal-Wallis + Dunn ---------
# Reset any open sink() (e.g. left by an IDE debugger) before printing below.
while (sink.number() > 0) sink()

long_track <- track_df %>%
  tidyr::pivot_longer(dplyr::all_of(final_vars),
                      names_to = "variable", values_to = "value") %>%
  tidyr::drop_na(value)

table1_meansd <- long_track %>%
  dplyr::group_by(variable, behavior) %>%
  dplyr::summarise(mean = mean(value), sd = stats::sd(value), .groups = "drop") %>%
  dplyr::mutate(mean_sd = sprintf("%.2f \u00b1 %.2f", mean, sd)) %>%
  dplyr::select(variable, behavior, mean_sd) %>%
  tidyr::pivot_wider(names_from = behavior, values_from = mean_sd)

kw_tbl <- purrr::map_dfr(final_vars, function(v) {
  d  <- track_df %>% dplyr::select(behavior, val = dplyr::all_of(v)) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(behavior = droplevels(behavior))   # drop now-empty levels
  if (dplyr::n_distinct(d$behavior) < 2) return(NULL) # need >= 2 groups
  kw <- stats::kruskal.test(val ~ behavior, data = d)
  tibble::tibble(variable = v, KW_chisq = unname(kw$statistic),
                 KW_df = unname(kw$parameter), KW_p = kw$p.value)
}) %>% dplyr::mutate(KW_p_BH = stats::p.adjust(KW_p, method = "BH"))

dunn_tbl <- purrr::map_dfr(final_vars, function(v) {
  d  <- track_df %>% dplyr::select(behavior, val = dplyr::all_of(v)) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(behavior = droplevels(behavior))    # drop now-empty levels
  if (dplyr::n_distinct(d$behavior) < 2) return(NULL)  # need >= 2 groups
  # rstatix::dunn_test ranks directly (no sink), so it is debugger-safe.
  rstatix::dunn_test(d, val ~ behavior, p.adjust.method = "BH") %>%
    dplyr::mutate(variable = v, .before = 1)
})

cat("\nTable 1 (track-level mean +/- SD):\n"); print(table1_meansd)
cat("\nKruskal-Wallis (track-level):\n");      print(kw_tbl)
readr::write_csv(table1_meansd, file.path(out_dir, "Table1_mean_sd.csv"))
readr::write_csv(kw_tbl,        file.path(out_dir, "Table1_kruskal_wallis.csv"))
readr::write_csv(dunn_tbl,      file.path(out_dir, "Table1_dunn_posthoc.csv"))


# ---- 8. Table 2 — track-level univariate OR (+1 SD) + BH -------------------
# Each row is one track -> independent -> ordinary logistic SE is valid.
compute_or_track <- function(data, focal_beh, ref_beh, vars = final_vars) {
  pair <- data %>%
    dplyr::filter(behavior %in% c(focal_beh, ref_beh)) %>%
    dplyr::mutate(y = as.integer(behavior == focal_beh))
  
  purrr::map_dfr(vars, function(v) {
    tmp <- pair %>% dplyr::select(y, x = dplyr::all_of(v)) %>% tidyr::drop_na()
    if (nrow(tmp) < 8 || stats::sd(tmp$x) == 0) return(NULL)
    tmp$z <- as.numeric(scale(tmp$x))                       # +1 SD (track-level)
    fit <- stats::glm(y ~ z, data = tmp, family = stats::binomial())
    ci  <- tryCatch(suppressMessages(stats::confint(fit)),  # profile-likelihood CI
                    error = function(e) stats::confint.default(fit))
    tibble::tibble(
      comparison = sprintf("%s vs %s", focal_beh, ref_beh),
      variable   = v,
      n_focal = sum(tmp$y == 1), n_ref = sum(tmp$y == 0),
      OR     = exp(stats::coef(fit)["z"]),
      OR_low = exp(ci["z", 1]), OR_high = exp(ci["z", 2]),
      p_value = broom::tidy(fit)$p.value[2]
    )
  }) %>%
    dplyr::mutate(
      p_BH  = stats::p.adjust(p_value, method = "BH"),
      stars = dplyr::case_when(p_BH < 0.001 ~ "***", p_BH < 0.01 ~ "**",
                               p_BH < 0.05  ~ "*",  TRUE ~ "")
    )
}

table2 <- dplyr::bind_rows(
  compute_or_track(track_df, "resting_day",   "traveling"),
  compute_or_track(track_df, "resting_night", "traveling"),
  compute_or_track(track_df, "resting_day",   "resting_night")
)
cat("\nTable 2 (track-level OR per +1 SD, BH-adjusted):\n")
print(table2, n = Inf)
readr::write_csv(table2, file.path(out_dir, "Table2_OR_track_level.csv"))

# ===========================================================================
# 9. Figure — three-panel forest plot (OR per 1 SD)
#    OR > 1 region shaded with the focal-behavior color; OR < 1 with reference.
# ===========================================================================
col_traveling     <- "#66BB6A"   # green
col_resting_day   <- "#42A5F5"   # light blue
col_resting_night <- "#3949AB"   # indigo

relabel_or <- function(or_data) {
  or_data %>%
    dplyr::mutate(variable_f = factor(dplyr::recode(variable, !!!var_labels),
                                      levels = label_order))
}

make_or_panel <- function(or_data, bg_grp1, bg_grp2, panel_title,
                          show_y = TRUE) {
  xr   <- range(c(or_data$OR_low, or_data$OR_high), na.rm = TRUE)
  x_lo <- xr[1] * 0.78
  x_hi <- xr[2] * 1.38
  
  p <- ggplot2::ggplot(or_data, ggplot2::aes(x = OR, y = variable_f)) +
    ggplot2::annotate("rect", xmin = x_lo, xmax = 1, ymin = -Inf, ymax = Inf,
                      fill = bg_grp2, alpha = 0.10) +   # OR < 1 = reference assoc.
    ggplot2::annotate("rect", xmin = 1, xmax = x_hi, ymin = -Inf, ymax = Inf,
                      fill = bg_grp1, alpha = 0.10) +   # OR > 1 = focal assoc.
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed",
                        linewidth = 0.45, color = "grey35") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = OR_low, xmax = OR_high),
                            height = 0.22, linewidth = 0.65, color = "grey20") +
    ggplot2::geom_point(size = 3.0, color = "grey10", shape = 16) +
    ggplot2::geom_text(ggplot2::aes(label = stars, x = OR_high),
                       hjust = -0.25, vjust = 0.5, size = 4.5,
                       color = "grey15", show.legend = FALSE) +
    ggplot2::scale_x_log10(limits = c(x_lo, x_hi),
                           expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::scale_y_discrete(limits = rev(label_order)) +
    ggplot2::labs(x = "Odds ratio (log scale)", y = NULL, title = panel_title) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(size = 11.5, face = "bold",
                                                 hjust = 0.5, color = "grey15",
                                                 margin = ggplot2::margin(b = 6)),
      axis.text.x        = ggplot2::element_text(size = 10),
      axis.title.x       = ggplot2::element_text(size = 11),
      panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.4),
      panel.grid.major.x = ggplot2::element_line(color = "grey88", linewidth = 0.4),
      panel.grid.minor   = ggplot2::element_blank(),
      plot.margin        = ggplot2::margin(8, 14, 5, 5)
    )
  
  if (show_y) {
    p + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10.5,
                                                           lineheight = 0.85))
  } else {
    p + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())
  }
}

or_a <- relabel_or(dplyr::filter(table2, comparison == "resting_day vs traveling"))
or_b <- relabel_or(dplyr::filter(table2, comparison == "resting_night vs traveling"))
or_c <- relabel_or(dplyr::filter(table2, comparison == "resting_day vs resting_night"))

pa <- make_or_panel(or_a, bg_grp1 = col_resting_day,   bg_grp2 = col_traveling,
                    panel_title = "(a) Traveling vs Resting (Daytime)",   show_y = TRUE)
pb <- make_or_panel(or_b, bg_grp1 = col_resting_night, bg_grp2 = col_traveling,
                    panel_title = "(b) Traveling vs Resting (Nighttime)", show_y = FALSE)
pc <- make_or_panel(or_c, bg_grp1 = col_resting_day,   bg_grp2 = col_resting_night,
                    panel_title = "(c) Resting (Nighttime) vs Resting (Daytime)", show_y = FALSE)

p_or <- (pa | pb | pc) +
  patchwork::plot_layout(ncol = 3, widths = c(1, 1, 1))

print(p_or)
ggplot2::ggsave(file.path(out_dir, "Fig_OR_pairwise_1x3.png"),
                p_or, width = 12, height = 6, dpi = 300)

## ===========================================================================
## Track-level descriptives + boxplots for the variables that differed
## significantly among behavioral states (track-level KW, BH < 0.05):
## distance to deck roads, understory density, elevation, eastness, northness.
## Unit = track, matching the KW/OR analysis. Distances/elevation in metres;
## aspect as northness/eastness in [-1, 1].
## Expects `df` and palette objects col_traveling / col_resting_day /
## col_resting_night defined upstream.
## ===========================================================================

# (packages already loaded in Section 0)

# Panel order: deck and understory density (management headlines) first.
sig_vars <- c("deck", "density0", "DEM", "eastness", "northness")
stopifnot(all(sig_vars %in% names(df)))

# ---- Track-level aggregation (raw, manager-interpretable units) -------------
desc_track <- df %>%
  dplyr::mutate(deck = dplyr::if_else(deck < 0, NA_real_, as.numeric(deck)),
                DEM  = dplyr::if_else(DEM  < -9000, NA_real_, as.numeric(DEM))) %>%
  dplyr::group_by(trackID, behavior) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(sig_vars), ~ mean(.x, na.rm = TRUE)),
                   .groups = "drop") %>%
  dplyr::mutate(behavior = factor(behavior,
                                  levels = c("traveling", "resting_day", "resting_night")))

# ---- Summary table (report median [Q1-Q3] in the manuscript text) -----------
summary_tbl <- desc_track %>%
  tidyr::pivot_longer(dplyr::all_of(sig_vars), names_to = "variable", values_to = "value") %>%
  tidyr::drop_na(value) %>%
  dplyr::group_by(variable, behavior) %>%
  dplyr::summarise(
    n_tracks = dplyr::n(),
    median   = stats::median(value),
    Q1       = stats::quantile(value, 0.25),
    Q3       = stats::quantile(value, 0.75),
    mean     = mean(value),
    sd       = stats::sd(value),
    min      = min(value),
    max      = max(value),
    .groups  = "drop") %>%
  dplyr::mutate(variable = factor(variable, levels = sig_vars)) %>%
  dplyr::arrange(variable, behavior) %>%
  dplyr::mutate(dplyr::across(c(median, Q1, Q3, mean, sd, min, max), ~ round(.x, 2)))

print(summary_tbl, n = Inf)
readr::write_csv(summary_tbl, file.path(out_dir, "Table_sig_descriptives_tracklevel.csv"))

# ---- Boxplots: one panel per variable --------------------------------------
n_beh    <- table(desc_track$behavior)
beh_cols <- c(traveling = col_traveling,
              resting_day = col_resting_day,
              resting_night = col_resting_night)
beh_labs <- c(traveling     = sprintf("Traveling\n(n=%d)",      n_beh["traveling"]),
              resting_day   = sprintf("Resting (Day)\n(n=%d)",   n_beh["resting_day"]),
              resting_night = sprintf("Resting (Night)\n(n=%d)", n_beh["resting_night"]))
var_titles <- c(deck = "Distance to deck roads (m)", density0 = "Understory density",
                DEM = "Elevation (m)", eastness = "Eastness", northness = "Northness")

# Overlay BH-adjusted Dunn pairwise brackets (matches Table 1 post-hoc).
use_brackets <- requireNamespace("rstatix", quietly = TRUE) &&
  requireNamespace("ggpubr",  quietly = TRUE)

make_box <- function(vv) {
  d <- desc_track %>%
    dplyr::select(behavior, value = dplyr::all_of(vv)) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(behavior = droplevels(behavior))
  
  p <- ggplot2::ggplot(d, ggplot2::aes(behavior, value, fill = behavior)) +
    ggplot2::geom_boxplot(width = 0.55, alpha = 0.55, outlier.shape = NA,
                          linewidth = 0.4, colour = "grey25") +
    ggplot2::geom_jitter(width = 0.12, height = 0, size = 1.5,
                         alpha = 0.5, colour = "grey15") +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 2.4,
                          fill = "white", colour = "grey15") +   # mean = diamond
    ggplot2::scale_fill_manual(values = beh_cols) +
    ggplot2::scale_x_discrete(labels = beh_labs) +
    ggplot2::labs(x = NULL, y = var_titles[[vv]], title = var_titles[[vv]]) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "none",
                   plot.title  = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
                   axis.text.x = ggplot2::element_text(size = 10),
                   panel.grid.major.x = ggplot2::element_blank())
  
  if (use_brackets && dplyr::n_distinct(d$behavior) >= 2) {
    st <- rstatix::dunn_test(d, value ~ behavior, p.adjust.method = "BH") %>%
      rstatix::add_xy_position(x = "behavior")
    p <- p + ggpubr::stat_pvalue_manual(st, label = "p.adj.signif",
                                        tip.length = 0.01, hide.ns = TRUE, size = 4)
  }
  p
}

p_box <- patchwork::wrap_plots(lapply(sig_vars, make_box), nrow = 2)

print(p_box)
ggplot2::ggsave(file.path(out_dir, "Fig_sig_boxplots.png"),
                p_box, width = 12, height = 7, dpi = 600)

# ---- Session info (record R and package versions for reproducibility) ------
sessionInfo()
