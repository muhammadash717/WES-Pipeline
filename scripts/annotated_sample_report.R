#!/usr/bin/env Rscript
# =============================================================================
#  Variant Annotation Report Generator
#  Usage:  Rscript variant_report.R <input.tsv>
#  Output: <input>.pdf  (same directory as input)
# =============================================================================

required_pkgs <- c("ggplot2", "dplyr", "tidyr", "gridExtra", "grid", "scales")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  message("Installing: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)
}
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(gridExtra); library(grid); library(scales)
})

# ── 0. Args & paths ───────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript variant_report.R <input.tsv>")
input_file  <- args[1]
if (!file.exists(input_file)) stop("File not found: ", input_file)
sample_name <- tools::file_path_sans_ext(basename(input_file))
output_pdf  <- file.path(dirname(input_file), paste0(sample_name, ".pdf"))

cat("Reading:", input_file, "\n")
df      <- read.delim(input_file, sep = "\t", stringsAsFactors = FALSE,
                      na.strings = c("", "NA", "N/A", "."))
n_total <- nrow(df)
cat("Variants loaded:", n_total, "\n")

# ── 1. Colours & theme ────────────────────────────────────────────────────────
C_BG      <- "#FFFFFF"
C_PANEL   <- "#F8FAFC"
C_GRID    <- "#E2E8F0"
C_BORDER  <- "#CBD5E1"
C_T1      <- "#0F172A"
C_T2      <- "#475569"
C_T3      <- "#94A3B8"
C_ACCENT  <- "#2563EB"
C_ACCENT2 <- "#0891B2"
C_HDR_BG  <- "#1E3A5F"
C_HDR_FG  <- "#FFFFFF"
C_ROW_A   <- "#FFFFFF"
C_ROW_B   <- "#F1F5F9"
C_DIV     <- "#E2E8F0"

PAL_ACMG <- c(
  "Pathogenic"               = "#DC2626",
  "Likely_pathogenic"        = "#EA580C",
  "Likely pathogenic"        = "#EA580C",
  "ModeratePathogenicSupport"= "#D97706",
  "Uncertain_significance"   = "#7C3AED",
  "Uncertain"                = "#7C3AED",
  "LowBenignSupport"         = "#059669",
  "Likely_benign"            = "#10B981",
  "LikelyBenign"             = "#10B981",
  "Benign"                   = "#0284C7",
  "drug response"            = "#6366F1",
  "Other"                    = "#94A3B8"
)
PAL_EFFECT   <- c("#0EA5E9","#F97316","#84CC16","#A855F7","#F43F5E","#14B8A6","#FB923C")
PAL_ZYG      <- c(heterozygous = "#3B82F6", homozygous = "#F43F5E")
PAL_CHROM    <- colorRampPalette(c("#1D4ED8","#60A5FA","#BAE6FD"))(25)
PAL_INSILICO <- c(Benign                = "#0284C7",
                  Uncertain_significance= "#7C3AED",
                  Pathogenic            = "#DC2626",
                  Other                 = "#94A3B8")

base_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace% theme(
    plot.background   = element_rect(fill = C_BG,    colour = NA),
    panel.background  = element_rect(fill = C_PANEL, colour = NA),
    panel.grid.major  = element_line(colour = C_GRID, linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    axis.line         = element_line(colour = C_BORDER, linewidth = 0.4),
    axis.ticks        = element_line(colour = C_BORDER, linewidth = 0.3),
    axis.text         = element_text(colour = C_T2, size = base_size * 0.82),
    axis.title        = element_text(colour = C_T1, size = base_size * 0.9),
    plot.title        = element_text(colour = C_T1, face = "bold",
                                     size = base_size * 1.05, hjust = 0,
                                     margin = margin(0, 0, 4, 0)),
    plot.subtitle     = element_text(colour = C_T2, size = base_size * 0.85,
                                     hjust = 0, margin = margin(0, 0, 8, 0)),
    plot.caption      = element_text(colour = C_T3, size = base_size * 0.72, hjust = 1),
    legend.background = element_rect(fill = C_BG, colour = NA),
    legend.key        = element_rect(fill = NA,   colour = NA),
    legend.text       = element_text(colour = C_T2, size = base_size * 0.80),
    legend.title      = element_text(colour = C_T1, size = base_size * 0.88, face = "bold"),
    plot.margin       = margin(10, 14, 10, 14),
    strip.background  = element_rect(fill = "#EFF6FF", colour = NA),
    strip.text        = element_text(colour = C_ACCENT, face = "bold", size = base_size * 0.88)
  )
}

theme_bar <- function(...) {
  base_theme(...) + theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = C_GRID, linewidth = 0.4)
  )
}

theme_pie <- function(...) {
  base_theme(...) + theme(
    axis.text  = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), axis.line  = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.42, "cm"),
    legend.text     = element_text(size = 8.5, colour = C_T2)
  )
}

# ── 2. Page assembly helpers ──────────────────────────────────────────────────

page_header <- function(title, subtitle = NULL) {
  title_g <- textGrob(title, x = 0, hjust = 0,
    gp = gpar(col = C_T1, fontface = "bold", fontsize = 15))
  rule_g  <- linesGrob(x = c(0, 1), y = c(0, 0),
    gp = gpar(col = C_ACCENT, lwd = 2))
  if (!is.null(subtitle)) {
    sub_g <- textGrob(subtitle, x = 0, hjust = 0,
      gp = gpar(col = C_T2, fontsize = 9, fontface = "italic"))
    arrangeGrob(title_g, rule_g, sub_g,
                heights = unit(c(0.52, 0.10, 0.38), "cm"), ncol = 1)
  } else {
    arrangeGrob(title_g, rule_g, heights = unit(c(0.52, 0.10), "cm"), ncol = 1)
  }
}

light_table <- function(df_tab, fontsize = 9.5) {
  nr  <- nrow(df_tab)
  nc  <- ncol(df_tab)
  
  # All data rows centered
  hjust_core <- rep(0.5, nr * nc)
  x_core     <- rep(0.5, nr * nc)
  
  # All headers centered
  hjust_hdr <- rep(0.5, nc)
  x_hdr     <- rep(0.5, nc)
  
  tt <- ttheme_minimal(
    core = list(
      fg_params = list(col = C_T1, fontsize = fontsize,
                       hjust = hjust_core, x = x_core),
      bg_params = list(fill = rep(c(C_ROW_A, C_ROW_B), length.out = nr),
                       col  = C_GRID)
    ),
    colhead = list(
      fg_params = list(col = C_HDR_FG, fontface = "bold", fontsize = fontsize + 0.5,
                       hjust = hjust_hdr, x = x_hdr),
      bg_params = list(fill = C_HDR_BG, col = C_HDR_BG)
    )
  )
  tableGrob(df_tab, rows = NULL, theme = tt)
}

draw_page <- function(grob, page_num) {
  grid.newpage()
  grid.rect(gp = gpar(fill = C_BG, col = NA))
  vp <- viewport(x = 0.5, y = 0.5, width = 0.94, height = 0.94)
  pushViewport(vp)
  grid.draw(grob)
  popViewport()
  grid.text(paste("Page", page_num), x = 0.97, y = 0.018,
            just = c("right","bottom"),
            gp = gpar(col = C_T3, fontsize = 7.5))
}

# ── 3. Data preparation ───────────────────────────────────────────────────────
ps     <- suppressWarnings(as.numeric(df$priority_score))
ps_min <- round(min(ps, na.rm = TRUE), 1)
ps_max <- round(max(ps, na.rm = TRUE), 1)

n_genes  <- length(unique(na.omit(df$franklin_gene)))
n_chrom  <- length(unique(na.omit(df$chr)))
n_hetero <- sum(df$zygosity == "heterozygous", na.rm = TRUE)
n_homo   <- sum(df$zygosity == "homozygous",   na.rm = TRUE)

summary_tbl <- data.frame(
  Metric = c("Total Variants","Unique Genes","Chromosomes",
             "Heterozygous","Homozygous","PS Range"),
  Value  = c(n_total, n_genes, n_chrom, n_hetero, n_homo,
             paste0(ps_min," \u2013 ",ps_max)),
  stringsAsFactors = FALSE
)

ps_tbl <- data.frame(
  `Score Range` = c("> 40","20 \u2013 40","10 \u2013 20","0 \u2013 10","< 0"),
  Count = c(sum(ps > 40,            na.rm=TRUE),
            sum(ps > 20 & ps<=40,   na.rm=TRUE),
            sum(ps > 10 & ps<=20,   na.rm=TRUE),
            sum(ps >= 0 & ps<=10,   na.rm=TRUE),
            sum(ps <  0,            na.rm=TRUE)),
  check.names = FALSE, stringsAsFactors = FALSE
)

chr_order <- c(as.character(1:22), "X","Y","M")
chr_df <- df %>% filter(!is.na(chr)) %>% count(chr) %>%
  mutate(chr = factor(chr, levels = chr_order)) %>%
  filter(!is.na(chr)) %>% arrange(chr)

effect_df <- df %>%
  filter(!is.na(franklin_effect), franklin_effect != "") %>%
  count(franklin_effect) %>%
  mutate(pct = n/sum(n), franklin_effect = reorder(franklin_effect, n))

zyg_df <- df %>% filter(!is.na(zygosity)) %>% count(zygosity) %>%
  mutate(pct = n/sum(n), pos = cumsum(pct) - pct/2,
         lbl = paste0(round(pct*100,1),"%"))

pred_cols <- c(AlphaMissense   = "genebe_alphamissense_prediction",
               REVEL           = "genebe_revel_prediction",
               SpliceAI        = "genebe_spliceai_max_prediction",
               `BayesDel-noAF` = "genebe_bayesdelnoaf_prediction",
               `PhyloP-100way` = "genebe_phylop100way_prediction")
pred_cols <- pred_cols[pred_cols %in% names(df)]

insilico_long <- do.call(rbind, lapply(names(pred_cols), function(nm) {
  col <- pred_cols[[nm]]
  tmp <- df[!is.na(df[[col]]) & df[[col]] != "", ]
  if (nrow(tmp) == 0) return(NULL)
  tbl <- as.data.frame(table(tmp[[col]]), stringsAsFactors = FALSE)
  colnames(tbl) <- c("prediction","n"); tbl$tool <- nm; tbl$pct <- tbl$n/sum(tbl$n); tbl
}))
insilico_sum <- NULL
if (!is.null(insilico_long) && nrow(insilico_long) > 0) {
  insilico_long$cat <- ifelse(
    grepl("enign",    insilico_long$prediction, ignore.case=TRUE), "Benign",
    ifelse(grepl("athog",   insilico_long$prediction, ignore.case=TRUE), "Pathogenic",
    ifelse(grepl("ncertain",insilico_long$prediction, ignore.case=TRUE),
           "Uncertain_significance", "Other")))
  insilico_sum <- aggregate(n ~ tool + cat, insilico_long, sum)
  insilico_sum <- do.call(rbind, lapply(split(insilico_sum, insilico_sum$tool), function(x) {
    x$pct <- x$n/sum(x$n); x }))
}

clinvar_df <- df %>%
  filter(!is.na(genebe_clinvar_classification), genebe_clinvar_classification != "") %>%
  count(genebe_clinvar_classification) %>%  rename(classification = 1) %>% mutate(pct = n/sum(n)) %>%
  slice_max(n, n = 15) %>% mutate(classification = reorder(classification, n))

top_genes <- df %>% filter(!is.na(franklin_gene)) %>%
  count(franklin_gene, name = "n") %>% arrange(desc(n)) %>%
  slice_head(n = 15) %>% mutate(franklin_gene = reorder(franklin_gene, n))

hpo_df <- df %>%
  mutate(hpo_n = suppressWarnings(as.integer(Matched_HPO_Count))) %>%
  filter(!is.na(hpo_n))
hpo_ratio_df <- df %>%
  mutate(ratio = suppressWarnings(as.numeric(Matched_HPO_Ratio)),
         ps_v  = suppressWarnings(as.numeric(priority_score))) %>%
  filter(!is.na(ratio), !is.na(ps_v))

# --- Additional table: Top priority variants (for quick review) --------------
top_priority_variants <- function(df, n = 20) {
  if (!"priority_score" %in% names(df)) return(NULL)
  top <- df %>%
    filter(!is.na(priority_score)) %>%
    arrange(desc(priority_score)) %>%
    head(n) %>%
    select(
      `Gene` = franklin_gene,
      `Variant` = franklin_c_dot,
      `Zygosity` = zygosity,
      `Franklin` = franklin_classification,
      `GeneBe` = genebe_acmg_classification,
      `OMIM` = OMIM
    ) %>%
    mutate(`OMIM` = gsub(" \\| ", "\n", `OMIM`)) %>%
    mutate(`Zygosity` = gsub("zygous", "", `Zygosity`))
  return(top)
}

# ── 4. Build ggplot objects ───────────────────────────────────────────────────

make_pie <- function(data, col, title) {
  tbl <- data %>%
    filter(!is.na(.data[[col]]), .data[[col]] != "") %>%
    count(.data[[col]], name = "n") %>% rename(label = 1) %>%
    mutate(pct = n/sum(n), label = factor(label, levels = label[order(-n)]),
           pos = cumsum(pct) - pct/2,
           lbl = ifelse(pct >= 0.04, paste0(round(pct*100,1),"%"), ""))
  lvls <- levels(tbl$label)
  cols <- PAL_ACMG[lvls]; cols[is.na(cols)] <- "#94A3B8"; names(cols) <- lvls
  ggplot(tbl, aes("", pct, fill = label)) +
    geom_col(width = 1, colour = "white", linewidth = 0.7) +
    coord_polar("y", start = 0) +
    geom_text(aes(label = lbl), position = position_stack(vjust = 0.5),
              colour = "white", fontface = "bold", size = 3.5) +
    scale_fill_manual(values = cols, name = NULL) +
    labs(title = paste0("\n", title)) + theme_void() +
    theme(legend.position = "bottom")
}

p_franklin_pie <- make_pie(df, "franklin_classification",   "Franklin Classification")
p_genebe_pie   <- make_pie(df, "genebe_acmg_classification","GeneBe Classification")

p_hist <- ggplot(data.frame(ps = ps[!is.na(ps)]), aes(x = ps)) +
  geom_histogram(aes(fill = after_stat(x)), bins = 28,
                 colour = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = C_ACCENT, mid = "#8B5CF6", high = "#DC2626",
                       midpoint = 0, guide = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = C_T2, linewidth = 0.7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Priority Score Distribution",
       subtitle = sprintf("Median %.1f  |  Mean %.1f  |  Range %s \u2013 %s",
                          median(ps,na.rm=TRUE), mean(ps,na.rm=TRUE), ps_min, ps_max),
       x = "Priority Score", y = "Count") +
  theme_bar()

p_chr <- ggplot(chr_df, aes(chr, n, fill = chr)) +
  geom_col(show.legend = FALSE, width = 0.75) +
  geom_text(aes(label = n), vjust = -0.35, colour = C_T2, size = 2.8) +
  scale_fill_manual(values = setNames(PAL_CHROM[seq_along(levels(chr_df$chr))],
                                      levels(chr_df$chr))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Variants per Chromosome", x = "Chromosome", y = "Count") +
  theme_bar() + theme(axis.text.x = element_text(size = 8.5))

p_effect <- ggplot(effect_df, aes(franklin_effect, n, fill = franklin_effect)) +
  geom_col(show.legend = FALSE, width = 0.7) +
  geom_text(aes(label = sprintf("%d  (%.0f%%)", n, pct*100)),
            hjust = -0.08, colour = C_T2, size = 3.1) +
  coord_flip() +
  scale_fill_manual(values = setNames(PAL_EFFECT[seq_len(nrow(effect_df))],
                                      levels(effect_df$franklin_effect))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.32))) +
  labs(title = "Variant Effect Distribution", x = NULL, y = "Count") +
  theme_bar() + theme(panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(colour = C_GRID))

p_zyg <- ggplot(zyg_df, aes("", pct, fill = zygosity)) +
  geom_col(width = 1, colour = "white", linewidth = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = lbl), position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 3.5) +
  scale_fill_manual(values = PAL_ZYG, name = NULL) +
  labs(title = "Zygosity") + theme_void()

if (!is.null(insilico_sum) && nrow(insilico_sum) > 0) {
  p_insilico <- ggplot(insilico_sum, aes(tool, pct, fill = cat)) +
    geom_col(width = 0.68, colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(pct >= 0.07, paste0(round(pct*100),"%"), "")),
              position = position_stack(vjust = 0.5),
              colour = "white", size = 3.2, fontface = "bold") +
    scale_fill_manual(values = PAL_INSILICO, name = "Prediction") +
    scale_y_continuous(labels = percent_format(),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = "In-Silico Prediction Tool Summary",
         subtitle = "Proportion of classified variants per tool",
         x = NULL, y = "Proportion") +
    theme_bar() + theme(axis.text.x = element_text(size = 9.5),
                        legend.position = "right")
} else {
  p_insilico <- ggplot() +
    annotate("text", 0.5, 0.5, label = "No in-silico data available",
             colour = C_T2, size = 5) + theme_bar()
}

p_genes <- ggplot(top_genes, aes(franklin_gene, n, fill = n)) +
  geom_col(show.legend = FALSE, width = 0.7) +
  geom_text(aes(label = n), hjust = -0.2, colour = C_T2, size = 3.2) +
  coord_flip() +
  scale_fill_gradient(low = "#BFDBFE", high = C_ACCENT) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "Top Genes by Variant Count", x = NULL, y = "Variants") +
  theme_bar() + theme(panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(colour = C_GRID))

p_clinvar <- if (nrow(clinvar_df) > 0) {
  cv_cols <- PAL_ACMG[as.character(levels(clinvar_df$classification))]
  cv_cols[is.na(cv_cols)] <- "#94A3B8"; names(cv_cols) <- levels(clinvar_df$classification)
  ggplot(clinvar_df, aes(classification, n, fill = classification)) +
    geom_col(show.legend = FALSE, width = 0.7) +
    geom_text(aes(label = sprintf("%d  (%.0f%%)", n, pct*100)),
              hjust = -0.08, colour = C_T2, size = 3.1) +
    coord_flip() +
    scale_fill_manual(values = cv_cols) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.32))) +
    labs(title = "ClinVar Classification", x = NULL, y = "Count") +
    theme_bar() + theme(panel.grid.major.y = element_blank(),
                        panel.grid.major.x = element_line(colour = C_GRID))
} else NULL

p_hpo <- if (nrow(hpo_df) > 0 && max(hpo_df$hpo_n, na.rm=TRUE) > 0) {
  ggplot(hpo_df, aes(x = hpo_n)) +
    geom_histogram(fill = C_ACCENT2, colour = "white",
                   bins = max(hpo_df$hpo_n)+1, alpha = 0.9) +
    scale_x_continuous(breaks = 0:max(hpo_df$hpo_n)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = "HPO Term Match Count per Variant",
         x = "Matched HPO Terms", y = "Variant Count") +
    theme_bar()
} else NULL

p_scatter <- if (nrow(hpo_ratio_df) > 3) {
  ggplot(hpo_ratio_df, aes(ratio, ps_v)) +
    geom_point(colour = C_ACCENT, alpha = 0.65, size = 2.2) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "#DC2626",
                se = FALSE, linewidth = 1) +
    labs(title = "HPO Ratio vs Priority Score",
         x = "HPO Match Ratio", y = "Priority Score") +
    theme_bar()
} else NULL

# ── 5. Render PDF ─────────────────────────────────────────────────────────────
cat("Writing:", output_pdf, "\n")
cairo_pdf(output_pdf, width = 11.7, height = 8.27, onefile = TRUE)

# ── PAGE 1: Cover ─────────────────────────────────────────────────────────────
grid.newpage()
grid.rect(gp = gpar(fill = C_BG, col = NA))
grid.rect(x = 0, y = 0, width = 0.007, height = 1,
          just = c("left","bottom"), gp = gpar(fill = C_ACCENT, col = NA))
grid.text("VARIANT ANNOTATION REPORT",
          x = 0.06, y = 0.70, just = c("left","bottom"),
          gp = gpar(col = C_T1, fontsize = 30, fontface = "bold"))
grid.text(sample_name, x = 0.06, y = 0.61, just = c("left","bottom"),
          gp = gpar(col = C_ACCENT, fontsize = 18))
grid.lines(x = c(0.06, 0.94), y = c(0.585, 0.585),
           gp = gpar(col = C_DIV, lwd = 1.5))
metrics <- list(
  list(val = n_total,  lbl = "Total Variants"),
  list(val = n_genes,  lbl = "Unique Genes"),
  list(val = n_chrom,  lbl = "Chromosomes"),
  list(val = n_hetero, lbl = "Heterozygous"),
  list(val = n_homo,   lbl = "Homozygous")
)
for (i in seq_along(metrics)) {
  cx <- 0.06 + (i-1)*0.18
  grid.rect(x=cx, y=0.37, width=0.16, height=0.17, just=c("left","bottom"),
            gp=gpar(fill="#EFF6FF", col="#BFDBFE", lwd=1))
  grid.text(as.character(metrics[[i]]$val),
            x=cx+0.08, y=0.455, just=c("centre","bottom"),
            gp=gpar(col=C_ACCENT, fontsize=22, fontface="bold"))
  grid.text(metrics[[i]]$lbl,
            x=cx+0.08, y=0.42, just=c("centre","bottom"),
            gp=gpar(col=C_T2, fontsize=10))
}
grid.text(paste("Generated:", format(Sys.Date(), "%B %d, %Y")),
          x=0.06, y=0.31, just=c("left","bottom"),
          gp=gpar(col=C_T3, fontsize=9.5))

# ── PAGE 2: Summary + Priority Score ─────────────────────────────────────────
hdr2 <- page_header("Summary Overview")

tbl_label_s <- textGrob("Sample Statistics", x=0, hjust=0,
  gp=gpar(col=C_T1, fontface="bold", fontsize=11))
tbl_label_p <- textGrob("Priority Score Tiers", x=0, hjust=0,
  gp=gpar(col=C_T1, fontface="bold", fontsize=11))
spacer <- rectGrob(gp=gpar(fill=NA, col=NA))

left_col <- arrangeGrob(
  tbl_label_s,
  spacer,
  light_table(summary_tbl),
  spacer,
  tbl_label_p,
  spacer,
  light_table(ps_tbl),
  heights = unit(c(0.5, 0.2, 4.6, 1, 0.5, 0.2, 3.8), "cm"), ncol = 1)

body2 <- arrangeGrob(left_col, ggplotGrob(p_hist),
                     ncol=2, widths=c(1, 1.6))
draw_page(arrangeGrob(hdr2, body2,
                      heights=unit(c(0.75, 12.0), "cm"), ncol=1), 2)

# # ── PAGE 3: ACMG Classifications ──────────────────────────────────────────────
# hdr3  <- page_header("ACMG Pathogenicity Classifications",
#   "\nFranklin and GeneBe classification distributions")
# body3 <- arrangeGrob(ggplotGrob(p_franklin_pie),
#                      ggplotGrob(p_genebe_pie), nrow=2)
# draw_page(arrangeGrob(hdr3, spacer, body3,
#                       heights=unit(c(2, 0.2, 16), "cm"), ncol=1), 3)

# # ── PAGE 4: Genomic Region & Effect ───────────────────────────────────────────
# hdr4  <- page_header("Variant Landscape")
# body4 <- arrangeGrob(ggplotGrob(p_region), ggplotGrob(p_effect), ncol=2)
# draw_page(arrangeGrob(hdr4, body4,
#                       heights=unit(c(0.75, 7.0), "cm"), ncol=1), 4)

# # ── PAGE 5: Zygosity & Chromosomes ────────────────────────────────────────────
# hdr5  <- page_header("Zygosity & Chromosomal Distribution")
# body5 <- arrangeGrob(ggplotGrob(p_zyg), ggplotGrob(p_chr),
#                      ncol=2, widths=c(0.7, 1.3))
# draw_page(arrangeGrob(hdr5, body5,
#                       heights=unit(c(0.75, 7.0), "cm"), ncol=1), 5)

# # ── PAGE 6: In-Silico Predictions ─────────────────────────────────────────────
# hdr6  <- page_header("In-Silico Pathogenicity Predictions",
#   "Proportion of classified variants across computational tools")
# body6 <- ggplotGrob(p_insilico)
# draw_page(arrangeGrob(hdr6, body6,
#                       heights=unit(c(0.9, 6.85), "cm"), ncol=1), 6)

# # ── PAGE 7: Top Genes & ClinVar ───────────────────────────────────────────────
# hdr7  <- page_header("Gene & ClinVar Overview")
# body7 <- if (!is.null(p_clinvar)) {
#   arrangeGrob(ggplotGrob(p_genes), ggplotGrob(p_clinvar), ncol=2)
# } else ggplotGrob(p_genes)
# draw_page(arrangeGrob(hdr7, body7,
#                       heights=unit(c(0.75, 7.0), "cm"), ncol=1), 7)

# # ── PAGE 8: HPO (conditional) ─────────────────────────────────────────────────
# if (!is.null(p_hpo)) {
#   hdr8  <- page_header("Phenotype (HPO) Matching Analysis")
#   body8 <- if (!is.null(p_scatter)) {
#     arrangeGrob(ggplotGrob(p_hpo), ggplotGrob(p_scatter), ncol=2)
#   } else ggplotGrob(p_hpo)
#   draw_page(arrangeGrob(hdr8, body8,
#                         heights=unit(c(0.75, 7.0), "cm"), ncol=1), 8)
# }

# ── PAGE 9: Top Priority Variants ─────────────────────────────────────────────────
hdr9 <- page_header("Top Priority Variants")

top_df <- top_priority_variants(df, n = 10)
body9 <- arrangeGrob(light_table(top_df), ncol=1)

draw_page(arrangeGrob(hdr9, body9, heights=unit(c(1, 15), "cm"), ncol=1), 3)

dev.off()
cat("Done! Report saved to:", output_pdf, "\n")