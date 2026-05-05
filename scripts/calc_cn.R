#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Script: calc_cn_fast.R
# Description: High-performance Copy Number estimation using bamsignals.
# Usage: Rscript calc_cn_fast.R --bam sample.bam --output "file prefix"
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install(c("Rsamtools", "GenomicRanges", "rtracklayer", "bamsignals", "DNACopy"))
# install.packages(c("optparse", "ggplot2", "dplyr", "tidyr", "data.table"))
# ------------------------------------------------------------------------------

# --- 0. Argument Parsing ---

suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-b", "--bam"), type = "character", default = NULL,
              help = "Path to input BAM file (Required). Must be sorted and indexed.", metavar = "file.bam"),
  make_option(c("-r", "--bed"), type = "character", default = NULL,
              help = "Path to BED file defining regions (Optional).", metavar = "regions.bed"),
  make_option(c("-k", "--interval"), type = "integer", default = 50,
              help = "Window size in kbp (default: 50). Ignored if BED is provided.", metavar = "50"),
  make_option(c("-o", "--output"), type = "character", default = "cn_report",
              help = "Output file prefix.", metavar = "output")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$bam) || !file.exists(opt$bam)) {
  print_help(opt_parser)
  stop("Error: Valid BAM file required.", call. = FALSE)
}
if (!file.exists(paste0(opt$bam, ".bai")) && !file.exists(sub(".bam$", ".bai", opt$bam))) {
  stop("Error: BAM index (.bai) not found.", call. = FALSE)
}

# --- 1. Package Loading ---

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(dplyr)
  library(Rsamtools)

  if (!require("bamsignals", quietly = TRUE)) {
    stop("Package 'bamsignals' is required for fast counting. Install via BiocManager::install('bamsignals').")
  }
})

# --- 2. Define Intervals (Bins) ---

message("Reading BAM file...")
bf <- BamFile(opt$bam)
seq_info <- seqinfo(bf)

# Filter to chromosomes (1-22, X, Y)
standard_chr <- grep("^chr[0-9XY]+$|^[0-9XY]+$", seqlevels(seq_info), value = TRUE)
if(length(standard_chr) > 0) seq_info <- seq_info[standard_chr]

if (!is.null(opt$bed)) {
  message(paste("Using regions from BED:", basename(opt$bed)))
  regions_gr <- import(opt$bed)
  regions_gr <- keepSeqlevels(regions_gr, seqlevels(seq_info), pruning.mode = "coarse")
} else {
  message(paste("Tiling genome (", opt$interval, "kbp windows)..."))
  regions_gr <- tileGenome(seq_info, tilewidth = opt$interval * 1000, cut.last.tile.in.chrom = TRUE)
}
regions_gr <- regions_gr[width(regions_gr) > 0]

# --- 3. Counting with bamsignals ---

message("Counting reads...")

counts_vec <- as.numeric(bamsignals::bamCount(
    opt$bam, 
    regions_gr, 
    mapq = 59, # filters out poor quality reads
    verbose = FALSE, # keeps the console clean
    ss = FALSE # Set TRUE if strand-specific counting is needed (rare for CNV)
))

# --- 4. Process and Calculate Copy Number ---

message("Calculating Copy Number...")

regions_df <- as.data.frame(regions_gr)
regions_df$counts <- counts_vec

# Calculate RPKM / Normalized Depth
total_reads <- sum(regions_df$counts)
regions_df$length_kb <- regions_df$width / 1000

# CPM (Counts Per Million mapped reads per kb)
regions_df$cpm <- (regions_df$counts / regions_df$length_kb) / (total_reads / 1e6)

# Estimate Copy Number
baseline_mean <- mean(regions_df$cpm[regions_df$cpm > 0]) 
regions_df$copy_number <- (regions_df$cpm / baseline_mean) * 2
regions_df$midpoint_Mb <- (regions_df$start + regions_df$width/2) / 1e6

# --- 5. Aggregation ---

chrom_summary <- regions_df %>%
  group_by(seqnames) %>%
  summarize(mean_cn = mean(copy_number, na.rm = TRUE))

write.table(x = regions_df, file = paste0(opt$output, "_data.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(x = chrom_summary, file = paste0(opt$output, "_chrom_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# --- 6. Visualization ---

message("Generating plots...")

# Plot 1: Detailed Scatter
p1 <- ggplot(regions_df, aes(x = midpoint_Mb, y = copy_number)) +
  geom_point(aes(color = copy_number > 5), alpha = 0.2, size = 0.1) + 
  scale_color_manual(values = c("FALSE" = "lightgray", "TRUE" = "yellow4"), guide = "none") +
  geom_hline(yintercept = 1.3, color = "red", linetype = "dashed", alpha = 0.9, linewidth = 0.3) +
  geom_hline(yintercept = 2.7, color = "red", linetype = "dashed", alpha = 0.9, linewidth = 0.3) +
  # geom_smooth(aes(x = midpoint_Mb, y = ifelse(copy_number < 0.3, NA, copy_number)),
  #         method = "lm", formula = y ~ 1, se = FALSE, color = "green", linewidth = 0.3) +
  geom_smooth(aes(x = midpoint_Mb, y = ifelse(copy_number < 0.3, NA, copy_number)),
          method = "loess", span = 1, se = FALSE, color = "blue", linewidth = 0.3) + # smaller span = more sensitive (< 0.3 = very wiggly, 0.3–0.5 = more sensitive, 0.75 = default (very smooth))
  facet_wrap(~seqnames, scales = "free_x") +
  theme_bw() +
  labs(title = paste("Copy Number Profile:", basename(opt$bam)), x = "Position (Mbp)", y = "Estimated Copy Number") +
  theme(strip.background = element_rect(fill="#ecf0f1")) +
  coord_cartesian(ylim = c(0, 6))
ggsave(paste0(opt$output, "_cn_profile.png"), plot = p1, dpi = 300)

# Plot 2: Chromosome Bar Chart
p2 <- ggplot(chrom_summary, aes(x = seqnames, y = mean_cn)) +
  geom_col(fill = "#3498db", color = "black") +
  geom_text(aes(label = round(mean_cn, 2)), vjust = -0.5, size = 2) +
  geom_hline(yintercept = 1.5, color = "red", linetype = "dashed", alpha = 0.9, linewidth = 0.3) +
  geom_hline(yintercept = 2.5, color = "red", linetype = "dashed", alpha = 0.9, linewidth = 0.3) +
  theme_bw() +
  labs(title = paste("Mean Copy Number per Chromosome: ", basename(opt$bam)), x = "Chromosome", y = "Estimated Copy Number") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(0, 4))
ggsave(paste0(opt$output, "_chromosome.png"), plot = p2, dpi = 300)

message("Done! Plots & Data tables saved with prefix: ", opt$output)

# # Plot 0: Detailed Scatter
# param <- ScanBamParam(what = "mapq")
# bam_data <- scanBam(opt$bam, param = param)

# mapq <- bam_data[[1]]$mapq
# mapq <- mapq[!is.na(mapq)]

# p0 <- ggplot(data.frame(mapq = mapq), aes(x = mapq)) +
#   geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
#   labs(
#     title = "Mapping Quality Distribution",
#     x = "MAPQ",
#     y = "Read Count"
#   ) +
#   theme_minimal()
# ggsave(paste0(opt$output, "_bam_qual.png"), plot = p0, dpi = 300)

