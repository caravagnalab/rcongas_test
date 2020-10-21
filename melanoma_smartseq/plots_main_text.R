require(Rcongas)
require(tidyverse)
library(cowplot)

# Final model fit
fit = Rcongas::GBM_smartseq_normal_vs_tumor
input_raw_counts_genes = Rcongas::get_input_raw_data(fit)

# DE analysis -- DESeq2 NON funziona
# fit = Rcongas::calculate_DE(
#   x = fit,
#   input = input_raw_counts_genes,
#   clone1 = 1,
#   clone2 = 2,
#   method = "wilcox"
#   # method = "DESeq2"
# )




# DE plot
DE_full_plot = Rcongas::plot_DE_volcano(fit, annotate_top = 7, cut_pvalue = 0.05) +
  labs(title = "Differential Expression")

# Mixing
mixing_plot = Rcongas::plot_mixing_proportions(fit)

# Whole-genome plot
CNA_wg_plot = Rcongas::plot_gw_cna_profiles(fit, whole_genome = TRUE, cutoff_p = 0.001) +
  labs(title = "Glioblastoma (Smart-Seq scRNAseq)")

# Chromosome 15 special counts plot
segments_ids = Rcongas::get_input_segmentation(fit) %>% Rcongas:::idify() %>% pull(segment_id)

counts_plot = Rcongas::plot_segment_density(
  fit,
  segments_ids = c(
    "chr2:272203:219146868",
    "chr5:102886512:162944709",
    "chr5:224633:180278393",
    "chr10:3154461:133786586",
    "chr7:44444122:158715219",
    "chr7:881668:44121936",
    "chr13:21086599:101294428"
  )
)

# plot_segment_density(fit, "chr8:1:67500000")

# z-score input RNA with clustering assignments
rna_plot_raw = Rcongas::plot_counts_rna_segments(fit, z_score = TRUE, cutoff_p = 0.001)

# Comparison with clonealign
congas_clusterings = Rcongas::get_clusters(fit) %>% select(cell, cluster) %>%
  rename(congas = cluster) %>% as_tibble()

load("~/Documents/GitHub/rcongas/data/intersection_congas_clonealign.rda")


# Strips
top_strip =  plot_grid(
  CNA_wg_plot,
  mixing_plot,
  nrow = 1,
  rel_widths = c(1, .3),
  align = 'h',
  labels = c("a", "b")
)

mid_strip =  plot_grid(
  rna_plot_raw,
  # counts_plot,
  DE_full_plot,
  nrow = 1,
  rel_widths = c(.7, .3),
  align = 'h',
  axis = 'b',
  labels = c("c", "d")
)

bottom_strip = CNAqc:::eplot()

bottom_strip = plot_grid(
  counts_plot$`chr15:1:102600000`,
  counts_plot$`chr18:7950001:55950000`,
  clonalign_comparison_plot,
  axis = 'b',
  nrow = 1,
  labels = c("e", "f", "g")
)


plot_grid(
  top_strip,
  mid_strip,
  bottom_strip,
  rel_heights = c(.8, 1, .8),
  ncol = 1
) %>%
  ggsave(filename = "GBM_SmartSeq_Main_text.pdf",
         width = 12,
         height = 10)

plot_grid(
  top_strip,
  mid_strip,
  bottom_strip,
  rel_heights = c(.8, 1, .8),
  ncol = 1
) %>%
  ggsave(filename = "breast_xeno_10x_Main_text.png",
         width = 12,
         height = 10)


# plot_DE_gw(fit)
# plot_latent_variables(fit)
