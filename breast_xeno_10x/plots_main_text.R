require(Rcongas)
require(tidyverse)
library(cowplot)

# Final model fit
fit = Rcongas::breast_xeno_10x_small_segments_total_CN_norm
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
DE_full_plot = Rcongas::plot_DE_volcano(fit, annotate_top = 7) +
  labs(title = "Differential Expression")

# Mixing
mixing_plot = Rcongas::plot_mixing_proportions(fit)

# Whole-genome plot
CNA_wg_plot = Rcongas::plot_gw_cna_profiles(fit, whole_genome = TRUE, cutoff_p = 0.001) +
  labs(title = "Breast cancer xenograft (10x scRNAseq)")

# Chromosome 15 special counts plot
segments_ids = Rcongas::get_clones_ploidy(fit) %>%
  Rcongas:::idify() %>%
  filter(highlight) %>%
  pull(segment_id) %>%
  unique()

counts_plot = Rcongas::plot_segment_density(
  fit,
  segments_ids = segments_ids
)
names(counts_plot) = segments_ids

# plot_segment_density(fit, "chr8:1:67500000")

# z-score input RNA with clustering assignments
rna_plot_raw = Rcongas::plot_counts_rna_segments(fit, z_score = TRUE, cutoff_p = 0.001)

# Comparison with clonealign
congas_clusterings = Rcongas::get_clusters(fit) %>%
  select(cell, cluster) %>%
  rename(congas = cluster) %>% as_tibble()

# Salvatore
load("../annealToolbox/data/intersection_congas_clonealign.rda")

# Giulio
intersection_congas_clonealign = Rcongas::intersection_congas_clonealign

intersection_congas_clonealign$cluster <- Rcongas:::get_cluster_assignments(fit)

intersection_congas_clonealign = intersection_congas_clonealign %>%
  rename(CONGAS = cluster, clonealign = clone)

# clonealign_clusterings = Rcongas::get_clusters(fit) %>% select(cell, cluster) %>%
#   rename(clonealign = cluster) %>% as_tibble()
joint_clusters = intersection_congas_clonealign %>%
  group_by(CONGAS, clonealign) %>%
  summarise(n = n()) %>% mutate(CONGAS = paste(CONGAS))


# require(infotheo)
# mutinformation(
#   intersection_congas_clonealign$CONGAS,
#   intersection_congas_clonealign$clonealign
# ) /
#   sqrt(
#     entropy(intersection_congas_clonealign$CONGAS) * entropy(intersection_congas_clonealign$clonealign)
#   )

# joint_clusters = reshape2::melt(joint_clusters)

require(ggalluvial)

require(pdfCluster)

clonalign_comparison_plot = ggplot(joint_clusters,
                                   aes(y = n, axis1 = CONGAS, axis2 = clonealign)) +
  geom_alluvium(aes(fill = CONGAS), width = 1 / 12) +
  geom_stratum(width = 1 / 12,
               fill = "black",
               color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("CONGAS", "clonealign"),
                   expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  CNAqc:::my_ggplot_theme() +
  guides(fill = guide_legend("cluster")) +
  labs(title = "Comparison with clonealign", y = 'Cell',
       caption = paste0('Adjusted Rand Index:', round(pdfCluster::adj.rand.index(intersection_congas_clonealign$CONGAS,
                                                                           intersection_congas_clonealign$clonealign), digits = 2), "."))



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


bottom_strip = plot_grid(
  counts_plot$`chr15:67050001:102600000`,
  counts_plot$`chr16:1:3750000`,
  counts_plot$`chr18:32400001:55950000`,
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
  ggsave(filename = "breast_xeno_10x_Main_text.pdf",
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
