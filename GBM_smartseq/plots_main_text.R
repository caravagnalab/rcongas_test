require(Rcongas)
require(tidyverse)
library(cowplot)

# Final model fit
fit_nt = Rcongas::GBM_smartseq_normal_vs_tumor

fit_t_14 = Rcongas::GBM_smartseq_tumor_only_df_14_no14
fit_t_5 = Rcongas::GBM_smartseq_tumor_only_df_5_no5

fit_t = fit_t_5

# fit = subclones
input_raw_counts_genes_nt = Rcongas::get_input_raw_data(fit_nt)
input_raw_counts_genes_t = Rcongas::get_input_raw_data(fit_t)

# DE plot
DE_full_plot_nt = Rcongas::plot_DE_volcano(fit_nt, annotate_top = 7, cut_pvalue = 0.05) +
  labs(title = "Differential Expression (Tumour/ normal)")

DE_full_plot_t = Rcongas::plot_DE_volcano(fit_t, annotate_top = 7, cut_pvalue = 0.05) +
  labs(title = "Differential Expression (Tumour)")

# Mixing
mixing_plot_nt = Rcongas::plot_mixing_proportions(fit_nt)
mixing_plot_t = Rcongas::plot_mixing_proportions(fit_t)

# Whole-genome plot
CNA_wg_plot_nt = Rcongas::plot_gw_cna_profiles(fit_nt, whole_genome = F, alpha = 0.1) +
  labs(title = "Glioblastoma (Smart-Seq scRNAseq), running CONGAS with normal (cluster 1) and tumour cells") +
  theme(legend.position = 'right')
CNA_wg_plot_t = Rcongas::plot_gw_cna_profiles(fit_t, whole_genome = F, alpha = 0.1) +
  labs(title = "Run with just tumour cells") +
  scale_color_brewer(palette = 'Set2')+
  theme(legend.position = 'right')

plot_highlights(fit_nt, alpha = 0.1)
plot_highlights(fit_t, alpha = 0.1)

# Mapping clustering assignments
first_run = data.frame(
  cell = Rcongas:::get_cluster_assignments(fit_nt) %>% names,
  first_cluster = Rcongas:::get_cluster_assignments(fit_nt)
)

second_run = data.frame(
  cell = Rcongas:::get_cluster_assignments(fit_t) %>% names,
  second_cluster = Rcongas:::get_cluster_assignments(fit_t)
)

joint_clusters = first_run %>%
  full_join(second_run) %>%
  mutate(second_cluster = ifelse(is.na(second_cluster), "Missing", second_cluster)) %>%
  group_by(first_cluster, second_cluster) %>%
  summarise(n = n())

require(ggalluvial)

runs_comparison_plot = ggplot(joint_clusters,
                                   aes(y = n, axis1 = first_cluster, axis2 = second_cluster)) +
  geom_alluvium(aes(fill = second_cluster), width = 1 / 12) +
  geom_stratum(width = 1 / 12,
               fill = "black",
               color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Tumour & Normal", "Tumour"),
                   expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  CNAqc:::my_ggplot_theme() +
  guides(fill = guide_legend("Clusters (second run)")) +
  labs(title = "Comparison across runs", y = 'Cell')


plot_grid(
  CNA_wg_plot_nt,
  plot_grid(CNA_wg_plot_t, runs_comparison_plot, nrow = 1, rel_widths = c(1, .6)),
  ncol = 1
)

# Chromosome special counts plot
segments_ids = Rcongas::get_clones_ploidy(fit_nt, alpha = 0.1) %>%
  Rcongas:::idify() %>%
  filter(highlight) %>%
  pull(segment_id)

# segments_ids = Rcongas::get_clones_ploidy(fit) %>%
#   Rcongas:::idify() %>%
#   filter(highlight) %>%
#   pull(segment_id) %>%
#   unique()

counts_plot = NULL
counts_plot = Rcongas::plot_segment_density(
  fit_nt,
  segments_ids = segments_ids
)

# plot_segment_density(fit, "chr8:1:67500000")

# z-score input RNA with clustering assignments
rna_plot_raw_nt = Rcongas::plot_counts_rna_segments(fit_nt,
                                                    normalised = FALSE,
                                                    sum_denominator = FALSE,
                                                    z_score = TRUE,
                                                    alpha = 0.1)

rna_plot_raw_t = Rcongas::plot_counts_rna_segments(fit_t,
                                                   normalised = FALSE,
                                                   sum_denominator = FALSE,
                                                   z_score = TRUE,
                                                   alpha = 0.1)


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

# bottom_strip = plot_grid(
#   counts_plot$`chr15:1:102600000`,
#   counts_plot$`chr18:7950001:55950000`,
#   clonalign_comparison_plot,
#   axis = 'b',
#   nrow = 1,
#   labels = c("e", "f", "g")
# )


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



# plot_DE_gw(fit)
# plot_latent_variables(fit)
