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

# Mixing
mixing_plot_nt = Rcongas::plot_mixing_proportions(fit_nt)

# Whole-genome plot
CNA_wg_plot_nt = Rcongas::plot_gw_cna_profiles(fit_nt, whole_genome = F, alpha = 0.1) +
  labs(title = "Glioblastoma (Smart-Seq scRNAseq), CONGAS with normal and tumour cells") +
  theme(legend.position = 'right')

CNA_wg_plot_t = Rcongas::plot_gw_cna_profiles(fit_t, whole_genome = F, alpha = 0.1) +
  labs(title = "CONGAS with only tumour cells (second run)") +
  scale_color_brewer(palette = 'Set2')+
  theme(legend.position = 'right')

plot_highlights(fit_nt, alpha = 0.1)
plot_highlights(fit_t, alpha = 0.1)

# Chromosome special counts plot
segments_ids_nt = Rcongas::get_clones_ploidy(fit_nt, alpha = 0.1) %>%
  Rcongas:::idify() %>%
  filter(highlight) %>%
  pull(segment_id) %>%
  unique

segments_ids_t = Rcongas::get_clones_ploidy(fit_t) %>%
  Rcongas:::idify() %>%
  filter(highlight) %>%
  pull(segment_id) %>%
  unique()

counts_plot = NULL

counts_plot_nt = Rcongas::plot_segment_density(
  fit_nt,
  segments_ids = segments_ids_nt,
  sum_denominator = FALSE
)

counts_plot_t = Rcongas::plot_segment_density(
  fit_t,
  segments_ids = segments_ids_t,
  sum_denominator = FALSE
)

counts_plot = append(
  lapply(counts_plot_nt, function(x) x + ggtitle("First run (Gaussian)")),
  lapply(counts_plot_t, function(x) x + ggtitle("Second run (Gaussian)"))
)

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

# z-score input RNA with clustering assignments
rna_plot_raw_nt = Rcongas::plot_counts_rna_segments(fit_nt,
                                                    normalised = FALSE,
                                                    sum_denominator = FALSE,
                                                    z_score = TRUE,
                                                    alpha = 0.1)

# Strips
top_strip =  plot_grid(
  plot_grid(CNA_wg_plot_nt, mixing_plot_nt, nrow = 1, rel_widths = c(1, .6), align = 'h', axis = 'bt', labels = c('a', 'b')),
  # CNA_wg_plot_nt,
  plot_grid(CNA_wg_plot_t, runs_comparison_plot, nrow = 1, rel_widths = c(1, .6), align = 'h', axis = 'bt', labels = c('c', 'd')),
  plot_grid(plotlist = append(list(rna_plot_raw_nt), counts_plot), rel_widths = c(1, .5, .5, .5), nrow = 1, labels = letters[5:10]),
  ncol = 1,
  labels = c("a", "b", '')
)


top_strip %>%
  ggsave(filename = "GBM_SmartSeq_Main_text.png",
         width = 12,
         height = 10)



# plot_DE_gw(fit)
# plot_latent_variables(fit)






# plot_segment_density(fit, "chr8:1:67500000")



rna_plot_raw_t = Rcongas::plot_counts_rna_segments(fit_t,
                                                   normalised = FALSE,
                                                   sum_denominator = FALSE,
                                                   z_score = TRUE,
                                                   alpha = 0.1)

