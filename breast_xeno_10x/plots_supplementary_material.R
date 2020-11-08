require(Rcongas)
require(tidyverse)
library(cowplot)

# Final model fit
fit = Rcongas::breast_xeno_10x_small_segments_total_CN_norm
input_raw_counts_genes = Rcongas::get_input_raw_data(fit)

# DE genes raw table
Rcongas::plot_raw_data(
  get_input_raw_data(fit) %>%  t(),
  genes = Rcongas::get_DE_table(fit, cut_pvalue = 0.001) %>% pull(gene),
  clusters = Rcongas::get_clusters(fit),
  prompt = F,
  fontsize_row = 6,
  filename = "SM_breast_xeno_10x_rawdata.png",
  width = 11,
  height = 7.5
)

gw_DE_plot = plot_DE_gw(fit)

counts_plot = Rcongas::plot_segment_density(
  fit,
  segments_ids = Rcongas::get_input_segmentation(fit) %>% Rcongas:::idify() %>% pull(segment_id)
)

all_counts_plot = ggarrange(
  plotlist = counts_plot,
  common.legend = TRUE,
  legend = 'bottom',
  nrow = 5,
  ncol = 6
)

ggarrange(
  gw_DE_plot,
  all_counts_plot,
  nrow = 2,
  heights = c(1.5, 4),
  labels = c('a', 'b')
) %>%
  ggsave(filename = "SM_breast_xeno_10x_counts.png",
         width = 13,
         height = 16)


# plot_latent_variables(fit)

Rcongas:::plot_inference_report(fit) %>%
  ggsave(filename = "SM_breast_xeno_10x_inference_report.png",
         width = 10,
         height = 12)

# Comparison clonealign
congas_clusterings = Rcongas::get_clusters(fit) %>% select(cell, cluster) %>%
  rename(congas = cluster) %>% as_tibble()

load("~/Documents/GitHub/rcongas/data/intersection_congas_clonealign.rda")

intersection_congas_clonealign = intersection_congas_clonealign %>%
  rename(CONGAS = cluster, clonealign = clone)

# clonealign_clusterings = Rcongas::get_clusters(fit) %>% select(cell, cluster) %>%
#   rename(clonealign = cluster) %>% as_tibble()
joint_clusters = intersection_congas_clonealign %>%
  group_by(CONGAS, clonealign) %>%
  summarise(n = n())

# Where we differ
plot_comparison_noi_clonealign = function(fit, intersection_congas_clonealign, chr)
{
  I = pio:::nmfy(
    intersection_congas_clonealign$cell,
    intersection_congas_clonealign$clonealign
  )

  M = Rcongas::get_counts(fit, chr) %>%
    mutate(cluster_clonealign = I[cell]) %>%
    Rcongas:::idify()

  plot_grid(
    ggplot(M, aes(n, fill = cluster_clonealign)) + geom_histogram(bins = 100) + CNAqc:::my_ggplot_theme() + scale_fill_brewer(palette = "Set2") + facet_wrap( ~
                                                                                                                                                                segment_id) + labs(title = paste0("Chromosome ", chr)),
    ggplot(M, aes(n, fill = cluster)) + geom_histogram(bins = 100) + CNAqc:::my_ggplot_theme() + scale_fill_brewer(palette = "Set1") + facet_wrap( ~
                                                                                                                                                     segment_id),
    ncol = 1
  )
}


ggarrange(
  plot_comparison_noi_clonealign(fit, intersection_congas_clonealign, 'chr1'),
  plot_comparison_noi_clonealign(fit, intersection_congas_clonealign, 'chr6'),
  plot_comparison_noi_clonealign(fit, intersection_congas_clonealign, 'chr8'),
  plot_comparison_noi_clonealign(fit, intersection_congas_clonealign, 'chr9'),
  plot_comparison_noi_clonealign(fit, intersection_congas_clonealign, 'chr15'),
  plot_comparison_noi_clonealign(fit, intersection_congas_clonealign, 'chr18'),
  labels = c('a', 'b', 'c', 'd', 'e', 'f')
) %>%
  ggsave(filename = "SM_breast_xeno_10x_comparison_clonealign.png",
         width = 10,
         height = 12)
