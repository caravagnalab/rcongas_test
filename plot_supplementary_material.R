require(Rcongas)
require(tidyverse)
library(cowplot)

# Final model fit
fit = Rcongas::GBM_smartseq_normal_vs_tumor
input_raw_counts_genes = Rcongas::get_input_raw_data(fit)

# DE genes raw table
Rcongas::plot_raw_data(
  get_input_raw_data(fit),
  genes = Rcongas::get_DE_table(fit, cut_pvalue = 0.05) %>% pull(gene),
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
