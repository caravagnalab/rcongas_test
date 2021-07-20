require(Rcongas)
require(tidyverse)
library(cowplot)

# Final model fit
fit_nt = Rcongas::GBM_smartseq_normal_vs_tumor

fit_t_14 = Rcongas::GBM_smartseq_tumor_only_df_14_no14
fit_t_5 = Rcongas::GBM_smartseq_tumor_only_df_5_no5

fit_t = fit_t_5

# input_raw_counts_genes = Rcongas::get_input_raw_data(fit)

# DE plot
DE_full_plot_nt = Rcongas::plot_DE_volcano(fit_nt, annotate_top = 7, cut_pvalue = 0.05) +
  labs(title = "Differential Expression (Tumour/ normal)")

DE_full_plot_t = Rcongas::plot_DE_volcano(fit_t, annotate_top = 7, cut_pvalue = 0.05) +
  labs(title = "Differential Expression (Tumour)")

# Tumour mixing
mixing_plot_t = Rcongas::plot_mixing_proportions(fit_t) +
  scale_fill_brewer(palette = 'Dark2')

# DE genes raw table
Rcongas::plot_raw_data(
  x = get_input_raw_data(fit_nt),
  # genes = (get_input_raw_data(fit_nt) %>% rownames())[174:233],
  genes = Rcongas::get_DE_table(fit_nt, cut_pvalue = 1) %>% pull(gene),
  lognormalise = FALSE,
  clusters = Rcongas::get_clusters(fit_nt),
  prompt = F,
  fontsize_row = 6,
  filename = "SM_gbm_rawdata.png",
  width = 11,
  height = 7.5
)

# Empty plot
# gw_DE_plot = plot_DE_gw(fit_t, cut_pvalue = 0.05)

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

# ggarrange(
#   gw_DE_plot,
#   all_counts_plot,
#   nrow = 2,
#   heights = c(1.5, 4),
#   labels = c('a', 'b')
# ) %>%
#   ggsave(filename = "SM_breast_xeno_10x_counts.png",
#          width = 13,
#          height = 16)


# plot_latent_variables(fit)

Rcongas:::plot_inference_report(fit_t) %>%
  ggsave(filename = "SM_GBM_inference_report.png",
         width = 10,
         height = 12)
