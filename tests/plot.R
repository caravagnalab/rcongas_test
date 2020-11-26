# Main_text score
main_text_score = "ARI"
supp_mat_scores = c("ARI", "NMI")

##################
## TEST  1
##################

load("test1.RData")

ls()

res = res %>%
  tidyr::separate(param_1, into = c('p1', 'Distance'), sep ='=') %>%
  tidyr::separate(param_2, into = c('p2', 'Clusters'), sep ='=')

# res = res %>% rename(Distance = param_1, Clusters = param_2)
# res %>%
#   reshape2::melt(id = c("Distance", "Clusters")) %>%
#   filter(variable %in% supp_mat_scores) %>%
#   as_tibble() %>%
#   rename(score = variable) %>%
#   ggplot(aes(
#     x = paste(Clusters),
#     y =  as.numeric(value),
#     fill = (score),
#     group = interaction(score, Clusters),
#   )) +
#   geom_boxplot() +
#   facet_wrap(~Distance, scales = 'free')
means_t1 = res %>%
  reshape2::melt(id = c("Distance", "Clusters")) %>%
  group_by(Distance, Clusters,variable) %>%
  summarize(score = mean(as.numeric(value), na.rm = T))
colnames(means_t1)[4] = 'mean'

sd_t1 = res %>%
  reshape2::melt(id = c("Distance", "Clusters")) %>%
  group_by(Distance, Clusters,variable) %>%
  summarize(score = sd(as.numeric(value), na.rm = T))
colnames(sd_t1)[4] = 'sd'

t1_data = full_join(means_t1, sd_t1)

test1_main_text = t1_data %>%
  filter(variable %in% main_text_score) %>%
  # sample_n(600) %>%
  ggplot(aes(
    x = Distance,
    y = Clusters,
    fill = as.numeric(mean)
  )) +
  geom_tile() +
  scale_fill_distiller("", palette = 'Blues', direction = 1) +
  facet_wrap(~ paste(variable, " (mean and standard deviation)"), scales = 'free') +
  geom_text(aes(label = format(sd, scientific = T, digits = 3)), color = 'white', size = 3) +
  CNAqc:::my_ggplot_theme() +
  labs(title = 'Different clonal architectures', x = "Evolutionary distance") +
  guides(fill = guide_colorbar(barheight = unit(4, 'cm'))) +
  theme(legend.position = 'right')


all_points_scatter = res %>%
  reshape2::melt(id = c("Distance", "Clusters")) %>%
  filter(variable %in% 'match') %>%
  as_tibble()

t1_scatter = ggplot(data = all_points_scatter,
       aes(
         x = as.numeric(Clusters),
         y = as.numeric(value),
         fill = Distance
       )) +
  labs(
    x = "Clusters",
    y = "% matched cluster labels"
  ) +
  geom_jitter(size = .5,
              alpha = .3, aes(color = Distance)) +
  facet_wrap(~ paste('Distance', Distance)) +
  geom_errorbar(
    data = t1_data %>% filter(variable == 'match'),
    aes(
      x = as.numeric(Clusters),
      color = Distance,
      ymin=mean-sd,
      ymax=mean+sd),
    width = .2,
    inherit.aes = FALSE,
    show.legend = F
  ) +
  geom_point(
    data = t1_data %>% filter(variable == 'match'),
    aes(
      x = as.numeric(Clusters),
      y = mean,
      fill = Distance
    ),
    color = 'black',
    inherit.aes = FALSE,
    show.legend = F
  ) +
  CNAqc:::my_ggplot_theme() +
  ggsci::scale_color_jama() +
  scale_x_continuous(breaks = 1:10) +
  guides(color = F, fill = F) +
  theme(legend.position = 'right')


library(ggridges)

t1_scatter =  ggplot(data = all_points_scatter %>% mutate(value = as.numeric(value)),
                     aes(
                       y = Clusters,
                       x = value,
                       fill = Distance
                     )) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 1.5, point_alpha = 1, alpha = 0.7
  ) +
  facet_wrap(~ paste('Distance', Distance)) +
  labs(
    y = "Clusters",
    x = "% matched cluster labels"
  ) +
  CNAqc:::my_ggplot_theme() +
  ggsci::scale_fill_jama()+
  theme(legend.position = 'right') +
  guides(fill=F)

# Supp Mat
s1f1 = res %>%
  ggplot(aes(match)) +
  facet_grid(Distance~Clusters, scales = 'free') +
  geom_histogram(binwidth = .1) +
  scale_x_continuous(breaks = scales::pretty_breaks(3)) +
  scale_fill_viridis_c("NMI") +
  CNAqc:::my_ggplot_theme() +
  labs(x = "% of matches",
       title = 'Percentage of strict matches for the main text')

s1f2 = res %>%
  ggplot(aes(NMI)) +
  facet_grid(Distance~Clusters, scales = 'free') +
  geom_histogram(binwidth = 0.01) +
  scale_x_continuous(breaks = scales::pretty_breaks(3)) +
  scale_fill_viridis_c("NMI") +
  CNAqc:::my_ggplot_theme() +
  labs(x = "Normalised Mutual Information (NMI)",
       title = 'Normalised Mutual Information for the main test')

cowplot::plot_grid(s1f1, s1f2, labels= c("a", "b"), nrow = 2) +
  ggsave("SM_test1.png", width = 8, height = 6)


##################
## TEST  2
##################

load("test2.RData")

ls()

res = res %>%
  mutate(param = factor(param, levels = sort(res$param %>% unique)))

res %>%
  reshape2::melt(id = "param") %>%
  ggplot(aes(
    x = as.factor(param),
    value,
    fill = variable,
    group = paste0(param)
  )) +
  geom_jitter(size = .4, color = 'steelblue') +
  geom_boxplot() +
  facet_wrap(~ variable, scales = 'free_y') +
  CNAqc:::my_ggplot_theme() +
  labs(x = "Parameter", title = 'Test 2') +
  ggsave(filename = "test_2_all.pdf",
         width = 15,
         height = 10)

res %>%
  reshape2::melt(id = "param") %>%
  ggplot(aes(
    x = as.factor(param),
    value,
    fill = variable,
    group = paste0(param)
  )) +
  geom_jitter(size = .4, color = 'steelblue') +
  geom_boxplot() +
  facet_wrap(~ variable, scales = 'free_y') +
  CNAqc:::my_ggplot_theme() +
  labs(x = "Parameter", title = 'Test 2') +
  ggsave(filename = "test_2_all.pdf",
         width = 15,
         height = 10)

# Main text
df_summary <- res %>%
  select(param, main_text_score) %>%
  group_by(param) %>%
  summarise_all(list(sd, mean),) %>%
  rename(sd = fn1, mean = fn2)

test2_main_text = ggplot() +
  geom_jitter(
    data = res,
    aes(x = param, y = eval(parse(text = main_text_score))),
    size = .4,
    color = 'darkgray'
  ) +
  geom_pointrange(data = df_summary,
                  aes(
                    x = param,
                    y = mean,
                    ymin = mean - sd,
                    ymax = mean + sd
                  ),
                  color = 'darkred') +
  scale_color_continuous() +
  CNAqc:::my_ggplot_theme() +
  labs(x = "Negative Binomial size",
       y = main_text_score,
       title = 'Sequencing dispersion') +
  geom_text(
    data = data.frame(
      x = c(2),
      y = c(1),
      label = c('high dispersion')
    ),
    aes(x = x, y = y, label = label),
    color = 'indianred3',
    size = 3
  ) +
  geom_text(
    data = data.frame(
      x = c(8),
      y = c(1),
      label = c('low dispersion')
    ),
    aes(x = x, y = y, label = label),
    # hjust = 0,
    color = 'steelblue',
    size = 3
  )

# Supp Mat
res %>%
  select(param, supp_mat_scores) %>%
  reshape2::melt(id = "param") %>%
  ggplot(aes(x = param, y = value)) +
  geom_jitter(size = .4,
              color = 'darkgray') +
  geom_boxplot(fill = NA) +
  scale_color_continuous() +
  facet_wrap( ~ variable) +
  CNAqc:::my_ggplot_theme() +
  labs(x = "Negative Binomial size",
       y = "Score",
       title = 'Sequencing dispersion (extended tests)') +
  ggsave("SM_test2.png", width = 8, height = 6)

##################
## TEST  3
##################

load("test3.RData")

ls()

res = res %>%
  mutate(param = factor(param, levels = sort(res$param %>% unique)))

res %>%
  reshape2::melt(id = "param") %>%
  ggplot(aes(
    x = as.factor(param),
    value,
    fill = variable,
    group = paste0(param)
  )) +
  geom_jitter(size = .4, color = 'steelblue') +
  geom_boxplot() +
  facet_wrap(~ variable, scales = 'free_y') +
  CNAqc:::my_ggplot_theme() +
  labs(x = "Parameter", title = 'Test 3') +
  ggsave(filename = "test_3_all.pdf",
         width = 15,
         height = 10)

# Main text
df_summary <- res %>%
  select(param, main_text_score) %>%
  group_by(param) %>%
  summarise_all(list(sd, mean),) %>%
  rename(sd = fn1, mean = fn2)

test3_main_text =
  ggplot() +
  geom_jitter(
    data = res,
    aes(x = param, y = eval(parse(text = main_text_score))),
    size = .4,
    color = 'darkgray'
  ) +
  geom_pointrange(data = df_summary,
                  aes(
                    x = param,
                    y = mean,
                    ymin = mean - sd,
                    ymax = mean + sd
                  ),
                  color = 'darkred') +
  scale_color_continuous() +
  CNAqc:::my_ggplot_theme() +
  labs(x = "% of genes with different CNA",
       y = main_text_score,
       title = 'Bad input segmentation') +
  geom_text(
    data = data.frame(
      x = c(2),
      y = c(1),
      label = c('few genes')
    ),
    aes(x = x, y = y, label = label),
    color = 'indianred3',
    size = 3
  ) +
  geom_text(
    data = data.frame(
      x = c(8),
      y = c(1),
      label = c('many genes')
    ),
    aes(x = x, y = y, label = label),
    color = 'steelblue',
    size = 3
  )

# Supp Mat
res %>%
  select(param, supp_mat_scores) %>%
  reshape2::melt(id = "param") %>%
  ggplot(aes(x = param, y = value)) +
  geom_jitter(size = .4,
              color = 'darkgray') +
  geom_boxplot(fill = NA) +
  scale_color_continuous() +
  facet_wrap( ~ variable) +
  CNAqc:::my_ggplot_theme() +
  labs(x = "% of genes with different CNA",
       y = "Score",
       title = 'Bad input segmentation (extended tests)') +
  ggsave("SM_test3.png", width = 8, height = 6)

require(cowplot)

plot_grid(
  plot_grid(test1_main_text, t1_scatter, nrow = 1, align = 'h', axis = 'bt',     labels = c('a', 'b')),
  plot_grid(
    test2_main_text,
    test3_main_text,
    align = 'h',
    axis = 'bt',
    nrow = 1,
    labels = c('c', 'd')),
  nrow = 2,
  labels = c('a')
)+
  ggsave("MT_Figura_tests.png", width = 10, height = 8)







