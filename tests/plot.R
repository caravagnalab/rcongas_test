# Main_text score
main_text_score = "ARI"
supp_mat_scores = c("ARI", "MI", "NMI", "match")

##################
## TEST  1
##################

load("test1.RData")

ls()

res = res %>%
  tidyr::separate(param_1, into = c('p1', 'Distance'), sep ='=') %>%
  tidyr::separate(param_2, into = c('p2', 'Clusters'), sep ='=')
# %>%
#   mutate(param = factor(param, levels = sort(res$param %>% unique)))

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

res %>%
  reshape2::melt(id = c("Distance", "Clusters")) %>%
  filter(variable %in% supp_mat_scores) %>%
  # sample_n(600) %>%
  ggplot(aes(
    x = paste(Distance),
    y = paste(Clusters),
    fill = as.numeric(value)
  )) +
  geom_tile() +
  scale_fill_distiller("Score", palette = 'Blues', direction = 1) +
  facet_wrap(~ variable, scales = 'free') +
  CNAqc:::my_ggplot_theme() +
  labs(title = 'Test 1') +
  ggsave(filename = "SM_test1.png",
         width = 15,
         height = 10)

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
      y = c(.4),
      label = c('high dispersion')
    ),
    aes(x = x, y = y, label = label),
    color = 'indianred3',
    size = 3
  ) +
  geom_text(
    data = data.frame(
      x = c(8),
      y = c(.4),
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
      y = c(.5),
      label = c('few genes')
    ),
    aes(x = x, y = y, label = label),
    color = 'indianred3',
    size = 3
  ) +
  geom_text(
    data = data.frame(
      x = c(8),
      y = c(.5),
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
  CNAqc:::eplot(),
  plot_grid(
    test2_main_text,
    test3_main_text,
    align = 'h',
    axis = 'b',
    nrow = 1,
    labels = c('b', 'c')),
  nrow = 2,
  labels = c('a')
)+
  ggsave("MT_Figura_tests.png", width = 8, height = 6)

