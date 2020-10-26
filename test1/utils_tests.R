
torch <- reticulate::import("torch")

run_example <-  function(arg_list_1, arg_list_2){

  cnv_df <- do.call(generate_cluster_ploidy_df, arg_list_1)
  simul <-  do.call(run_simulation_generative, c(list(cnv_df = cnv_df), arg_list_2))

  res <- run_inference(simul, "MixtureGaussianDMP", param_list = list("theta_scale" = 16, "theta_rate" = 16, "cnv_var" = 0.95,"alpha" = 0.0001, 'T' = 5L), lr = 0.1)

  return(list(inference = res, simulation = simul))

}


subset_sim <- function(simul, nsegs) {
  pld_real <- simul$data$cnv %>% select(matches("ploidy")) %>% select(-ploidy_real)
  equal_segs <- which(apply(pld_real,1, function(x) (var(x) == 0)))
  diff_segs <-   which(apply(pld_real,1, function(x) (var(x) != 0)))
  equal_segs <- sample(equal_segs, size = nsegs)
  simul <- simul[,c(diff_segs, equal_segs)]

  print(simul[,c(diff_segs, equal_segs)])

  return(simul)
}

run_example_fixed <-  function(arg_list_1, arg_list_2, IC = "BIC", nsegs = NULL, posteriors = FALSE, steps = 300L){

  cnv_df <- do.call(generate_cluster_ploidy_df, arg_list_1)
  simul <-  do.call(run_simulation_generative, c(list(cnv_df = cnv_df), arg_list_2))
  mixture <-  list(NULL, torch$tensor(c(0.7,0.3)), torch$tensor(c(0.4,0.3,0.3)), torch$tensor(c(0.3,0.3,0.2,0.2)) )

  if(!is.null(nsegs)){
    simul <- subset_sim(simul, nsegs)
  }
  simul$cnv$ploidy_real <-  round(simul$data$cnv$ploidy_real)

  res <- best_cluster(simul, "MixtureGaussian", 1:4,param_list = list("theta_scale" = 32, "theta_rate" = 32, "cnv_var" = 0.65),steps = steps, lr = 0.1, mixture = mixture, method = IC, posteriors = posteriors)

  #res <-  NULL

  return(list(inference = res$models[[which.min(res$IC)]], simulation = simul, IC = res$IC))

}



plot_test_1 <- function(test_res, p2s = 1:2,font_size = 10) {
  p1 <- pheatmap::pheatmap(test_res$CNVs, cluster_rows = F, cluster_cols = F, breaks=seq(0, 1, length.out=101),fontsize = font_size )[[4]]
  #p2 <- pheatmap::pheatmap(test_res$clusters[,p2s], cluster_rows = F, cluster_cols = F)[[4]]
  #do.call(gridExtra::grid.arrange, list(p1,p2, ncol = 2))
  return(p1)
}


plot_test_5 <- function(test_res) {

  data <- reshape2::melt(test_res$CNVs)
  data <- cbind(data, sd = as.vector(test_res$table_pred_sd))
  colnames(data) <- c("Zeta", "Changes", "NMI", "sd")
  ggplot(data = data, aes(y = NMI, x = factor(Zeta), color = Changes, group = Changes)) + geom_point(position=position_dodge(width = 0.6)) +
    geom_errorbar(aes(ymin=NMI-sd, ymax=NMI+sd), width=.1,
                  position=position_dodge(0.6)) + xlab("Size") + ylab("NMI") +
    ggtitle("NMI vs NB overdispersion") +  theme_bw() + scale_color_discrete("CNVs", labels = c("Amplification", "Duplication"))


}



plot_test_2 <- function(test_res) {

  data <- reshape2::melt(test_res$CNVs)
  data <- cbind(data, sd = as.vector(test_res$table_pred_sd))
  colnames(data) <- c("Zeta","Changes", "NMI", "sd")
  print(data)
  ggplot(data = data, aes(y = NMI, x = factor(Zeta))) + geom_point(color = "blue") +
    geom_errorbar(aes(ymin=NMI-sd, ymax=NMI+sd), width=.1,
                  position=position_dodge(0.6)) + xlab("Size") + ylab("NMI") + ylim(c(0,1.2)) +
    ggtitle("NMI vs #neutral segment") +  theme_bw()


}

plot_test_3 <- function(test_res) {

  data <- reshape2::melt(test_res$CNVs)
  data <- cbind(data, sd = as.vector(test_res$table_pred_sd))
  colnames(data) <- c("Zeta", "Changes", "NMI", "sd")
  ggplot(data = data, aes(y = NMI, x = factor(Zeta), color = Changes, group = Changes)) + geom_point(position=position_dodge(width = 0.6)) +
    geom_errorbar(aes(ymin=NMI-sd, ymax=NMI+sd), width=.1,
                  position=position_dodge(0.6)) + xlab("%genes") + ylab("NMI") +
    ggtitle("NMI vs subclonal gene %") +  theme_bw() + scale_color_discrete("CNVs", labels = c("Amplification", "Duplication"))


}


plot_test_6 <- function(test_res) {

  data <- reshape2::melt(test_res$CNVs)
  data <- cbind(data, as.vector(test_res$table_pred_sd))
  colnames(data) <- c("Minor_freq", "Changes", "NMI", "sd")

  ggplot(data = data, aes(y = NMI, x = factor(Minor_freq), color = factor(Changes))) + geom_point(position = position_dodge(width=0.75)) +
    geom_errorbar(aes(ymin=NMI-sd, ymax=NMI+sd), width=.1,
                  position = position_dodge(width=0.75)) + xlab("Minor cluster size") + ylab("NMI") + scale_color_discrete("CN difference")  +
    ggtitle("NMI vs size of minor clusers") + theme_bw() + scale_color_discrete("CNVs", labels = c("Amplification", "Duplication"))

}


plot_test_4 <- function(test_res) {

  data <- reshape2::melt(test_res$CNVs)
  data <- cbind(data, as.vector(test_res$table_pred_sd))
  colnames(data) <- c("Minor_freq", "Changes", "NMI", "sd")

  ggplot(data = data, aes(y = NMI, x = factor(Minor_freq), color = factor(Changes))) + geom_point(position = position_dodge(width=0.75)) +
    geom_errorbar(aes(ymin=NMI-sd, ymax=NMI+sd), width=.1,
                  position = position_dodge(width=0.75)) + xlab("Size and rate") + ylab("NMI") + scale_color_discrete("CN difference") +
    ggtitle("NMI vs library size varaince") + theme_bw() + scale_color_discrete("CNVs", labels = c("Amplification", "Duplication"))


}


rmserr <- function(pred, real) {

  res <-  sum((real-pred)**2)
  return(sqrt(res))


}

real_CNV_conv <-  function(sim) {
  res <- sim$cnv %>% select(matches("ploidy[1-9]")) %>%   as.matrix(.) %>% t(.)
  colnames(res) <- colnames(sim$cnv_mat)
  return(res)
}

iterate_rmserr <- function(pred_CNV, real_CNV){

  res <- vector(length = nrow(real_CNV))

  for(i in seq_len(nrow(real_CNV))){

    res[i] <- rmserr(as.numeric(round(pred_CNV[i,])), as.numeric(real_CNV[i,]))

  }

  return(mean(res))

}

test4_core <-  function(parameter,nsegs = NULL) {
  tmp <- run_example_fixed(list('K' = 2, 'spots' = 2, 'changes' = c(-1,1)),
                           list('K' = 2, 'props' = c(0.8,0.2),'theta_shape' = parameter, 'theta_rate' = parameter ), IC = "BIC", steps = 1)
  return(tmp)
}

test2_core <-  function(parameter,nsegs = NULL) {
  tmp <- run_example_fixed(list('K' = 3, 'spots' = 2, 'changes' = c(1,-1)), list('K' = 3, 'props' = c(0.4,0.4,0.2)), IC = "BIC", nsegs = parameter, steps = 1)
  return(tmp)
}


test3_core <-  function(parameter,nsegs = NULL) {
  tmp <- run_example_fixed(list('K' = 2, 'spots' = 2, 'changes' = c(1)), list('K' = 2, 'props' = c(0.7,0.3), 'perc_genes' = parameter), IC = "BIC", nsegs = nsegs, steps = 1)
  return(tmp)
}

test6_core <-  function(parameter,nsegs = NULL) {
  tmp <- run_example_fixed(list('K' = 2, 'spots' = 5, 'changes' = c(1)),
                           list('K' = 2, 'props' = c(1-parameter,parameter)), IC = "BIC", steps = 1)
  return(tmp)
}

test5_core <-  function(parameter,nsegs = NULL) {
  tmp <- run_example_fixed(list('K' = 2, 'spots' = 2, 'changes' = c(1)),
                           list('K' = 2, 'props' = c(0.8,0.2), 'nbinom' = TRUE, "size" = parameter ), IC = "BIC", steps = 1)
  return(tmp)
}

test_cases_simulation <- function(parameters, test_function, inf_res, plot_fun,plot_2=NULL, nsegs = NULL, measure_unit = "size=", test4 = F, gen_plot2=T, fsize = 20, title_word = "overdispersion") {
  pres <- list()
  parameters2 = parameters
  if(length(parameters %% 2 == 0))
    parameters2 = parameters[-1]
  cases <- c(min(parameters), max(parameters))
  titles <-  paste0(c("High ", "Low "), title_word)
  sims <-  list()

  for(i in seq_along(cases)){
    tmp <- test_function(cases[i])
    if(gen_plot2){
      sims[[i]] <- tmp$simulation
    }

    pres[[i]] <- plot_counts.CNVSimulation(tmp$simulation, plot = F, title = titles[i], font_size = fsize)

  }
  if(gen_plot2){
    segmenti_pld <- lapply(sims, function(x) x$cnv %>% select(matches("ploidy")) %>% select(-ploidy_real))
    segmenti_diversi <- lapply(segmenti_pld, function(x) sample(which(apply(x,1, function(x) (var(x) != 0))), size = 1 ))
    if(test4){
      test_ridge <- sapply(sims, FUN = function(x) x$params$theta)

    } else{
      test_ridge <- mapply(sims, segmenti_diversi, FUN = function(x,y) (x$counts[,y] / rowSums(x$counts)) / x$cnv$mu[y])
    }
    test_ridge <-  as.data.frame(test_ridge)
    colnames(test_ridge) <- titles
    test_ridge_long <- reshape2::melt(test_ridge)
    test_ridge_long$clust <- c(sims[[1]]$clust_ids$cluster_id, sims[[2]]$clust_ids$cluster_id)
    print(head(test_ridge_long))
    plot_2 <- ggplot(test_ridge_long, aes(y = variable, fill = paste(clust))) +
      geom_density_ridges2(aes(x = value,  height = ..count..), stat = "density", ) + scale_fill_discrete("clusters") + ylab("") +
      theme_ridges(fsize) +theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

  }
  p1 <- cowplot::plot_grid(pres[[1]], pres[[2]], align = "hv", axis = "lrtb", ncol = 1, labels = c("A", "B"))
  p2 <-  cowplot::plot_grid(plot_2, (plot_fun(inf_res) + theme(text = element_text(size=fsize))), align = "hv", axis = "lrtb", ncol = 2, labels = c("C", "D"))
  res <-cowplot::plot_grid( p1, p2,
    ncol = 1, align = "hv", axis = "lrtb", rel_heights = c(2,1))


  return(res)

}




test_cases_simulation_test1 <- function(distance, clusters ) {

  for(i in seq_along(distance))
  {
    tmp <- run_example(list('K' = j, 'spots' = i + 1, 'changes' = c(-1,1) ), list('K' = j, 'props' = c(0.8,0.2), 'ncells' = 1000))
  }

}