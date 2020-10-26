source("./utils_tests.R")
devtools::load_all("../../annealToolbox/")
MAX_CLUST <-  5
MAX_DIST <- 5
EXAMPLES_FOR_POINT <- 50

set.seed(3)

#props <- list(c(1), c(0.7,0.3), c(0.5,0.3,0.2), c(0.3,0.3,0.2,0.2))


table_clusters <- matrix(NA, nrow = MAX_DIST, ncol = MAX_CLUST)
table_pred <- matrix(NA, nrow = MAX_DIST, ncol = MAX_CLUST)
table_pred_sd <- matrix(NA, nrow = MAX_DIST, ncol = MAX_CLUST)

res <- list()
counter <-  1
for(i in seq_len(MAX_DIST)){
  for(j in seq_len(MAX_CLUST)){
    pred_temp <- vector(length = EXAMPLES_FOR_POINT)
    clusters_temp <-  vector(length = EXAMPLES_FOR_POINT)
      for(k in seq_len(EXAMPLES_FOR_POINT)) {



        props <-  gtools::rdirichlet(1, rep(5,j))
        props <-  round(sort(props, decreasing = T), 1)

        if(!all(props >= 0.1)) {
          props <-  gtools::rdirichlet(1, rep(5,j))
          props <-  round(sort(props, decreasing = T), 1)
        }

        if(sum(props) != 1) {
          rest <- 1 - sum(props)
          props[1] <- props[1] + rest
        }

        print(paste0("Doing cluster ", j, " distance ", i, " replicate ", k))
        try({
        tmp <- run_example_fixed(list('K' = j, 'spots' = i, 'changes' = c(-1,1) ), list('K' = j, 'props' = props, 'ncells' = 1000))
        res[counter] <- list(tmp)
        NMI <-  aricode::NMI(tmp$simulation$clust_ids$cluster_id, tmp$inference$parameters$assignement)
        if(j == 1 & length(tmp$inference$parameters$mixture_weights) == 1)
          {NMI <- 1}
        else if (j == 1)
          {NMI <- sum(tmp$inference$parameters$assignement == tmp$simulation$clust_ids$cluster_id) / length(tmp$simulation$clust_ids$cluster_id)}
        print(tmp$IC)
        print(NMI)

        pred_temp[k] <- NMI
        clusters_temp[k] <- iterate_rmserr(tmp$inference$parameters$cnv_probs, real_CNV_conv(tmp$simulation))

        })
        counter <-  counter + 1
      }

    table_pred[i,j] <- mean(pred_temp, na.rm = T)
    table_pred_sd[i,j] <- sd(pred_temp, na.rm = T)
    table_clusters[i,j] <- mean(clusters_temp, na.rm = T)

  }


}

colnames(table_clusters) <- paste0("Clusters=", seq(1, MAX_CLUST))
rownames(table_clusters) <-  paste0("Distance=", seq(1, MAX_DIST))

colnames(table_pred) <- paste0("Clusters=", seq(1, MAX_CLUST))
rownames(table_pred) <-  paste0("Distance=", seq(1, MAX_DIST))

clusters <- rep(seq(1,MAX_CLUST), rep(MAX_DIST * EXAMPLES_FOR_POINT ,MAX_CLUST))
distances <- rep(seq(1,MAX_DIST), rep(MAX_CLUST * EXAMPLES_FOR_POINT ,MAX_DIST))
examples <- rep(seq(1,EXAMPLES_FOR_POINT), rep(MAX_DIST * MAX_CLUST ,EXAMPLES_FOR_POINT))
names(res) <- paste0("D:", distances, ",C:", clusters,  ",R:",examples)[-200]

to_save <-  list(simulations = res,CNVs = table_pred, clusters = table_clusters, table_pred_sd = table_pred_sd)
save(to_save, file = "test1_BIC_10x.rda")
