source("tests/utils_tests.R")
devtools::load_all(".")

K <- c(1,2)
PGENES <- seq(0.1,1, by = 0.1)
REP <- 10

props <- c(0.7,0.3)


table_clusters <- matrix(NA, nrow = length(PGENES), ncol = length(K))
table_pred <- matrix(NA, nrow = length(PGENES), ncol = length(K))
table_pred_sd <- matrix(NA, nrow = length(PGENES), ncol = length(K))

for( j in seq_along(K)){
  for (i in seq_along(PGENES)){
    pred_temp <- vector(length = REP)
    clusters_temp <-  vector(length = REP)
    for (k in seq_len(REP)){



      print(paste0("Doing cluster ", 2, " %genes ", PGENES[i], " replicate ", k))
      print(props)
      try({
      tmp <- run_example_fixed(list('K' = 2, 'spots' = 2, 'changes' = c(j)), list('K' = 2, 'props' = props, 'perc_genes' = PGENES[i]), IC = "BIC", nsegs = 2)

        print(tmp$inference$parameters$assignement)
        NMI <-  aricode::NMI(tmp$simulation$clust_ids$cluster_id, tmp$inference$parameters$assignement)

        print(tmp$IC)
        print(NMI)

        pred_temp[k] <- NMI
        clusters_temp[k] <- iterate_rmserr(tmp$inference$parameters$cnv_probs, real_CNV_conv(tmp$simulation))

      })
      table_pred[i,j] <- mean(pred_temp, na.rm = T)
      table_pred_sd[i,j] <- sd(pred_temp, na.rm = T)
      table_clusters[i,j] <- mean(clusters_temp, na.rm = T)



    }
  }
}


colnames(table_clusters) <- paste0("Changes=", K)
rownames(table_clusters) <-  paste0("", PGENES)

colnames(table_pred) <- paste0("Changes=", K)
rownames(table_pred) <-  paste0("", PGENES)


to_save <-  list(CNVs = table_pred, clusters = table_clusters, table_pred_sd = table_pred_sd)
save(to_save, file = "test3_BIC.rda")

