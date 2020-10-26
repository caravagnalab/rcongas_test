source("./utils_tests.R")
devtools::load_all("../../annealToolbox/")

K <- 1
size <- c(1,2,4,8,16,32,50,100,150,200)
REP <- 50

props <- c(0.7,0.3)


table_clusters <- matrix(NA, nrow = length(size), ncol = length(K))
table_pred <- matrix(NA, nrow = length(size), ncol = length(K))
table_pred_sd <- matrix(NA, nrow = length(size), ncol = length(K))

for( j in seq_along(K)){
  for (i in seq_along(size)){
    pred_temp <- vector(length = REP)
    clusters_temp <-  vector(length = REP)
    for (k in seq_len(REP)){



      print(paste0("Doing cluster ", 2, " size ", size[i], " replicate ", k))
      print(props)
      try({
        tmp <- run_example_fixed(list('K' = 2, 'spots' = 3, 'changes' = c(j)),
                                 list('K' = 2, 'props' = props, 'nbinom' = TRUE, "size" = size[i] ), IC = "BIC")

        NMI <-  aricode::NMI(tmp$simulation$clust_ids$cluster_id, tmp$inference$parameters$assignement)

        print(NMI)
        print(tmp$IC)

        pred_temp[k] <- NMI
        clusters_temp[k] <- iterate_rmserr(tmp$inference$parameters$cnv_probs, real_CNV_conv(tmp$simulation))

      })

      table_pred[i,j] <- mean(pred_temp, na.rm = T)
      table_pred_sd[i,j] <- sd(pred_temp, na.rm = T)
      table_clusters[i,j] <- mean(clusters_temp, na.rm = T)



    }
  }
}


colnames(table_clusters) <- paste0("Clusters=", K)
rownames(table_clusters) <-  paste0("", size)

colnames(table_pred) <- paste0("Clusters=", K)
rownames(table_pred) <-  paste0("", size)


to_save <-  list(CNVs = table_pred, clusters = table_clusters, table_pred_sd = table_pred_sd)
save(to_save, file = "test5_BIC_10x.rda")
