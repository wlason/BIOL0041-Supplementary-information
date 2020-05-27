library(pheatmap)

load("Fisher_data.RData")
load("hmap.genes.cancer.RData")
load("robust_1000.RData")



for (k in 4:6){
hc <- clust_robust[[k]][["consensusMatrix"]]

rnaseq_annot$Cluster_ID <- NULL  
clusters <- as.data.frame(clust_robust[[k]][["consensusClass"]])
colnames(clusters) <- "Cluster_ID"
rnaseq_annot <- merge(x = rnaseq_annot, y = clusters, by = "row.names")
rownames(rnaseq_annot) <- rnaseq_annot$Row.names
rnaseq_annot$Row.names <- NULL
rnaseq_annot$Cluster_ID <- as.character(rnaseq_annot$Cluster_ID)

for(ext in c(".pdf", ".png")){
pheatmap(hmap_genes_cancer, 
         show_colnames = FALSE, 
         filename = paste0("ccplus_cancer_with_clusters_k_", k, ext), 
         annotation_col = rnaseq_annot_cancer[, c(top_drivers, "Cluster_ID")],
         annotation_row = gene_annot,
         clustering_distance_cols = as.dist(hc),
         annotation_names_row = F,
         treeheight_col = 100, 
         treeheight_row = 0,
         cutree_cols = k,
         fontsize_row = 6,
         annotation_colors = c(mut_colours, prog_colours))
}

cluster_ids <- c(sort(unique(rnaseq_annot$Cluster_ID)), "pan")    #Pan is a placeholder for "all clusters"
combinations <- combn(cluster_ids, m = 2, simplify = FALSE)   #Finding all the cluster combinations that need to be tested
for (i in 1:length(combinations)){     
  
  cluster_1 <- combinations[[i]][1]
  cluster_2 <- combinations[[i]][2]
  
  for (gene in cosmic_list){
  
    if (cluster_2 == "pan"){

      fisher_data <- rnaseq_annot[which(rnaseq_annot$Sample_type ==  "Cancer"), c(gene, "Cluster_ID")]
      
      fisher_matrix <- matrix(c(length(which(fisher_data[, gene] == 1 & fisher_data[, "Cluster_ID"] == cluster_1)),
                                length(which(fisher_data[, gene] == 0 & fisher_data[, "Cluster_ID"] == cluster_1)),
                                length(which(fisher_data[, gene] == 1 & fisher_data[, "Cluster_ID"] != cluster_1)),
                                length(which(fisher_data[, gene] == 0 & fisher_data[, "Cluster_ID"] != cluster_1))), nrow = 2)
      
      rownames(fisher_matrix) <- c("cluster_1", "rest")
      colnames(fisher_matrix) <- c("MU", "WT")
      
      f <- fisher.test(fisher_matrix)

    } else {

      fisher_data <- rnaseq_annot[which(rnaseq_annot$Sample_type ==  "Cancer"), c(gene, "Cluster_ID")]
    
      fisher_matrix <- matrix(c(length(which(fisher_data[, gene] == 1 & fisher_data[, "Cluster_ID"] == cluster_1)),
                                length(which(fisher_data[, gene] == 0 & fisher_data[, "Cluster_ID"] == cluster_1)),
                                length(which(fisher_data[, gene] == 1 & fisher_data[, "Cluster_ID"] == cluster_2)),
                                length(which(fisher_data[, gene] == 0 & fisher_data[, "Cluster_ID"] == cluster_2))), nrow = 2)
      rownames(fisher_matrix) <- c("cluster_1", "cluster_2")
      colnames(fisher_matrix) <- c("MU", "WT")
    
      f <- fisher.test(fisher_matrix)

    }
    
    if (gene == cosmic_list[1] && i == 1){
        results <- data.frame()
        n <- 1
      }

    results[n, "gene"] <- gene
    results[n, "cluster_1"] <- cluster_1
    results[n, "cluster_2"] <- cluster_2
    results[n, "c1_WT"] <- fisher_matrix["cluster_1", "WT"]
    results[n, "c1_MUT"] <- fisher_matrix["cluster_1", "MU"]
    results[n, "odds"] <- f$estimate
    results[n, "p_val"] <- f$p.value
    
    n <- n + 1
    
  }
}

results$p_adj <- p.adjust(results$p_val, method = "BH")
results_x <- subset(results, c1_MUT>c1_WT & p_adj<0.05)
results_x$p_adj_star <- stars.pval(results_x$p_adj)

write.csv(results_x, file = paste0("results_", k, ".RData"))

