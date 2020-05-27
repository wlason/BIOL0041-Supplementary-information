library(pheatmap)

load("Fisher_data.RData")
load("hmap.genes.cancer.RData")
load("robust_1000.RData")

rnaseq_annot$Cluster_ID <- NULL

for (k in 4:6){
hc <- clust_robust[[k]][["consensusMatrix"]]
  
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

for (cluster in unique(clusters$Cluster_ID)){
  
  for (gene in cosmic_list){
    fisher_data <- rnaseq_annot_cancer[which(rnaseq_annot_cancer$Sample_type ==  "Cancer"), c(gene, "Cluster_ID")]
    
    fisher_matrix <- matrix(c(length(which(fisher_data[, gene] == 1 & fisher_data[, "Cluster_ID"] == cluster)),
                              length(which(fisher_data[, gene] == 0 & fisher_data[, "Cluster_ID"] == cluster)),
                              length(which(fisher_data[, gene] == 1 & fisher_data[, "Cluster_ID"] != cluster)),
                              length(which(fisher_data[, gene] == 0 & fisher_data[, "Cluster_ID"] != cluster))), nrow = 2)
    rownames(fisher_matrix) <- c("c_interest", "rest")
    colnames(fisher_matrix) <- c("MU", "WT")
    
    f <- fisher.test(fisher_matrix)
    
    if (gene == cosmic_list[1] && cluster == clusters[1]){
      l <- length(unique(clusters$Cluster_ID))
      results <- vector(mode = "list", length = l)
      for (i in 1:l){
        results[[i]] <- data.frame()
      }
    }
    
    results[[cluster]][gene, "p_val"] <- f$p.value
    results[[cluster]][gene, "odds"] <- f$estimate
    results[[cluster]][gene, "WT"] <- fisher_matrix["c_interest", "WT"]
    results[[cluster]][gene, "MU"] <- fisher_matrix["c_interest", "MU"]
    results[[cluster]][gene, "gene"] <- gene
    results[[cluster]][gene, "cluster"] <- cluster
  }
}

results_x <- bind_rows(results)
results_x$p_adj <- p.adjust(results_x$p_val, method = "BH")
results_x <- subset(results_x, MU>WT & p_adj<0.05)
results_x$p_adj_star <- stars.pval(results_x$p_adj)

write.csv(results_x, file = paste0("results_", k, ".RData"))
}
