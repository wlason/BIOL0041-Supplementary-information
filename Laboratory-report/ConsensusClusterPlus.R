library(ConsensusClusterPlus)
load("hmap.genes.cancer.RData")
clust_robust <- ConsensusClusterPlus(as.matrix(hmap_genes_cancer),
                                     maxK=10,
                                     reps=1000,
                                     pItem=0.80,
                                     pFeature=1.00,
                                     title="Pancancer robust 1000", plot="png",
                                     clusterAlg="hc",
                                     distance="euclidean",
                                     finalLinkage="complete", innerLinkage="complete",
                                     seed=2019,
                                     verbose = TRUE)
save(clust_robust, file = "robust_1000.RData")