library(Seurat)
region_filenames <- dir("./MRI_dissections/")
cluster_pseudobulk <- list()
for (i in 1:28) {
  print(i)
  MRI_data <- readRDS(paste0("./MRI_dissections/",region_filenames[i]))
  MRI_data$donor_id <- as.character(MRI_data$donor_id)
  MRI_data <- subset(MRI_data,subset = (donor_id %in% c("H18.30.002","H19.30.001","H19.30.002")))
  cluster_pseudobulk[[i]] <- Seurat:::PseudobulkExpression(object = MRI_data,group.by=c("cluster_id"),
                                                          pb.method = 'aggregate', slot = 'counts')$RNA
}
save(cluster_pseudobulk,file = "./interata/cluster_pseudobulk.RData") 

#####
load("./interata/cluster_pseudobulk.RData")
use_cluster <- lapply(cluster_pseudobulk,colnames)
use_cluster <- Reduce(intersect,use_cluster)
cluster_pseudobulk <- lapply(cluster_pseudobulk,function(x){x[,use_cluster]})

cluster_pseudobulk2 <- list()
for (cluster in use_cluster) {
  cluster_pseudobulk2[[cluster]] <- sapply(cluster_pseudobulk,function(x){x[,cluster]})
}

cluster_pseudobulk2 <- lapply(cluster_pseudobulk2,function(x){colnames(x) <- cortical_regions
return(x)})


supercluster_anno <- read.csv("./interdata/cluster_anno.csv")
supercluster_anno <- supercluster_anno[1:461,]
supercluster_anno$Cluster.ID <- as.character(supercluster_anno$Cluster.ID)
rownames(supercluster_anno) <- supercluster_anno$Cluster.ID

all(names(cluster_pseudobulk2) %in% rownames(supercluster_anno))
supercluster_anno2 <- supercluster_anno[names(cluster_pseudobulk2),]
any(duplicated(supercluster_anno2$Cluster.name))
all(rownames(supercluster_anno2) == names(cluster_pseudobulk2))
names(cluster_pseudobulk2) <- supercluster_anno2$Cluster.name


asso_cluster <- read.delim("./interdata/degree_cluster_netdata.txt")
research_cluster <- asso_cluster$Cluster
research_cluster <- intersect(research_cluster,names(cluster_pseudobulk2))

cluster_pseudobulk2 <- cluster_pseudobulk2[research_cluster]

filter_genes <- function(df){
  nonzero_num <- sapply(df,function(x){sum(x > 0)})
  keep_genes <- (nonzero_num >= 3)
  df <- df[,keep_genes]
  return(colnames(df))
}
cluster_pseudobulk2 <- lapply(cluster_pseudobulk2,t)
cluster_pseudobulk2 <- lapply(cluster_pseudobulk2,as.data.frame)
filtered_genes <- lapply(cluster_pseudobulk2,filter_genes)
keep_genes <- Reduce(intersect,filtered_genes)
length(keep_genes)


cluster_pseudobulk2 <- lapply(cluster_pseudobulk2,function(x){x[,keep_genes]})
lapply(cluster_pseudobulk2,dim)

library(edgeR)

cluster_pseudobulk2 <- lapply(cluster_pseudobulk2,function(x){t(x)})



cluster_pseudobulk3 <- lapply(cluster_pseudobulk2,function(x){
  y <- DGEList(counts = x)
  y <- calcNormFactors(y, method = "TMM")
  normalized_counts <- cpm(y, log = TRUE)
  return(normalized_counts)
})



cluster_pseudobulk3 <- lapply(cluster_pseudobulk3,t)

connectome_degree <- read.csv("./interdata/degree_centrality.csv",row.names = 1)
all(sapply(cluster_pseudobulk3,function(x){all(rownames(x) == rownames(connectome_degree))}))

cluster_pseudobulk3 <- lapply(cluster_pseudobulk3,function(x){x[rownames(connectome_degree),]})
all(sapply(cluster_pseudobulk3,function(x){all(rownames(x) == rownames(connectome_degree))}))

filenames <- names(cluster_pseudobulk3)
for (i in 1:length(filenames)) {
  print(i)
  temp <- cluster_pseudobulk3[[i]]
  write.csv(temp,file = paste0("./data/",filenames[i],".csv"),quote = FALSE)
}


