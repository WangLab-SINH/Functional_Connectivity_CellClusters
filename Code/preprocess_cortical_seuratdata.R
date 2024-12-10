setwd("F:/brain_data/human/dissections/")
library(Seurat)
library(stringr)



use_dissections <- read.csv("./interdata/cortical_regions.csv")
use_dissections2 <- use_dissections$region
use_dissections2 <- paste0(use_dissections2,".rds")


for (i in 1:length(use_dissections2)) {
  dissection_data <- readRDS(use_dissections2[i])
  diss_region <- substr(use_dissections2[i],1,nchar(use_dissections2[i])-4)
  if(all(dissection_data$dissection == diss_region)){
    print(i)
    mycount <-dissection_data@assays$RNA@data
    meta.features <- dissection_data@assays$RNA@meta.features
    
    genenames <- str_split(meta.features$Gene,"[.]",simplify = T)[,1]
    dup_genes <- unique(genenames[duplicated(genenames)])
    keepgenes <- !(genenames %in% dup_genes)
    
    mycount <- mycount[keepgenes,]
    meta.features <- meta.features[keepgenes,]
    genenames <- genenames[keepgenes]
    
    rownames(mycount) <- genenames
    
    metadata <- dissection_data@meta.data
    reductions <- dissection_data@reductions
    misc <- dissection_data@misc
    
    dissection_data <- CreateSeuratObject(counts = mycount, min.cells = 10)
    dissection_data@meta.data <- cbind(dissection_data@meta.data,metadata)
    dissection_data@reductions <- reductions
    dissection_data@misc <- misc
    dissection_data@assays$RNA@meta.features <- meta.features
    
    
    minGene=500
    maxGene=9000
    pctMT=0.05
    
    dissection_data <- subset(dissection_data, subset = total_genes > minGene & total_genes < maxGene & fraction_mitochondrial < pctMT)
    saveRDS(dissection_data,file = paste0("./MRI_dissections/",use_dissections2[i])) 
  } else{
    print(paste0("No success: ",use_dissections2[i]))
    next;
  }
}


renana <- c("Cerebral cortex (Cx) - Perirhinal gyrus (PRG) - A35-A36.rds",
            "Cerebral cortex (Cx) - Rostral gyrus (RoG) - Dorsal division of MFC - A32.rds")

dissection_data <- readRDS(renana[1])
diss_region <- substr(renana[1],1,nchar(renana[1])-4)
table(dissection_data$dissection)
all(dissection_data$dissection == diss_region)
as.character(dissection_data$dissection[dissection_data$dissection != diss_region][1:10])
dim(dissection_data)
dissection_data <- subset(dissection_data, subset = (dissection == diss_region))
dim(dissection_data)

mycount <-dissection_data@assays$RNA@data
meta.features <- dissection_data@assays$RNA@meta.features
all(rownames(dissection_data) == rownames(meta.features))

genenames <- str_split(meta.features$Gene,"[.]",simplify = T)[,1]
dup_genes <- unique(genenames[duplicated(genenames)])
keepgenes <- !(genenames %in% dup_genes)

mycount <- mycount[keepgenes,]
meta.features <- meta.features[keepgenes,]
genenames <- genenames[keepgenes]

rownames(mycount) <- genenames

metadata <- dissection_data@meta.data
reductions <- dissection_data@reductions
misc <- dissection_data@misc

dissection_data <- CreateSeuratObject(counts = mycount, min.cells = 10)
dissection_data@meta.data <- cbind(dissection_data@meta.data,metadata)
dissection_data@reductions <- reductions
dissection_data@misc <- misc
dissection_data@assays$RNA@meta.features <- meta.features


minGene=500
maxGene=9000
pctMT=0.05

dissection_data <- subset(dissection_data, subset = total_genes > minGene & total_genes < maxGene & fraction_mitochondrial < pctMT)
saveRDS(dissection_data,file = paste0("./MRI_dissections/",renana[1])) 


############################################################################
dissection_data <- readRDS(renana[2])
diss_region <- substr(renana[2],1,nchar(renana[2])-4)
table(dissection_data$dissection)
all(dissection_data$dissection == diss_region)
as.character(dissection_data$dissection[dissection_data$dissection != diss_region][1:10])
dim(dissection_data)
dissection_data <- subset(dissection_data, subset = (dissection == diss_region))
dim(dissection_data)

mycount <-dissection_data@assays$RNA@data
meta.features <- dissection_data@assays$RNA@meta.features
all(rownames(dissection_data) == rownames(meta.features))


genenames <- str_split(meta.features$Gene,"[.]",simplify = T)[,1]
dup_genes <- unique(genenames[duplicated(genenames)])
keepgenes <- !(genenames %in% dup_genes)

mycount <- mycount[keepgenes,]
meta.features <- meta.features[keepgenes,]
genenames <- genenames[keepgenes]

rownames(mycount) <- genenames

metadata <- dissection_data@meta.data
reductions <- dissection_data@reductions
misc <- dissection_data@misc

dissection_data <- CreateSeuratObject(counts = mycount, min.cells = 10)
dissection_data@meta.data <- cbind(dissection_data@meta.data,metadata)
dissection_data@reductions <- reductions
dissection_data@misc <- misc
dissection_data@assays$RNA@meta.features <- meta.features


minGene=500
maxGene=9000
pctMT=0.05

dissection_data <- subset(dissection_data, subset = total_genes > minGene & total_genes < maxGene & fraction_mitochondrial < pctMT)
saveRDS(dissection_data,file = paste0("./MRI_dissections/",renana[2])) 

################################### Proportion of cell clusters in different brain region
region_filenames <- dir("./MRI_dissections/")
meta <- list()
for (i in 1:28) {
  print(i)
  MRI_data <- readRDS(paste0("./MRI_dissections/",region_filenames[i]))
  meta[[i]] <- MRI_data@meta.data
}
names(meta) <- allregionnames
save(meta,file = "./interdata/cortical_regions_metadata.RData")
##
load("./interdata/cortical_regions_metadata.RData")
meta <- lapply(meta,function(x){x$donor_id <- as.character(x$donor_id)
return(x)})
meta <- lapply(meta,function(x){x[x$donor_id %in% c("H18.30.002","H19.30.001","H19.30.002"),]})
cluster_density <- lapply(meta,function(x){table(x$cluster_id)/length(x$cluster_id)})

cluster_density_df <- data.frame()
for (i in 1:length(cluster_density)) {
  temp1 <- data.frame(cluster_density[[i]])
  colnames(temp1) <- c('clusterID','density')
  temp1$type <- names(cluster_density)[i]
  cluster_density_df <- rbind(cluster_density_df,temp1)
  rm(temp1)
}
cluster_density_df <- reshape2::dcast(cluster_density_df,
                                      clusterID~type,value.var = 'density')

rownames(cluster_density_df) <- as.character(cluster_density_df$clusterID)
cluster_density_df <- cluster_density_df[,-1]
cluster_density_df[is.na(cluster_density_df)] <- 0

supercluster_anno <- read.csv("./interdata/cluster_anno.csv")
supercluster_anno <- supercluster_anno[1:461,]
rownames(supercluster_anno) <- as.character(supercluster_anno$Cluster.ID)
clustername <- supercluster_anno[rownames(cluster_density_df),2]
rownames(cluster_density_df) <- clustername
cluster_density_df <- as.data.frame(t(cluster_density_df))

cluster_density_df2 <- cluster_density_df

nonzero_num <- sapply(cluster_density_df2,function(x){sum(x != 0)})
cluster_density_sd <- sapply(cluster_density_df2,sd)
quantile(cluster_density_sd,seq(0,1,0.1))
cluster_density_mean <- sapply(cluster_density_df2,mean)
quantile(cluster_density_mean,seq(0,1,0.1))

keep_cluster <- (nonzero_num >= 14 & cluster_density_sd >= quantile(cluster_density_sd,seq(0,1,0.1))[3] & cluster_density_mean >= quantile(cluster_density_mean,seq(0,1,0.1))[3])
sum(keep_cluster)

cluster_density_df_filtered <- cluster_density_df2[,keep_cluster]



