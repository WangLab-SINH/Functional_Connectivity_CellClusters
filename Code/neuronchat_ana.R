library(Seurat)
library(NeuronChat)
###### Construction of pairwise cell clusters communication network
load("./interdata//cortical_regions_metadata.RData")
common_cluster <- lapply(meta,function(x){unique(as.character(x$cluster_id))})
common_cluster <- Reduce(intersect,common_cluster)

use_region <- names(meta)
filenames <- paste0(use_region,".rds")
all(filenames %in% dir("./data"))

region_pair <- sapply(use_region,function(x){sapply(use_region,function(y){paste(x,y,sep = ",")})})
region_pair <- region_pair[lower.tri(region_pair)]
region_pair <- strsplit(region_pair,",")

region_pair_filename <- lapply(region_pair,function(x){ paste0(x,".rds")})


region_chat <- list()
for (i in 1:378) {
  MRIdata1 <- readRDS(paste0("./data",region_pair_filename[[i]][1]))
  MRIdata1$donor_id <- as.character(MRIdata1$donor_id)
  MRIdata1 <- subset(MRIdata1,subset = (donor_id %in% c("H18.30.002","H19.30.001","H19.30.002")))
  MRIdata1 <- subset(MRIdata1,subset = (cluster_id %in% common_cluster))
  MRIdata1$cluster_id <- paste(region_pair[[i]][1],MRIdata1$cluster_id,sep = "_")
  
  
  MRIdata2 <- readRDS(paste0("./data",region_pair_filename[[i]][2]))
  MRIdata2$donor_id <- as.character(MRIdata2$donor_id)
  MRIdata2 <- subset(MRIdata2,subset = (donor_id %in% c("H18.30.002","H19.30.001","H19.30.002")))
  MRIdata2 <- subset(MRIdata2,subset = (cluster_id %in% common_cluster))
  MRIdata2$cluster_id <- paste(region_pair[[i]][2],MRIdata2$cluster_id,sep = "_")
  
  MRIdata <- merge(MRIdata1,MRIdata2)
  MRIdata <- NormalizeData(MRIdata,assay ="RNA")

  x <- createNeuronChat(MRIdata@assays$RNA@data,DB='human',group.by = as.character(MRIdata@meta.data$cluster_id))
  x <- run_NeuronChat(x,M=100)
  region_chat[[i]] <- net_aggregation(x@net,method = 'weight')
}

region_pair <- sapply(use_region,function(x){sapply(use_region,function(y){paste(x,y,sep = ",")})})
region_pair <- region_pair[lower.tri(region_pair)]

names(region_chat) <- region_pair
save(region_chat,file = "./interdata/net_aggregated_interregion.RData")


################### Intercellular communication between regions and functional connectivity
library(stringr)
library(ggplot2)
library(ggpubr)
library(cowplot)

load("./interdata/net_aggregated_interregion.RData")
asso_cluster <- read.delim("./interdata/degree_cluster_netdata.txt")
pos_asso_cluster <- asso_cluster$Cluster[asso_cluster$association > 0]
asso_clusterid <- strsplit(pos_asso_cluster,"_")
asso_clusterid <- sapply(asso_clusterid,function(x){x[2]})
asso_clusterid <- intersect(asso_clusterid,use_cluster)

region_chat_assoclusters <- matrix(nrow = 28,ncol = 28)

for (i in 1:27) {
  for (j in (i + 1):28) {
    temp1 <- region_chat[[region_pair[i,j]]]
    cur_regions <- strsplit(region_pair[i,j],",")
    region1 <- cur_regions[[1]][1]
    region2 <- cur_regions[[1]][2]
    chat1 <- temp1[str_detect(rownames(temp1),region1),str_detect(colnames(temp1),region2)]
    chat2 <- temp1[str_detect(rownames(temp1),region2),str_detect(colnames(temp1),region1)]
    chat2 <- t(chat2)
    final_chat <- (chat1 + chat2)
    use_clusters1 <- paste(region1,asso_clusterid,sep = "_")
    use_clusters2 <- paste(region2,asso_clusterid,sep = "_")
    use_clusters <- c(use_clusters1,use_clusters2)
    final_chat <- final_chat[intersect(rownames(final_chat),use_clusters),intersect(colnames(final_chat),use_clusters)]
    region_chat_assoclusters[i,j] <- mean(final_chat)
  }
}

rownames(region_chat_assoclusters) <- use_region
colnames(region_chat_assoclusters) <- use_region

region_chat_assoclusters1 <- t(region_chat_assoclusters)
region_chat_assoclusters[lower.tri(region_chat_assoclusters)] <- region_chat_assoclusters1[lower.tri(region_chat_assoclusters1)]

isSymmetric(region_chat_assoclusters)
diag(region_chat_assoclusters) <- 0


chat_degree_assoclusters <- apply(region_chat_assoclusters,1,sum)
chat_degree_assoclusters <- chat_degree_assoclusters[rownames(connectome_degree)]
all(rownames(connectome_degree) == names(chat_degree_assoclusters))


plotdata <- data.frame(connectome_degree = connectome_degree$cutoff0.3,chat_degree_assoclusters= chat_degree_assoclusters)
plotdata <- plotdata[plotdata$chat_degree_assoclusters > 11,]
p1 <- ggplot(plotdata,aes(chat_degree_assoclusters,connectome_degree)) +
  geom_point(color = "#a48cbe") +theme_classic() +
  stat_cor(label.x = min(plotdata$chat_degree_assoclusters), label.y = max(plotdata$connectome_degree), size = 5) +
  labs(y="Degree of functional network",x="Degree of region communication network")+
  theme(axis.title.y = element_text(size = 14),
        axis.title.x =element_text(size = 14) )+
  theme(axis.text.x = element_text(size = 12,color="black"),
        axis.text.y = element_text(size = 12,color="black"),
        plot.title = element_text(hjust = 0.5,size=16))


all(rownames(region_chat_assoclusters) == rownames(connect_mat))
all(colnames(region_chat_assoclusters) == colnames(connect_mat))


plotdata2 <- data.frame(functional_connectome = connect_mat[lower.tri(connect_mat)],region_chat = region_chat_assoclusters[lower.tri(region_chat_assoclusters)])

p2 <- ggplot(plotdata2,aes(region_chat,functional_connectome)) +
  geom_point(color = "#a48cbe") +theme_classic() +
  stat_cor(label.x = min(plotdata2$region_chat), label.y = max(plotdata2$functional_connectome), size = 5) +
  labs(y="Functional network",x="Region communication network")+
  theme(axis.title.y = element_text(size = 14),
        axis.title.x =element_text(size = 14) )+
  theme(axis.text.x = element_text(size = 12,color="black"),
        axis.text.y = element_text(size = 12,color="black"),
        plot.title = element_text(hjust = 0.5,size=16))
plot_grid(p1,p2,nrow = 1)


##

neg_asso_cluster <- asso_cluster$Cluster[asso_cluster$association < 0]
asso_clusterid <- strsplit(neg_asso_cluster,"_")
asso_clusterid <- sapply(asso_clusterid,function(x){x[2]})
asso_clusterid <- intersect(asso_clusterid,use_cluster)

region_chat_assoclusters <- matrix(nrow = 28,ncol = 28)

for (i in 1:27) {
  for (j in (i + 1):28) {
    temp1 <- region_chat[[region_pair[i,j]]]
    cur_regions <- strsplit(region_pair[i,j],",")
    region1 <- cur_regions[[1]][1]
    region2 <- cur_regions[[1]][2]
    chat1 <- temp1[str_detect(rownames(temp1),region1),str_detect(colnames(temp1),region2)]
    chat2 <- temp1[str_detect(rownames(temp1),region2),str_detect(colnames(temp1),region1)]
    chat2 <- t(chat2)
    final_chat <- (chat1 + chat2)
    use_clusters1 <- paste(region1,asso_clusterid,sep = "_")
    use_clusters2 <- paste(region2,asso_clusterid,sep = "_")
    use_clusters <- c(use_clusters1,use_clusters2)
    final_chat <- final_chat[intersect(rownames(final_chat),use_clusters),intersect(colnames(final_chat),use_clusters)]
    region_chat_assoclusters[i,j] <- mean(final_chat)
  }
}

rownames(region_chat_assoclusters) <- use_region
colnames(region_chat_assoclusters) <- use_region

region_chat_assoclusters1 <- t(region_chat_assoclusters)
region_chat_assoclusters[lower.tri(region_chat_assoclusters)] <- region_chat_assoclusters1[lower.tri(region_chat_assoclusters1)]

isSymmetric(region_chat_assoclusters)
diag(region_chat_assoclusters) <- 0


chat_degree_assoclusters <- apply(region_chat_assoclusters,1,sum)

chat_degree_assoclusters <- chat_degree_assoclusters[rownames(connectome_degree)]
all(rownames(connectome_degree) == names(chat_degree_assoclusters))

plotdata <- data.frame(connectome_degree = connectome_degree$cutoff0.3,chat_degree_assoclusters= chat_degree_assoclusters)

p3 <- ggplot(plotdata,aes(chat_degree_assoclusters,connectome_degree)) +
  geom_point(color = "#a9dce6") +theme_classic() +
  stat_cor(label.x = min(plotdata$chat_degree_assoclusters), label.y = max(plotdata$connectome_degree), size = 5) +
  labs(y="Degree of functional network",x="Degree of region communication network")+
  theme(axis.title.y = element_text(size = 14),
        axis.title.x =element_text(size = 14) )+
  theme(axis.text.x = element_text(size = 12,color="black"),
        axis.text.y = element_text(size = 12,color="black"),
        plot.title = element_text(hjust = 0.5,size=16))


all(rownames(region_chat_assoclusters) == rownames(connect_mat))
all(colnames(region_chat_assoclusters) == colnames(connect_mat))


cor.test(region_chat_assoclusters[lower.tri(region_chat_assoclusters)],connect_mat[lower.tri(connect_mat)])
cor.test(region_chat_assoclusters[lower.tri(region_chat_assoclusters)],connect_mat[lower.tri(connect_mat)],method = "spearman")

plotdata2 <- data.frame(functional_connectome = connect_mat[lower.tri(connect_mat)],region_chat = region_chat_assoclusters[lower.tri(region_chat_assoclusters)])

p4 <- ggplot(plotdata2,aes(region_chat,functional_connectome)) +
  geom_point(color = "#a9dce6") +theme_classic() +
  stat_cor(label.x = min(plotdata2$region_chat), label.y = max(plotdata2$functional_connectome), size = 5) +
  labs(y="functional network",x="Region communication network")+
  theme(axis.title.y = element_text(size = 14),
        axis.title.x =element_text(size = 14) )+
  theme(axis.text.x = element_text(size = 12,color="black"),
        axis.text.y = element_text(size = 12,color="black"),
        plot.title = element_text(hjust = 0.5,size=16))
plot_grid(p3,p4,nrow = 1)

##
all_cluster <- read.csv("./interdata/cluster_density_df.csv",row.names = 1)
all_cluster <- colnames(all_cluster)
other_cluster <- setdiff(all_cluster,asso_cluster$Cluster)

asso_clusterid <- strsplit(other_cluster,"_")
asso_clusterid <- sapply(asso_clusterid,function(x){x[2]})
asso_clusterid <- intersect(asso_clusterid,use_cluster)

region_chat_assoclusters <- matrix(nrow = 28,ncol = 28)

for (i in 1:27) {
  for (j in (i + 1):28) {
    temp1 <- region_chat[[region_pair[i,j]]]
    cur_regions <- strsplit(region_pair[i,j],",")
    region1 <- cur_regions[[1]][1]
    region2 <- cur_regions[[1]][2]
    chat1 <- temp1[str_detect(rownames(temp1),region1),str_detect(colnames(temp1),region2)]
    chat2 <- temp1[str_detect(rownames(temp1),region2),str_detect(colnames(temp1),region1)]
    chat2 <- t(chat2)
    final_chat <- (chat1 + chat2)
    use_clusters1 <- paste(region1,asso_clusterid,sep = "_")
    use_clusters2 <- paste(region2,asso_clusterid,sep = "_")
    use_clusters <- c(use_clusters1,use_clusters2)
    final_chat <- final_chat[intersect(rownames(final_chat),use_clusters),intersect(colnames(final_chat),use_clusters)]
    region_chat_assoclusters[i,j] <- mean(final_chat)
  }
}

rownames(region_chat_assoclusters) <- use_region
colnames(region_chat_assoclusters) <- use_region

region_chat_assoclusters1 <- t(region_chat_assoclusters)
region_chat_assoclusters[lower.tri(region_chat_assoclusters)] <- region_chat_assoclusters1[lower.tri(region_chat_assoclusters1)]

isSymmetric(region_chat_assoclusters)
diag(region_chat_assoclusters) <- 0


chat_degree_assoclusters <- apply(region_chat_assoclusters,1,sum)

chat_degree_assoclusters <- chat_degree_assoclusters[rownames(connectome_degree)]
all(rownames(connectome_degree) == names(chat_degree_assoclusters))

cor(chat_degree_assoclusters,connectome_degree$cutoff0.3)
cor.test(chat_degree_assoclusters,connectome_degree$cutoff0.3)
cor(chat_degree_assoclusters,connectome_degree)
cor(chat_degree_assoclusters,connectome_degree,method = "spearman")

plotdata <- data.frame(connectome_degree = connectome_degree$cutoff0.3,chat_degree_assoclusters= chat_degree_assoclusters)

p5 <- ggplot(plotdata,aes(chat_degree_assoclusters,connectome_degree)) +
  geom_point(color = "#e2b159") +theme_classic() +
  stat_cor(label.x = min(plotdata$chat_degree_assoclusters), label.y = max(plotdata$connectome_degree), size = 5) +
  labs(y="Degree of functional network",x="Degree of region communication network")+
  theme(axis.title.y = element_text(size = 14),
        axis.title.x =element_text(size = 14) )+
  theme(axis.text.x = element_text(size = 12,color="black"),
        axis.text.y = element_text(size = 12,color="black"),
        plot.title = element_text(hjust = 0.5,size=16))


all(rownames(region_chat_assoclusters) == rownames(connect_mat))
all(colnames(region_chat_assoclusters) == colnames(connect_mat))


cor.test(region_chat_assoclusters[lower.tri(region_chat_assoclusters)],connect_mat[lower.tri(connect_mat)])
cor.test(region_chat_assoclusters[lower.tri(region_chat_assoclusters)],connect_mat[lower.tri(connect_mat)],method = "spearman")

plotdata2 <- data.frame(functional_connectome = connect_mat[lower.tri(connect_mat)],region_chat = region_chat_assoclusters[lower.tri(region_chat_assoclusters)])

p6 <- ggplot(plotdata2,aes(region_chat,functional_connectome)) +
  geom_point(color = "#e2b159") +theme_classic() +
  stat_cor(label.x = min(plotdata2$region_chat), label.y = max(plotdata2$functional_connectome), size = 5) +
  labs(y="functional network",x="Region communication network")+
  theme(axis.title.y = element_text(size = 14),
        axis.title.x =element_text(size = 14) )+
  theme(axis.text.x = element_text(size = 12,color="black"),
        axis.text.y = element_text(size = 12,color="black"),
        plot.title = element_text(hjust = 0.5,size=16))
plot_grid(p5,p6,nrow = 1)


################### cell clusters communication within regions
load("./interdata/cortex_dissection.RData")

cortex_filename <- paste0(cortex_dissection,".rds")
for (i in 1:28) {
  print(i)
  MRI_data <- readRDS(cortex_filename[i])
  MRI_data <- NormalizeData(MRI_data,assay ="RNA")
  MRI_data <- subset(MRI_data,subset = (donor_id %in% c("H18.30.002","H19.30.001","H19.30.002")))
  x <- createNeuronChat(MRI_data@assays$RNA@data,DB='human',group.by = as.character(MRI_data@meta.data$cluster_id))
  x <- run_NeuronChat(x,M=100)
  net_aggregated_x <- net_aggregation(x@net,method = 'weight')
  saveRDS(net_aggregated_x,file = paste0("./neuronchat/","net_aggregated_cluster_",cortex_filename[i]))
}


allfiles <- dir("./neuronchat")
use_files <- allfiles[str_detect(allfiles,"_cluster_")]

cluster_chat <- list()
for (i in 1:28) {
  cluster_chat[[i]] <- readRDS(use_files[[i]])
}

regions <- strsplit(use_files,"_")
regions <- sapply(regions,function(x){x[4]})
regions <- strsplit(regions,"[.]")
regions <- sapply(regions,function(x){x[1]})

names(cluster_chat) <- regions

degree_cluster_mat <- read.delim("./interdata/degree_cluster_netdata.txt")
clusterid <- strsplit(degree_cluster_mat$Cluster,"_")
clusterid <- sapply(clusterid,function(x){x[2]})

clusterid_pos <- clusterid[degree_cluster_mat$association > 0]
clusterid_neg <- clusterid[degree_cluster_mat$association < 0]

comparer_data <- lapply(cluster_chat,function(x){
  allcluster <- rownames(x)
  use_cluster1 <- intersect(clusterid_pos,allcluster)
  use_cluster2 <- intersect(clusterid_neg,allcluster)
  chatmat1 <- x[use_cluster1,use_cluster1]
  chatmat2 <- x[use_cluster2,use_cluster2]
  chatmat3 <- x[setdiff(allcluster,clusterid),setdiff(allcluster,clusterid)]
  out <- list(pos_asso_cluster = chatmat1[lower.tri(chatmat1)],
              neg_asso_cluster = chatmat2[lower.tri(chatmat2)],
              other_cluster = chatmat3[lower.tri(chatmat3)])
  return(out)
})


plotdata <- lapply(comparer_data,function(x){data.frame(chat_weight = c(x[[1]],x[[2]],x[[3]]),
                                                        gr = c(rep("Pos_asso_cluster",length(x[[1]])),
                                                               rep("Neg_asso_cluster",length(x[[2]])),
                                                               rep("Other",length(x[[3]]))))})

my_comparisons <- list(c("Pos_asso_cluster", "Neg_asso_cluster"), c("Pos_asso_cluster", "Other"), c("Neg_asso_cluster", "Other"))


plot_list <- list()
for (i in 1:28) {
  plot_list[[i]] <- ggplot(plotdata[[i]], aes(x = gr, y = chat_weight, fill = gr)) +
    geom_boxplot(fill = c("#a9dce6","#e2b159","#a48cbe"),outlier.shape = NA) + 
    stat_compare_means(label = "p.signif",comparisons = my_comparisons) + theme_classic() +
    labs(y="Cell clusters chat weight",x = NULL,title = regions[i])+
    theme(axis.title.y = element_text(size = 10,color="black"))+
    theme(axis.text.x = element_text(size = 8,color="black",angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8,color="black"),
          plot.title = element_text(hjust = 0.5,size=12,color="black")) +
    theme(legend.position = "none")
}
p1 <- plot_grid(plotlist = plot_list,nrow = 4)


plotdata2 <- data.frame(
  cluster_type = rep(c("pos_asso_cluster","neg_asso_cluster","other_cluster"),each = 28),
  region_chat_weight = c(sapply(comparer_data,function(x){median(x[["pos_asso_cluster"]])}),
                         sapply(comparer_data,function(x){median(x[["neg_asso_cluster"]])}),
                         sapply(comparer_data,function(x){median(x[["other_cluster"]])})),
  region = rep(names(comparer_data),3)
)
plotdata2$cluster_type <- factor(plotdata2$cluster_type,levels = c("neg_asso_cluster","pos_asso_cluster","other_cluster"))

plotdata2$paired <- plotdata2$region
df_p_val1 <- plotdata2 %>%
  t_test(region_chat_weight  ~ cluster_type,paired = TRUE,p.adjust.method = "fdr") %>%
  add_significance(p.col = "p") %>% 
  add_xy_position(x = "group", dodge = 0.8) 



colors <-  c("#64cccf","#7ee7bb","#f7905a","#e187cb","#fb948d",
             "#e2b159","#b2db87","#a9dce6","#a48cbe","#ea5c6f",
             "#1a476f","#90353b","#55752f","#e37e00","#6e8e84",
             "#c10534","#938dd2","#cac27e","#a0522d","#7b92a8",
             "#2d6d66","#9c8847","#bfa19c","#ffd200","#d9e6eb",
             "#F28E2B","#e187cb","#59A14F", "#EDC948", "#B07AA1", "#FF9DA7")
p2 <- ggplot(plotdata2,aes(x=cluster_type,y=region_chat_weight)) +
  stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.1)+
  geom_boxplot(position=position_dodge(width =0.2),width=0.4,outlier.shape = NA,fill =c("#AECDE1","#E69F84","#7DBFA7"))+
  geom_line(aes(group=region),position = position_dodge(0.2),color="grey80") +
  geom_point(aes(fill = region),pch=21,
             position = position_dodge(0.2),size = 1.5) + ylim(0,1.1) +
  stat_pvalue_manual(df_p_val1,label = "p.signif",label.size=4,hide.ns = F)+
  scale_fill_manual(values = colors)+
  scale_x_discrete(guide = "prism_bracket")+
  labs(x=NULL,y=NULL)+
  theme_prism(base_line_size =0.5)+
  theme(axis.line = element_line(color = "black",size = 0.4),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y = element_text(color="black",size=10),
        axis.text.x = element_text(margin = margin(t = -5),color="black",size=10),
        panel.grid = element_blank()) +
  NoLegend()
