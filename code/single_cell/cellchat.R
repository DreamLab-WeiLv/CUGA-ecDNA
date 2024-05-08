library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(CellChat)

seurat_data<-readRDS('all.rds')
d<- subset(seurat_data,ecDNA == 'Positive')
#d<- subset(seurat_data,ecDNA == 'Negative')

data.input = d@assays$RNA@data 
meta = d@meta.data
meta$sub_cell_type <- factor(meta$sub_cell_type)
meta$sub_cell_type = droplevels(meta$sub_cell_type, exclude = setdiff(levels(meta$sub_cell_type),unique(meta$sub_cell_type)))

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "sub_cell_type")
cellchat <- setIdent(cellchat, ident.use = "sub_cell_type")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat,raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)  
write.csv(df.net, "cell-cell_communications.all.csv")

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

save(cellchat, file = "cellchat.RData")
#save(cellchat_nonec, file = "cellchat_nonec.RData")

load('cellchat.RData')
cellchat.ec <- cellchat
cellchat.ec <- netAnalysis_computeCentrality(cellchat.ec, slot.name = "netP")
load('cellchat_nonec.RData')
cellchat.nonec <- cellchat_nonec
cellchat.nonec <- netAnalysis_computeCentrality(cellchat.nonec, slot.name = "netP")
object.list <- list(ec = cellchat.ec, nonec = cellchat.nonec)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
save(cellchat, file = 'cellchat_merge.RData')

df.net <- subsetCommunication(cellchat)
pairLR.use <- as.data.frame(c('HLA-A_CD8A','HLA-A_CD8B','HLA-B_CD8A','HLA-B_CD8B','HLA-C_CD8A','HLA-C_CD8B',
                              'HLA-E_CD8A','HLA-E_CD8B','HLA-F_CD8A','HLA-F_CD8B'))
colnames(pairLR.use) <- 'interaction_name'

p <- netVisual_bubble(cellchat,
                      sources.use = c("Tumor"),
                      targets.use = c('CD8','INF+T','NK-1'),
                      comparison = c(2, 1),
                      pairLR.use = pairLR.use)
ggsave('netVisual_bubble.pdf',p,w=3.5,h=3)


