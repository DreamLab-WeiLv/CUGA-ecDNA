library(Seurat)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(DoubletFinder)

setwd('~/PROJECT/raw_rds')
data<-read.table('samples',header=F,sep='\t')
sclist <- list()
for (j in 1:nrow(data)){
  i<-as.character(data[j,1])
  object <- readRDS(i)
  name <- gsub('.rds','',basename(i))
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
  object <- object %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
  sweep.res.list <- paramSweep_v3(object, PCs = 1:30, sct = F)     
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  p <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  homotypic.prop <- modelHomotypic(object@meta.data$seurat_clusters)
  DoubletRate <- ncol(object)*8*1e-6
  nExp_poi <- round(DoubletRate*ncol(object))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  object <- doubletFinder_v3(object, PCs = 1:30, pN = 0.25, pK = p, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(object@meta.data)[ncol(object@meta.data)] <- "doublet_info"
  object <- subset(object, subset = doublet_info=="Singlet")
  sclist[[i]] <- object
}
seurat_data <- merge(sclist[[1]],unlist(sclist[-1]))

seurat_data <- seurat_data %>% NormalizeData() %>%
  FindVariableFeatures( selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>% 
  RunPCA( npcs = 30,verbose = FALSE) %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.1)

saveRDS(seurat_data,'all.rds')

pan_marker <- c("VWF", "PECAM1", "ENG","CDH5","DCN","ACTA2","COL1A1","COL1A2","CD68", "LYZ", "CD163", "CD14","CD79A", "CD79B", "MS4A1", "MZB1","CD2","CD3D","CD3E","NKG7","GZMA","EPCAM","KRT18","KRT7","KRT19")
p1 <- DimPlot(seurat_data, reduction = "umap", group.by = "orig.ident",raster = T) + theme(legend.position="none")
p2 <- DotPlot(seurat_data, features = pan_marker,group.by = "RNA_snn_res.0.1")+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5))

############################################
##re-clustering of subtype cells  
library(harmony)
seurat_data<- readRDS('B.rds') #
seurat_data <- seurat_data %>%
  NormalizeData() %>%
  FindVariableFeatures( selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 30,verbose = FALSE) %>%
  RunHarmony("orig.ident") %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 1)   # epi resolution=0.5

seurat.markers <- FindAllMarkers(object = seurat_data,
                                 only.pos = T,
                                 min.pct = 0.25,
                                 logfc.threshold = 0.25,
                                 assay = "RNA",slot = "data")

