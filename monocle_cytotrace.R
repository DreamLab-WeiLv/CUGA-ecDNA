library(Seurat)
library(dplyr)
library(monocle)
library(CytoTRACE)

seurat_data<- readRDS('prx.rds')
Idents(seurat_data) <- "seurat_clusters"

expr_matrix <- as(as.matrix(seurat_data@assays$RNA@counts),'sparseMatrix')
p_data <- seurat_data@meta.data
f_data <- data.frame(data.frame(gene_short_name = row.names(seurat_data),
                                row.names = row.names(seurat_data)))
#---- constract CDS object ----
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)

monocle_cds <- newCellDataSet(expr_matrix,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

#---- HVGs guide tradj ----
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

diff <- FindAllMarkers(seurat_data,logfc.threshold = 0.25,only.pos = T,min.pct = 0.25)
diff <- subset(diff,!grepl(pattern = 'RP[LS]',gene))
diff <- subset(diff,!grepl(pattern = 'MT-',gene))
diff <- diff %>% filter(p_val_adj<0.05)

ordergene <- unique(diff$gene)
monocle_cds <- setOrderingFilter(monocle_cds, ordergene)

monocle_cds <- reduceDimension(monocle_cds, max_components = 2,method = "DDRTree")
monocle_cds <- orderCells(monocle_cds)

#---- cytoTrace ----
results <- CytoTRACE(mat = expr, enableFast = T)

# ---- add cytoTrace to monocle2 ----
monocle_cds@phenoData@data$CytoTrace <- results$CytoTRACE

saveRDS(monocle_cds,'prx_cds.rds')
data <- pData(monocle_cds)
write.table(data,'allmetadata.tsv',sep= '\t',quote=F,row.names=T)
