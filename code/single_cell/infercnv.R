library(infercnv)
library(Seurat)

#Epi-separate by sample
epi <- readRDS('prx_Epithelial.rds')
seurat_Fibro <- readRDS('ref_Fibro.rds')
seurat_Endo <- readRDS('ref_Endo.rds')

seurat_object <- merge(x=epi,y=c(seurat_Fibro,seurat_Endo))
Idents(seurat_object) <- seurat_object$cell_type
counts <- as.data.frame(GetAssayData(seurat_object,assay ="RNA", slot = 'counts'))
anno <- data.frame(Idents(seurat_object))
gene_order <- "/home/projects/ku_00009/people/congli/EC/human_genes_pos.txt"
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts,
                                    annotations_file = anno,
                                    delim="\t",
                                    gene_order_file = gene_order,
                                    min_max_counts_per_cell = c(100, +Inf),
                                    ref_group_names = c("Fibroblast", "Endothelial"))
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,
                             out_dir = outdir,
                             cluster_by_groups = F,
                             k_obs_groups = 8,
                             HMM = T,
                             denoise = T,
                             write_expr_matrix=TRUE,
                             num_threads = 10)

##############################################################
##CNV score
setwd('/inferCNV')
dir_list <- list.dirs(recursive = F)

score <- data.frame()
for (i in dir_list){
  setwd(i)
  expr <- read.table("infercnv.observations.txt", header=T) %>% as.matrix()
  expr.scale <- scale(t(expr))
  tmp1 <- sweep(expr.scale, 2, apply(expr.scale, 2, min),'-')
  tmp2 <- apply(expr.scale, 2, max) - apply(expr.scale,2,min)
  expr_1 <- t(2*sweep(tmp1, 2, tmp2, "/")-1)
  cnv_score <- as.data.frame(colSums(expr_1 * expr_1))
  colnames(cnv_score)="cnv_score"
  cnv_score <- rownames_to_column(cnv_score, var='cell')
  score <- rbind(score,cnv_score)
  setwd("..")
}

write.table(score,'cnv_score.tsv',sep='\t',quote=F,row.names=F)
