library(dndscv)
library(data.table)
library(tidyverse)
df <- fread("prx.maf.hg19.gz",sep = "\t",header = T)
df_1 <- df[,c(13,5,6,11,12)]
colnames(df_1) = c("sampleID","chr","pos","ref","mut")
df_1$chr = str_replace(df_1$chr,"chr","")
Combine_448 = dndscv(df_1)
sel_cv = Combine_448$sel_cv
genes <- sel_cv[,c("gene_name","qglobal_cv")]
signif_genes = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
