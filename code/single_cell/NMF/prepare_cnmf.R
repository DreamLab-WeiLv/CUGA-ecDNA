library(Seurat)
library(dplyr)
library(stringr)

epi <- readRDS('Epithelial.rds')
epi$CB <- rownames(epi@meta.data)
### export sample matrix 
for (si in as.character(unique(epi@meta.data$sample_name))) {
  small.meta.data <- epi@meta.data %>% filter(sample_name == si)
  small.count <- as.data.frame(epi@assays$RNA@counts[,small.meta.data$CB])
  small.count <- small.count[rowSums(small.count) > 0,]
  small.count <- small.count[!str_detect(rownames(small.count), "^MT-"),]
  small.count <- data.frame(t(small.count))
  write.table(small.count,file = paste0("home/NMF",si,"count.txt"),quote = F,sep = "\t",row.names = T,col.names = T)
}
