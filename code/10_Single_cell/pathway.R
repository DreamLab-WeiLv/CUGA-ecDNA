library(Seurat)

seurat_data<- readRDS('Epithelial.rds')
##findmarker
Idents(seurat_data) <- 'ecDNA'
findmarkers <- Findfindmarkers(seurat_data,ident.1='Positive',ident.2='Negative',group.by="ecDNA",min.pct = 0,logfc.threshold = 0)
up1 <- rownames(findfindmarkers[which(findfindmarkers$avg_log2FC > 0 & findfindmarkers$p_val < 0.05),])
down1 <- rownames(findfindmarkers[which(findfindmarkers$avg_log2FC < 0 & findfindmarkers$p_val < 0.05),])

##pseudobulk
library(Libra)
pseudobulk <- run_de(seurat_data,de_family='pseudobulk',replicate = 'sample_name',cell_type_col = 'orig_cell_type', label_col = 'ecDNA')  #edgeR
pseudobulk <- na.omit(pseudobulk)
up2 <- pseudobulk[which(pseudobulk$avg_logFC < 0 & pseudobulk$p_val < 0.05),'gene']
down2 <- pseudobulk[which(pseudobulk$avg_logFC > 0 & pseudobulk$p_val < 0.05),'gene']

up <- intersect(up1,up2) 
down <- intersect(down1,down2) 

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tidyverse)

hsets <- read.gmt("hallmark_cancersea.gmt")

tmp <- enricher(down, TERM2GENE = hsets)
down.result <- tmp@result 
down.result$group <- 'down'
tmp <- enricher(up, TERM2GENE = hsets)
up.result <- tmp@result
up.result$group <- 'up'
enrich.result <- rbind(up.result,down.result)
enrich.result <- enrich.result[enrich.result$pvalue <= 0.05,]
enrich.result <- enrich.result %>% 
  group_by(group) %>% 
  top_n(n = 5,wt = -qvalue)

enrich.result$group <- factor(enrich.result$group, levels = c("up","down"))
enrich.result <- enrich.result[order(enrich.result$group), ]
enrich.result$Description <- factor(enrich.result$Description, levels = enrich.result$Description)
enrich.result$geneID  <- sapply(strsplit(enrich.result$geneID , "/"), function(x) paste(x[1:5], collapse = "/"))

enrich.result.up <- enrich.result[enrich.result$group=='up',]
enrich.result.down <- enrich.result[enrich.result$group=='down',]
p1<-ggplot(enrich.result.up, aes(x = -log10(qvalue), y = rev(Description),fill=group))+
  geom_bar(stat = "identity", width = 0.5)+
  geom_text(aes(x = 0.1,y=rev(Description),label = Description),size=3.5, hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(colour = 'black', size = 12),
        axis.line = element_line(colour = 'black', linewidth =0.5),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 12),
        legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_manual(values = '#AA3538')+ #'#5891BF'
  geom_text(data = enrich.result.up,
            aes(x = 0.1, y = rev(Description), label = geneID,color=group),
            size = 4,
            fontface = 'italic', 
            hjust = 0,
            vjust = 2.3)+
  scale_color_manual(values = c('#AA3538'))+ #'#5891BF'
  scale_y_discrete(expand = c(0.12,0))+
  labs(title = "Enrichment of genes",
       y=c(" Up"))
p1
p2<-ggplot(enrich.result.down, aes(x = -log10(qvalue), y = rev(Description),fill=group))+
  geom_bar(stat = "identity", width = 0.5)+
  geom_text(aes(x = 0,y=rev(Description),label = Description),size=3.5, hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(colour = 'black', size = 12),
        axis.line = element_line(colour = 'black', linewidth =0.5),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 12),
        legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_manual(values = '#5891BF')+ #
  geom_text(data = enrich.result.down,
            aes(x = 0.17, y = rev(Description), label = geneID,color=group),
            size = 4,
            fontface = 'italic', 
            hjust = 0.1,
            vjust = 2.3)+
  scale_color_manual(values = c('#5891BF'))+ #
  scale_y_discrete(expand = c(0.12,0))+
  labs(title = "Enrichment of genes",
       y=c("Down"))
p2

p<- p1+p2
ggsave('pathway_hallmark.pdf',p,w=10,h=4.5)
