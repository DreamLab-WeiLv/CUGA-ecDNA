# k_selection #
setwd('/share/home/luoylLab/zengyuchen/PROJECT/chromothripsis/36.ec_scRNA/06.UC1/064.cnmf')
library(tidyverse)
dir_input <- "./res1" #
dir_output <- "./res2" #
dir.create(dir_output)
dirs <- setdiff(dir(dir_input),c("k_selection","k_selection.txt")) #
ref.file <- read.table(paste0(dir_input,"/k_selection.txt"), header = F,sep = " ",stringsAsFactors = F)
colnames(ref.file)=c("sample","k")
rownames(ref.file)=ref.file$sample

for (i in dirs) {
  system(paste("python -W ignore cnmf.py consensus --output-dir ",dir_input," --name ",i," --components ",ref.file[i,"k"]," --local-density-threshold 0.02 --show-clustering",sep = ""))
  
  path1 <- paste0(dir_input,"/",i,"/",i,".usages.k_",ref.file[i,"k"],".dt_0_02.consensus.txt")
  path2 <- paste0(dir_input,"/",i,"/",i,".gene_spectra_score.k_",ref.file[i,"k"],".dt_0_02.txt")
  
  if( file.exists(path1) & file.exists(path2) ){
    system(paste("cp ",path1," ",path2," ",dir_output,"/",sep = ""))
  }
}
top_gene <- 50 
for (i in dirs) {
  usage.file <- dir(dir_output,pattern = paste0(i,".usages"))
  usage.df <- read.table(paste(dir_output,"/",usage.file,sep = ""),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
  colnames(usage.df) <- paste(i,1:dim(usage.df)[2],sep = "_")
  
  #normalize
  usage.df <- usage.df / rowSums(usage.df) 
  write.table(usage.df,file = paste(dir_output,"/",i,"_program.usage.norm.txt",sep = ""),quote = F,sep = "\t",row.names = T,col.names = T)
  
  #QC1
  tmpdf1 <- gather(usage.df,"program","ratio")
  tmpdf1 %>% ggplot(aes(x=program,y=ratio))+geom_boxplot(outlier.shape = NA)+geom_jitter(color="red",alpha=0.4,width = 0.2)+
    labs(title = i)+
    theme(
      axis.text.x.bottom = element_text(angle = 45,hjust = 1),
      plot.title = element_text(hjust = 0.5,size=20)
    )
  ggsave(paste(dir_output,"/",i,"_program.usage.norm.QC.png",sep = ""),device = "png",width = 20,height = 16,units = c("cm"))
  
  #score
  score.file <- dir(dir_output,pattern = paste0(i,".gene_spectra_score"))
  score.df <- read.table(paste(dir_output,"/",score.file,sep = ""),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
  score.df <- as.data.frame(t(score.df))
  colnames(score.df) <- paste(i,1:dim(score.df)[2],sep = "_")
  
  topn.df <- as.data.frame(matrix(nrow = top_gene,ncol = ncol(score.df)))
  colnames(topn.df) <- colnames(score.df)
  
  for (k in colnames(score.df)) {
    tmpv <- score.df[,k]
    names(tmpv) <- rownames(score.df)
    topn.df[,k] <- names(rev(tail(sort(tmpv),top_gene)))
  }
  
  #save
  write.table(topn.df, file = paste(dir_output,"/",i,"_program.Zscore.top",top_gene,"gene.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = T)
  score.df$gene=rownames(score.df)
  write.table(score.df,file = paste(dir_output,"/",i,"_program.Zscore.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = T)
}

# check usages program (mean ratio > 0.03) ###
check.usage <- data.frame()
usage_filter <- 0.03
for (i in dirs) {
  usage.file <- paste(dir_output,"/",i,"_program.usage.norm.txt",sep = "")
  usage.df <- read.table(usage.file,header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
  check.usage <- rbind(check.usage,as.data.frame(colMeans(usage.df)))
}
colnames(check.usage) <- c("mean_ratio")
check.usage$sample_programs <- rownames(check.usage)
check.usage <- check.usage %>% arrange(mean_ratio)
check.usage$sample_programs <- factor(check.usage$sample_programs,levels = check.usage$sample_programs)
linex <- sum(check.usage$mean_ratio < usage_filter)
check.usage %>% ggplot(aes(x=sample_programs,y=mean_ratio))+geom_point()+
  geom_hline(yintercept = usage_filter,color="red")+
  geom_vline(xintercept = linex+0.5,color="red")+
  theme(
    axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)
  )
ggsave(paste(dir_output,"/","check.usage.png",sep = ""),width = 30,height = 16,device = "png",units = "cm")
maybe.bg <- as.character(check.usage$sample_programs[check.usage$mean_ratio < usage_filter])

library(pheatmap)
library(RColorBrewer)
library(scales)

all.score.df <- data.frame()
all.score.topn.df <- data.frame()
for (i in dirs) {
  score.file <- paste(dir_output,"/",i,"_program.Zscore.txt",sep = "")
  score.df <- read.table(score.file,header = T,sep = "\t",stringsAsFactors = F)
  if (i==dirs[1]) {all.score.df=score.df}
  if (i!=dirs[1]) {
    all.score.df <- all.score.df %>% inner_join(score.df,by="gene")
  } # merge all score files 
  
  score.topn.file <- paste(dir_output,"/",i,"_program.Zscore.top",top_gene,"gene.txt",sep = "")
  score.topn.df <- read.table(score.topn.file,header = T,sep = "\t",stringsAsFactors = F)
  if (i==dirs[1]) {all.score.topn.df=score.topn.df}
  if (i!=dirs[1]) {
    all.score.topn.df <- cbind(all.score.topn.df,score.topn.df)
  }
}

rownames(all.score.df) <- all.score.df$gene
all.score.df$gene <- NULL
all.score.df <- all.score.df[rowSums(is.na(all.score.df)) == 0,] #可能有空值，需要去掉
all.score.rm.df <- all.score.df[,setdiff(colnames(all.score.df),maybe.bg)] #在质控这一步检测出来的噪声
all.score.rm.df.cor <- cor(all.score.rm.df,method = "pearson")

cor_min <- 0
cor_max <- 0.6
all.score.rm.df.cor[all.score.rm.df.cor < cor_min] <- cor_min
all.score.rm.df.cor[all.score.rm.df.cor > cor_max] <- cor_max

colanno <- as.data.frame(colnames(all.score.rm.df.cor))
colnames(colanno) <- "colnames"
colanno$sample <- str_replace(colanno$colnames,"_.*","")
rownames(colanno) <- colanno$colnames
colanno$colnames <- NULL

rowanno <- as.data.frame(rownames(all.score.rm.df.cor))
colnames(rowanno) <- "rownames"
rowanno$sample <- str_replace(rowanno$rownames,"_.*","")
rownames(rowanno) <- rowanno$rownames
rowanno$rownames <- NULL

cluster_method <- "complete"
scale_min <- -2
scale_max <- 2
cluster_method2 <- "complete"
color <- NULL
#color
if (is.null(color)){
  color_v <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(dirs))
  names(color_v) <- dirs
}else{
  color_v <- color
}
ann_colors  <- list(sample = color_v)
#visual
tmpp <- pheatmap(all.score.rm.df.cor,cluster_rows = T,cluster_cols = T,
              clustering_method = cluster_method, 
              show_colnames = F,
              treeheight_row=30,treeheight_col=0,
              border_color=NA,
              annotation_row = rowanno,annotation_col = colanno,
              annotation_names_row = F,annotation_names_col = F,
              #annotation_colors = ann_colors,
              color = colorRampPalette(c("white","yellow", "red","#67001F"))(50),
              fontsize_row=12,
              width = 11.5,height = 9,
              filename = paste0(dir_output,"/","program_pearson_cor.",cluster_method,".heatmap.pdf"))

#save visual data
write.table(all.score.rm.df.cor,file = paste0(dir_output,"/","cor_heatmap_data.txt"),quote = F,sep = "\t",row.names = T,col.names = T)
all.score.topn.rm.df <- all.score.topn.df[,setdiff(colnames(all.score.topn.df),maybe.bg)]#
write.table(all.score.topn.rm.df,file = paste0(dir_output,"/","program_topngene.txt"),quote = F,sep = "\t",row.names = F,col.names = T)

# annotation
library(clusterProfiler)
library(org.Hs.eg.db)
hsets <- read.gmt("hallmark_cancersea.gmt")
enrich.result <- data.frame()
pathway_v <- c()
program_v <- c()
program_topn <- read.table(paste(dir_output,"/","program_topngene.txt",sep = ""),header = T,sep = "\t",stringsAsFactors = F)
for (i in 1:dim(program_topn)[2]) {
  tmp <- enricher(program_topn[,i], TERM2GENE = hsets)
  if (is.null(tmp)) {
    next
  }
  tmp1 <- head(tmp@result)
  tmp1$program <- colnames(program_topn)[i]
  rownames(tmp1) <- NULL
  enrich.result <- rbind(enrich.result,tmp1)
  
  program_v <- append(program_v,colnames(program_topn)[i])
  pathway_v <- append(pathway_v,paste(tmp1$Description,collapse = ","))
}

write.table(enrich.result,file = paste0(dir_output,"/","program_topngene_enrichment.tsv"),row.names = F)

enrich.df <- data.frame(program=program_v,pathway=pathway_v)
enrich.df$program <- factor(enrich.df$program,levels = tmpp$tree_row$labels[tmpp$tree_row$order])
enrich.df <- enrich.df%>%arrange(program)
write.csv(enrich.df,file = paste0(dir_output,"/","program_topngene_enrichment_order.csv"),row.names = F,quote = F)

dir_count <- getwd()
for (i in dirs) {
  
  one_matrix <- read.table(paste(dir_count,"/",i,".count.txt",sep = ""),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
  RowSum <- rowSums(one_matrix)
  one_matrix <- log1p((one_matrix / RowSum) * 10000)
  one_matrix <- as.data.frame(t(one_matrix))
  
  plot_gene <- read.table(paste(dir_output,"/","program_topngene.txt",sep = ""),header = T,sep = "\t",stringsAsFactors = F)
  plot_gene <- plot_gene[,str_detect(colnames(plot_gene),i)]
  plot_gene <- gather(plot_gene,key = "program",value = "gene")
  plot_gene <- plot_gene%>%arrange(program)
  
  one_matrix <- one_matrix[plot_gene$gene,]
  one_matrix <- t(scale(t(one_matrix)))
  one_matrix[one_matrix < scale_min]  <- scale_min
  one_matrix[one_matrix > scale_max] <- scale_max
  
  if(length(unique(plot_gene$gene)) < length(plot_gene$gene)){
    plot_gene$gene <- rownames(one_matrix)
  }
  pn <- length(unique(plot_gene$program))
  tmpp2 <- pheatmap(one_matrix,cluster_rows = F,cluster_cols = T,
                 color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
                 treeheight_col=0,
                 clustering_method = cluster_method2,
                 show_colnames=F,
                 border_color=NA,
                 gaps_row=as.numeric(cumsum(table(plot_gene$program))[- pn])
  )
  write.table(one_matrix,paste(dir_output,"/",i,"_data_heatmap.txt",sep = ""),quote = F,sep = "\t",row.names = T,col.names = T)
  
  cell_sort <- tmpp2$tree_col$labels[tmpp2$tree_col$order]
  gene_sort <- plot_gene$gene
  matrix_new <- as.data.frame(one_matrix)
  matrix_new$gene <- rownames(matrix_new)
  matrix_new <- matrix_new%>%reshape2::melt(id="gene")
  colnames(matrix_new)[c(2,3)] <- c("cell","exp")
  matrix_new$gene <- factor(matrix_new$gene,levels = rev(gene_sort))
  matrix_new$cell <- factor(matrix_new$cell,levels = cell_sort)
  
  plot1 <- matrix_new%>%ggplot(aes(x=cell,y=gene,fill=exp))+geom_tile()+
    geom_hline(yintercept = as.numeric(cumsum(table(plot_gene$program))[- pn])+0.5,color="black",linetype=5)+
    labs(title = paste(i,": ",length(cell_sort)," cells; ",pn," programs",sep = ""))+
    scale_x_discrete("",expand = c(0,0))+
    scale_y_discrete("",expand = c(0,0))+
    scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x.bottom = element_blank(),axis.text.y.left = element_blank(),
      plot.title = element_text(hjust = 0.5,size = 20),
      legend.position = "left"
    )
  
  dim1 <- 5*pn
  dim2 <- (top_gene / 10) *2
  gene.text <- t(matrix(plot_gene$gene,nrow = dim2,ncol = dim1))
  rownames(gene.text) <- seq(-0.5,-(dim1-0.5),-1)
  colnames(gene.text) <- seq(0.5,(dim2-0.5),1)
  gene.text <- as.data.frame(gene.text)
  gene.text$dim1 <- rownames(gene.text)
  gene.text <- reshape2::melt(gene.text,id="dim1")
  colnames(gene.text)[2:3] <- c("dim2","gene")
  gene.text$dim1 <- as.numeric(as.character(gene.text$dim1))
  gene.text$dim2 <- as.numeric(as.character(gene.text$dim2))
  
  plot2 <- gene.text%>%ggplot(aes(x=dim2,y=dim1))+geom_text(aes(label=gene))+
    geom_hline(yintercept = seq(-dim1,0,5)[-c(1,length(seq(-dim1,0,5)))],color="black",linetype=5)+
    scale_x_continuous("",expand = c(0,0),limits = c(0,10))+
    scale_y_continuous("",expand = c(0,0),limits = c(-dim1,0))+
    labs(title = paste(unique(plot_gene$program),collapse = "; "))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      plot.title = element_text(hjust = 0.5,size = 20)
    )
  
  library(patchwork)
  plot3 <- plot1+plot2+plot_layout(widths = c(1,2))
  ggsave(filename = paste(dir_output,"/",i,"_program_gene.heatmap.pdf",sep = ""),plot = plot3,width = 46,height = 16,units = "cm")
}

