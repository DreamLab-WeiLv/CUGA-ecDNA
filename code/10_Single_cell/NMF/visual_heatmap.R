cor <- read.table('./res2/cor_heatmap_data.txt',sep='\t',header = T,row.names = 1)

annotation_col <- data.frame(ID=colnames(cor))
annotation_col$ID <- factor(annotation_col$ID,levels=unique(annotation_col$ID))

mycol1 <- c('#4bb3d2','#fac6e0','#87cab9','#344981',"#f78171",'#c3a8d3','#e54034',"#f7ac5c",'#01967a','#c6e6c0','#a7d95c','#5f368b','#b376b1')
pati_col <-setNames(mycol1,c('P1','P10','P11','P13','P14','P15','P16','P17','P19','P3','P6','P7','P9'))

p <- pheatmap(as.matrix(cor),cluster_rows = T,cluster_cols = T,
                 clustering_method = 'complete', 
                 show_colnames = F,
                 treeheight_row=30,treeheight_col=0,
                 border_color = NA,
                 border=TRUE,
                 annotation_col = annotation_col,
                 annotation_colors = pati_col,
                 annotation_names_row = F,annotation_names_col = F,
                 color = colorRampPalette(c("white","yellow", "red","#67001F"))(50),
                 fontsize_row=12)
ggsave('program_pearson_cor.heatmap.pdf',p,width = 11, height = 9)

##########################################
##module signature gene
library(tidyverse)
mat.files=dir("./res2/",pattern='dt_0_02.txt$')
all.mat <- data.frame()
for (fi in mat.files){
  tmp.mat <- read.table(paste0('./res2/',fi),header = T,sep='\t',row.names = 1,stringsAsFactors = F)
  tmp.mat <- as.data.frame(t(tmp.mat))
  sampleid <- str_replace(fi,"_.*$","")
  colnames(tmp.mat) <- paste(sampleid,colnames(tmp.mat),sep='_')
  tmp.mat$gene <- rownames(tmp.mat)
  
  if (sampleid == '001TB'){
    all.mat=tmp.mat
  }else{
    all.mat <- all.mat %>% full_join(tmp.mat,by="gene")
  }
}

write.table(all.mat,'all_program_signature.tsv',sep='\t',quote=F,row.names = F)

#program gene
MP<- as.list(readxl::read_xlsx('MP_programs.xlsx'))
MP_gene <- list()
for (i in 1:5){
  signature.programs <- na.omit(MP[[i]])
  signature.loading <- all.mat.tmp[,c("gene",signature.programs)]
  
  used.gene=c()
  for (pi in signature.programs){
    tmp.df = signature.loading[,c('gene',pi)]
    tmp.loading = tmp.df[,2]
    names(tmp.loading)=tmp.df[,1]
    
    tmp.loading=tmp.loading[!is.na(tmp.loading)]
    used.gene=append(used.gene,names(tail(sort(tmp.loading),100)))
  }
  used.gene=unique(used.gene) 
  
  signature.loading=signature.loading[signature.loading$gene %in% used.gene,]
  rownames(signature.loading)=signature.loading$gene
  signature.loading$gene=NULL
  signature.loading$total_loading=rowSums(signature.loading)
  signature.loading$average_loading=signature.loading$total_loading/length(signature.programs)
  signature.loading=signature.loading %>% arrange(desc(average_loading))
  #mp_top <- list(rownames(signature.loading)[1:30])
  mp_top <- list(rownames(signature.loading))
  names(mp_top) <- paste0('MP',i)
  MP_gene <- c(MP_gene,mp_top)
}
MP_gene <- as.data.frame(MP_gene)

max_length <- max(sapply(MP_gene, length))

#
filled_list <- lapply(MP_gene, function(x) {
  if (length(x) < max_length) {
    c(x, rep(NA, max_length - length(x)))
  } else {
    x
  }
})

df <- as.data.frame(filled_list)
write.table(df,'MP_signature.tsv',sep='\t',quote=F,row.names = F) #MP_signature_top30.tsv

########################################################
MP_signature <- read.table('MP_signature.tsv',header=T,sep='\t')
gene <- unique(MP_signature[,2]) #
path='/home/projects/ku_00009/people/congli/EC/noAT/subtype/filter1/epi_harmony/new_13/NMF/res2/'
files=dir(path,'program.Zscore.txt')
df <- data.frame(gene=gene)
for (i in files) {
  data <- read.table(paste0(path,i),sep='\t',header=T)
  data <- data[data$gene %in% gene,]
  df <- merge(df,data,by='gene')
}
write.table(df,'NMF_score_heatmap.csv',sep=',')

#############################################
#NMF score heatmap
NMF_score <- read.table('NMF_score_heatmap.csv',sep=',',header = T)
colanno<- read.table('program_sample_name.tsv',sep='\t',header = T)
colnames(NMF_score) <- colanno[match(colnames(NMF_score),colanno$colnames),'new_name']

or <- read.table('./res2/program_topngene_enrichment_order.csv',sep=',',header = T)
matrix <- subset(NMF_score,select=or$program_new)
matrix$gene <- NMF_score[,1]

pdf('NMF_score_heatmap.pdf')
Heatmap(as.matrix(matrix[,1:50]),
        col = colorRamp2(c(-0.005,0,0.005),c("#2268ad",'#F0F0F0',"#b92732")),
        cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
        row_split = factor(rep(paste0("MP",1:5),c(510,417,517,148,565)),levels = paste0("MP",1:5)), 
        column_gap = unit(0, "mm"),
        row_gap = unit(0, "mm"),
        column_split = factor(rep(paste0("MP",1:5),c(16,8,14,4,8)),levels = paste0("MP",1:5)),
        heatmap_legend_param = list(title = "NMF score",direction = "vertical",title_position = "leftcenter-rot",at=c(-0.01,0,0.01),legend_height = unit(3, "cm")),  
        use_raster = TRUE,
        border = TRUE,
        row_title = NULL,column_title = NULL)
dev.off()

###############################################################
##pathway
MP <- read.table('MP_signature.tsv',sep='\t',header = T)
colnames(MP) <- c('group','gene')
MP_gene <- as.list(split(MP$gene,MP$group))
hsets <- read.gmt("D:/work/pan cancer_EC/EC/new0229/noAT/epi/NMF/hallmark_cancersea.gmt") #hall +cancersea(14)
enrich.result <- data.frame()
for (i in 1:5){
  gene_list <- MP_gene[[i]]
  tmp <- enricher(gene_list, TERM2GENE = hsets)
  tmp1 <- head(tmp@result)
  tmp1$Description <- tolower(tmp1$Description)
  tmp1$Description <- capitalize(tmp1$Description)
  tmp1$program <- names(MP_gene)[i]
  enrich.result <- rbind(enrich.result,tmp1)
}
write.table(enrich.result,"MP_pathway_hallmark.csv",sep=',',row.names = F)

#msigdbr
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(msigdbr)

hs <- msigdbr(species = "Homo sapiens")
hs[1:5,1:5]
table(hs$gs_subcat)

hs_kegg = msigdbr(species = "Homo sapiens",
                  category = "C2", 
                  subcategory = "CP:KEGG") %>%
                  dplyr::select(gs_name,gene_symbol)
hs_go_bp = msigdbr(species = "Homo sapiens",
                category = "C5", 
                subcategory = "GO:BP" ) %>%
                dplyr::select(gs_name,gene_symbol)
hs_go_cc = msigdbr(species = "Homo sapiens",
                category = "C5",
                subcategory = "GO:CC" ) %>%
                dplyr::select(gs_name,gene_symbol)

library(Hmisc)
for (i in 1:5){
  gene_list <- MP_gene[[i]]
  em_kegg <- as.data.frame(enricher(gene_list, TERM2GENE=hs_kegg))
  em_go_bp <- as.data.frame(enricher(gene_list, TERM2GENE=hs_go_bp))
  em_go_cc <- as.data.frame(enricher(gene_list, TERM2GENE=hs_go_cc))
  result <- rbind(em_kegg,em_go_bp,em_go_cc)
  result$ID <- tolower(result$Description) %>% gsub('_',' ',.)  %>% capitalize()
  result_split <- separate(result, ID, into = c("subcategory", "name"), sep = " ", extra = "merge",remove = T) 
  result_split$name <-  capitalize(result_split$name)
  write.table(result_split,paste0('MP',i,'_enrichment.csv'),sep = ',',row.names = F)
}





