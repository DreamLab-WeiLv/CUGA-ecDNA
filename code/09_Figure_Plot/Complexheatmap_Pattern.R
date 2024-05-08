library(readxl)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(data.table)
library(maftools)
library(ComplexHeatmap)
library(circlize)
## 1. Input maf files ##
ddwes_maf <- as.data.frame(fread('maf.tsv',sep = '\t',header = T,check.names = F))
info <- as.data.frame(read_xlsx('cli.xlsx'))
ddwes_maf <- ddwes_maf %>%
  filter(ddwes_maf$Tumor_Sample_Barcode%in%info$SampleName)
table(ddwes_maf$Tumor_Sample_Barcode)
ddwes_maf$CaseID <- NA
matches <- match(ddwes_maf$Tumor_Sample_Barcode, info$SampleName)
ddwes_maf$CaseID <- info$CaseID[matches]
maf_MR <- fread('maf.tsv',sep = '\t',data.table = F,header = T)
length(table(maf_MR$Tumor_Sample_Barcode)) 
single_sample <- c('sample1','sample2','sample3')
maf_MR <- maf_MR %>%
  filter(!maf_MR$Tumor_Sample_Barcode%in%single_sample)
length(table(maf_MR$Tumor_Sample_Barcode)) 
maf_merge <- rbind(maf_MR,ddwes_maf) 
table(maf_merge$Tumor_Sample_Barcode)
length(table(maf_merge$Tumor_Sample_Barcode))
maf_1 <- read.maf(maf_merge)
maf_2 <- as.data.frame(maf_1@data[,c(1,9,13)])
mut <- reshape(maf_2,idvar = "Hugo_Symbol",timevar = "Tumor_Sample_Barcode",direction = "wide")
rownames(mut) <- mut[,1]
mut <- mut[,-1]
colnames(mut) <- str_replace(colnames(mut),"Variant_Classification.","")
mut_1 <- as.data.frame(lapply(mut,as.character))
mut_1[is.na(mut_1)] = ""
rownames(mut_1) <- rownames(mut) 
genelist=read_xlsx('genelist.xlsx',sheet = 4)
snv_list = genelist[1:10,1]
snv <- mut_1[rownames(mut_1) %in% snv_list$Gene,]
sample_order <- info$SampleName 
snv[,198:211] <- ""
colnames(snv)[198:211] <- info$SampleName[198:211]
snv <- snv[,sample_order]
info <- as.data.frame(read_xlsx('cli.xlsx',sheet = 8))
info$CaseID <- factor(info$CaseID,levels = c("CUGA-006","CUGA-009","CUGA-044",
                                             "CUGA-MR-001","CUGA-MR-002","CUGA-MR-003","CUGA-MR-004","CUGA-MR-005","CUGA-MR-006","CUGA-MR-007","CUGA-MR-008",
                                             "CUGA-MR-009","CUGA-MR-010","CUGA-MR-012","CUGA-MR-013","CUGA-MR-014","CUGA-MR-015","CUGA-MR-016","CUGA-MR-017",
                                             "CUGA-MR-018","CUGA-MR-019","CUGA-MR-020","CUGA-MR-022","CUGA-MR-023","CUGA-MR-024","CUGA-MR-026","CUGA-MR-027",
                                             "CUGA-MR-028","CUGA-MR-029","CUGA-MR-030","CUGA-MR-031","CUGA-MR-032","CUGA-MR-033","CUGA-MR-034","CUGA-MR-035",
                                             "CUGA-MR-036","CUGA-MR-037","CUGA-MR-038","CUGA-MR-040","CUGA-MR-042","CUGA-MR-044","CUGA-MR-045","CUGA-MR-046",
                                             "CUGA-MR-047","CUGA-MR-048","CUGA-MR-049","CUGA-MR-050","CUGA-MR-051","CUGA-MR-052","CUGA-MR-053","CUGA-MR-054",
                                             "CUGA-MR-055","CUGA-MR-056","CUGA-MR-058","CUGA-MR-039","CUGA-MR-043","CUGA-MR-057"))
are_equal <- identical(colnames(snv), info$SampleName)

if (are_equal) {
  cat("same")
} else {
  cat("no same")
}

col_fun_1 = c("Frame_Shift_Del" = "#7f83bf",
              "Frame_Shift_Ins" = "#7f83bf",
              "In_Frame_Del" = "#7f83bf",
              "In_Frame_Ins" = "#7f83bf",
              "Missense_Mutation" = "#7f83bf",
              "Nonsense_Mutation" = "#7f83bf",
              "Splice_Site" = "#7f83bf",
              "Nonstop_Mutation" = "#7f83bf")
alter_fun_1 = list(
  background = alter_graphic('rect',horiz_margin = unit(0,'pt'),
                             vertical_margin = unit(0,'pt'),col= '#f0f0f0',fill = 'white',lwd=0),
  Frame_Shift_Del = alter_graphic("rect",horiz_margin = unit(0,"pt"),
                                  vertical_margin = unit(0,"pt"),
                                  fill = col_fun_1["Frame_Shift_Del"]),
  Frame_Shift_Ins = alter_graphic("rect",horiz_margin = unit(0,"pt"),
                                  vertical_margin = unit(0,"pt"),
                                  fill = col_fun_1["Frame_Shift_Ins"]),
  In_Frame_Del = alter_graphic("rect",horiz_margin = unit(0,"pt"),
                               vertical_margin = unit(0,"pt"),
                               fill = col_fun_1["In_Frame_Del"]),
  In_Frame_Ins = alter_graphic("rect",horiz_margin = unit(0,"pt"),
                               vertical_margin = unit(0,"pt"),
                               fill = col_fun_1["In_Frame_Ins"]),
  Missense_Mutation = alter_graphic("rect",horiz_margin = unit(0,"pt"),
                                    vertical_margin = unit(0,"pt"),
                                    fill = col_fun_1["Missense_Mutation"]),
  Nonsense_Mutation = alter_graphic("rect",horiz_margin = unit(0,"pt"),
                                    vertical_margin = unit(0,"pt"),
                                    fill = col_fun_1["Nonsense_Mutation"]),
  Splice_Site = alter_graphic("rect",horiz_margin = unit(0,"pt"),
                              vertical_margin = unit(0,"pt"),
                              fill = col_fun_1["Splice_Site"]),
  Nonstop_Mutation = alter_graphic("rect",horiz_margin = unit(0,"pt"),
                                   vertical_margin = unit(0,"pt"),
                                   fill = col_fun_1["Nonstop_Mutation"]))

location_color <- c("Bladder" = "#c5996d","LN" = "#ce2724",'Pelvis'='#f0b8b8','Ureter'='#b89bc0','Urethra'='#bd6b85')
ecDNA_color <- c("Positive" = "#263366","Negative"="white")
chromo_color <- c("Yes" = "#76cbea", "No" = "white","NA" = "lightgrey")
BFB_color <- c("Positive" = "#e98e29","Negative"="white")
Linear_color <- c("Positive" = "#aeb9c5","Negative" = "white")
Complex_color <- c("Positive" = "#8283a7","Negative" = "white")
WGD_color <- c("TRUE" = "#d1a4ac", 'FALSE' = "white", 'NA' = 'lightgrey')

clinical_1 <- HeatmapAnnotation(#CaseID = info$CaseID,
  Location = info$Location,
  ecDNA = info$ecDNA,
  BFB = info$BFB,
  Linear = info$Linear,
  Complex = info$Complex,
  Chromothripsis = info$Chromothripsis,
  WGD = info$WGD,
  col = list(ecDNA = ecDNA_color,
             Location = location_color,
             BFB = BFB_color,
             Linear = Linear_color,
             Complex = Complex_color,
             Chromothripsis = chromo_color,
             WGD = WGD_color),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize =10),
  simple_anno_size = unit(0.5,"cm"),
  show_legend = T,
  border = T,
  gp = gpar(col='#f6f6f6',lwd=0.1),
  gap = unit(1.2, "mm"),
  annotation_legend_param = list(
    # CaseID = list(nrow=3),ecDNA = list(nrow=1),
    Location = list(nrow=1),BFB = list(nrow=1),
    Linear = list(nrow=1),Complex = list(nrow=1),
    Chromothripsis = list(nrow=1),WGD = list(nrow=1)))
rownames(snv) <- factor(rownames(snv),levels = c('TP53','FGFR3','CREBBP','EP300','ERBB3',
                                                 'ERCC2','HRAS','KDM6A','KMT2D','PIK3CA'))
pdf('Figure.pdf',width = 20,height = 8)
pf <- oncoPrint(snv,
                border = T,
                row_split = factor(rownames(snv),levels = c('TP53','FGFR3','CREBBP','EP300','ERBB3',
                                                            'ERCC2','HRAS','KDM6A','KMT2D','PIK3CA')),
                column_order = sample_order,
                alter_fun = alter_fun_1,
                col = col_fun_1,
                na_col = 'lightgrey',
                heatmap_legend_param = list(title = 'Alterations',nrow=1),
                alter_fun_is_vectorized = F,
                show_pct = F,
                height = unit(5,"cm"),
                width = unit(45,'cm'),
                row_names_side = "left",
                top_annotation = clinical_1,
                column_split = factor(info$CaseID),
                row_names_gp = gpar(fontface="italic",fontsize = 10),
                right_annotation = NULL)
draw(pf,heatmap_legend_side="bottom",annotation_legend_side="bottom")
dev.off()
