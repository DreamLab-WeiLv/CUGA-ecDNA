library(data.table)
library(dplyr)
library(readxl)
library(tidyverse)
# Input all patients info, caculate counts and tumor type
tcga_meta <- read.table("~/../../Volumes/ed/TCGA-PCAWG/TCGA-filtered.tar/results/other_files/manifest/manifest-TCGA-aaSuite_somatic_ss.csv",sep = ",",header = T)
pcawg_meta <- read.table("~/../../Volumes/ed/TCGA-PCAWG/PCAWG-filtered/other_files/manifest/manifest-PCAWG-aaSuite_somatic_ss_specimen_id_added_07192023.csv",sep = ",",header = T)
# Filter PCAWG and TCGA duplicate patients 
info <- fread('~/../../Volumes/ed/TCGA-filter-results/icgc-dataset-1702364647893/sample.tsv.gz',sep = '\t',data.table = F,header = T)
tcga_meta <- tcga_meta[,c(9:11,15,16)]
colnames(tcga_meta)[4:5] <- c('Tumor_barcode','Normal_barcode')
pcawg_meta <- pcawg_meta[,c(9:11,13,14)]
colnames(pcawg_meta)[4:5] <- c('Tumor_barcode','Normal_barcode')
length(table(tcga_meta$patient_barcode))
# add a column with TCGA ID in PCAWG files
matched <- match(pcawg_meta$patient_barcode,info$icgc_donor_id)
pcawg_meta$TCGAid <- info$submitted_donor_id[matched]
# remove dup-samples with TCGA and PCAWG 
need_rm <- pcawg_meta[pcawg_meta$TCGAid%in%tcga_meta$patient_barcode,]
pcawg_meta <- pcawg_meta[!pcawg_meta$patient_barcode%in%need_rm$patient_barcode,]
length(table(pcawg_meta$patient_barcode)) # 2033 patients/ 2156 samples
all_meta <- rbind(tcga_meta,pcawg_meta[,c(1:5)]) 
length(unique(all_meta$patient_barcode)) # 3973个patients
num_pcawg <- as.data.frame(table(pcawg_meta[!duplicated(pcawg_meta$patient_barcode),]$project))
num_tcga <- as.data.frame(table(tcga_meta[!duplicated(tcga_meta$patient_barcode),]$project))
all_type <- as.data.frame(table(all_meta$project[!duplicated(all_meta$patient_barcode)]))

# Input AA results 
aa_pcawg <- read.table('~/../../Volumes/ed/TCGA-PCAWG/PCAWG-filtered/aggregated_results.csv',sep = ',',header = T,row.names = 1,check.names = F)
aa_tcga <- read.table('~/../../Volumes/ed/TCGA-PCAWG/TCGA-filtered.tar/results/aggregated_results.csv',sep = ',',header = T,row.names = 1,check.names = F)
# PCAWG AA output
aa_pcawg <- aa_pcawg[!is.na(aa_pcawg$`AA amplicon number`),]
tmp1 <- strsplit(aa_pcawg$`Sample name`,"-")
aa_pcawg$patient_barcode <- sapply(tmp1,function(x) paste0(x[1]))
aa_pcawg <- aa_pcawg[,c(25,2,4)]
matched <- match(aa_pcawg$patient_barcode,info$icgc_donor_id)
aa_pcawg$New_patient <- info$submitted_donor_id[matched]
aa_pcawg$patient_barcode[!is.na(aa_pcawg$New_patient)] <- aa_pcawg[!is.na(aa_pcawg$New_patient),]$New_patient
aa_pcawg <- aa_pcawg[,-4]
# TCGA AA output
aa_tcga <- aa_tcga[!is.na(aa_tcga$`AA amplicon number`),]
aa_tcga$patient_barcode <- sapply(strsplit(aa_tcga$`Sample name`, "-"), function(x) paste(x[1],x[2],x[3],sep = "-"))
aa_tcga <- aa_tcga[,c(25,2,4)]
# merge PCAWG and TCGA AA results, and duplicate different amplicon results in same patients with priority(ecDNA>BFB>Complex>Linear)
aa_merge <- rbind(aa_tcga,aa_pcawg)
aa_merge$Feature <- NA
unique_pt <- unique(aa_merge$patient_barcode)
for (pt in unique_pt) {
  tmp <- aa_merge[aa_merge$patient_barcode == pt, ]
  if ('ecDNA' %in% tmp$Classification) {
    aa_merge[aa_merge$patient_barcode == pt, 4] <- 'ecDNA'
  } else if ('BFB' %in% tmp$Classification) {
    aa_merge[aa_merge$patient_barcode == pt, 4] <- 'BFB'
  } else if ('Complex non-cyclic' %in% tmp$Classification) {
    aa_merge[aa_merge$patient_barcode == pt, 4] <- 'Complex'
  } else if ('Linear amplification' %in% tmp$Classification) {
    aa_merge[aa_merge$patient_barcode == pt, 4] <- 'Linear'
  }
}
aa_merge <- aa_merge[,c(1,2,4)]
aa_merge <- aa_merge %>%
  mutate(Cohort = case_when(
    grepl('TCGA',patient_barcode) ~ 'TCGA',
    .default = 'PCAWG'
  ))
# Annote tumor type with on the samples
aa_merge$Type <- NA
aa_merge_t <- aa_merge[grepl('TCGA',aa_merge$patient_barcode),]
aa_merge_p <- aa_merge[!grepl('TCGA',aa_merge$patient_barcode),]
matched1 <- match(aa_merge_t$patient_barcode,tcga_meta$patient_barcode)
aa_merge_t$Type <- tcga_meta$project[matched1]
matched2 <- match(aa_merge_p$patient_barcode,pcawg_meta$patient_barcode)
aa_merge_p$Type <- pcawg_meta$project[matched2]
aa_merge_t_dup <- aa_merge_t %>%
  distinct(patient_barcode, Feature, .keep_all = TRUE)
aa_merge_p_dup <- aa_merge_p %>%
  distinct(patient_barcode, Feature, .keep_all = TRUE)

# caculate Feature cases counts in every tumor type with TCGA and PCAWG 
tcga_res <- as.data.frame(table(aa_merge_t_dup$Feature,aa_merge_t_dup$Type))
pcawg_res <- as.data.frame(table(aa_merge_p_dup$Feature,aa_merge_p_dup$Type))
colnames(tcga_res) <- c('Feature','Type','Count')
colnames(pcawg_res) <- c('Feature','Type','Count')

tongji_res <- as.data.frame(read_xlsx('TCGA-PCAWG-ecDNAFreq.xlsx',sheet = 1))
tongji_tcga_res <- tongji_res[tongji_res$Project=='TCGA',]
tongji_pcawg_res <- tongji_res[tongji_res$Project=='PCAWG',]
for (i in 1:nrow(tongji_tcga_res)) {
  type <- tongji_tcga_res[i,2]
  tmp <- tcga_res[tcga_res$Type==type,]
  if (nrow(tmp)>0) {
    tongji_tcga_res[i,7] <- tmp[tmp$Feature=='ecDNA',3]
    tongji_tcga_res[i,8] <- tmp[tmp$Feature=='BFB',3]
    tongji_tcga_res[i,9] <- tmp[tmp$Feature=='Complex',3]
    tongji_tcga_res[i,10] <- tmp[tmp$Feature=='Linear',3]
  }
}
for (i in 1:nrow(tongji_pcawg_res)) {
  type <- tongji_pcawg_res[i,3]
  tmp <- pcawg_res[pcawg_res$Type==type,]
  if (nrow(tmp)>0) {
    tongji_pcawg_res[i,7] <- tmp[tmp$Feature=='ecDNA',3]
    tongji_pcawg_res[i,8] <- tmp[tmp$Feature=='BFB',3]
    tongji_pcawg_res[i,9] <- tmp[tmp$Feature=='Complex',3]
    tongji_pcawg_res[i,10] <- tmp[tmp$Feature=='Linear',3]
  }
}
tongji_res <- rbind(tongji_tcga_res,tongji_pcawg_res)
tongji_res[is.na(tongji_res)] <- 0

# Plot
tongji_res <- as.data.frame(read_xlsx('TCGA-PCAWG-ecDNAFreq.xlsx',sheet = 3))
fin_res <- tongji_res %>%
  group_by(Tumor_Type) %>%
  summarize(No.Patients = sum(No.Patients), 
            ecDNA = sum(ecDNA),
            BFB = sum(BFB),
            Complex = sum(Complex),
            Linear = sum(Linear),
            `Non-fSCNA` = sum(`Non-fSCNA`))

# Figure S1a
plot_data <- fin_res %>%
  gather(Index,Value,ecDNA,BFB,Complex,Linear,`Non-fSCNA`)
plot_data$Freq <- plot_data$Value/plot_data$No.Patients*100
plot_data$Index <- factor(plot_data$Index,levels = rev(c('ecDNA','BFB','Complex','Linear','Non-fSCNA')))
plot_data$Newname <- paste0(plot_data$Tumor_Type," (",plot_data$No.Patients,")")
plot_data_order <- plot_data[order(-plot_data$Freq[plot_data$Index=='ecDNA']), ]

tumor_types <- c("Glioblastoma (51)", "Sarcoma (39)", "Lung SCC (50)", "In-House Urothelial carcinoma (493)",
                 "BIG Urothelial carcinoma (102)", "Esophageal AdenoCA (99)", "Breast cancer (237)",
                 "Ovarian cancer (167)", "Biliary Tract cancer (12)", "Bone cancer (73)", "Skin Melanoma (207)",
                 "Pancreatic AdenoCA (243)", "Cervical SCC (69)", "Endometrial cancer (139)",
                 "Head and neck SCC (143)", "Hepatocellular carcinoma (346)", "Gastric cancer (170)",
                 "Urothelial carcinoma (136)", "Malignant Lymphoma (108)", "Pediatric Brain cancer (248)",
                 "Lower grade glioma (89)", "Lung AdenoCA (147)", "Colorectal AdenoCA (74)",
                 "Pancreatic endocrine neoplasms (86)", "Kidney cancer (227)", "Prostate cancer (340)",
                 "In-House Kidney cancer (122)", "Oral cancer (13)","Acute Myeloid Leukemia (92)", "Chronic Lymphocytic Leukemia (97)",
                 "Esophageal SCC (51)", "Papillary thyroid carcinoma (137)", "Uveal Melanoma (51)")



p <- ggplot(plot_data,aes(x=factor(Newname,levels = tumor_types),y=Freq,fill=Index))+
  geom_bar(stat = 'identity',position = 'stack',color='black',linewidth = 0.1)+
  scale_fill_manual(values = c('ecDNA'='#c23831','BFB'='#66ac55','Complex'='#efaf6f','Linear'='#7089bf','Non-fSCNA'='#fcfdf5'))+
  theme_classic()+
  labs(x='Cancer Type',y='Frequency')+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(angle=90,hjust = 1, vjust = 0.5))+
  coord_flip()

#dev.off()
pdf('Figure.pdf',width = 8,height = 10)
print(p)
dev.off()

# Figure S1b (Bladder cancer Plot) ----
bladder <- c("In-House Urothelial carcinoma (493)","BIG Urothelial carcinoma (102)","Urothelial carcinoma (136)")
plot_data1 <- plot_data[plot_data$Newname%in%bladder,]
plot_data1$Newname <- factor(plot_data1$Newname,levels = bladder)
p1 <- ggplot(plot_data1,aes(x=factor(Newname,levels = c('Urothelial carcinoma (136)','BIG Urothelial carcinoma (102)','In-House Urothelial carcinoma (493)')),y=Freq,fill=Index))+
  geom_bar(stat = 'identity',position = 'stack',color='black',linewidth = 0.1)+
  scale_fill_manual(values = c('ecDNA'='#c23831','BFB'='#66ac55','Complex'='#efaf6f','Linear'='#7089bf','Non-fSCNA'='#fcfdf5'))+
  theme_classic()+
  labs(x='Cancer Type',y='Frequency')+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(angle=90,hjust = 1, vjust = 0.5))

pdf('~/Desktop/ecDNA-revision-Figure/FigureS1b.pdf',width = 4,height = 6)
print(p1)
dev.off()
# pvalue
contingency_table <- matrix(c(186,307,30,72), ncol = 2, byrow = TRUE)
rownames(contingency_table) <- c("G1", "G2")
colnames(contingency_table) <- c("Pos", "Neg")
contingency_table
fisher_result <- fisher.test(contingency_table)
print(fisher_result)



