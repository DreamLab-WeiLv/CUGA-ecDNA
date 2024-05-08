 # library R package
library(readxl)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(bedtoolsr)
library(ggplot2)
library(ggsci)
library(ggrepel)
# get ref chromosome region
df <- data.frame(chromName = seqnames(BSgenome.Hsapiens.UCSC.hg38), 
                 chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38))
df$chromNum <- 1:length(df$chromName)
df <- df[1:22,] 
df$chromlengthCumsum <- cumsum(as.numeric(df$chromlength)) 
df$chormStartPosFrom0 <- c(0,df$chromlengthCumsum[-nrow(df)])
tmp_middle <- diff(c(0,df$chromlengthCumsum)) / 2
df$chromMidelePosFrom0 <- df$chormStartPosFrom0 + tmp_middle
# Input segment files
d <- as.data.frame(read_xlsx('All_tumor_samples_segcn.xlsx',sheet = 1))
d <- na.omit(d)
filter_chromo <- c("chrX","chrY")
d <- d %>%
  filter(!d$V1%in%filter_chromo)
d$V1 <- gsub('chr','',d$V1)
d$V1 <- as.numeric(d$V1)
chromID <- d$V1
d$StartPos <- d$V2 + df$chormStartPosFrom0[chromID]
d$EndPos <- d$V3 + df$chormStartPosFrom0[chromID]
# get ecDNA and other fSCNA matrix
decDNA <- d %>%
  filter(grepl('ecDNA',d$Feature))
table(decDNA$Feature)
dother <- d %>%
  filter(!grepl('ecDNA',d$Feature)&!Feature=='unknown_1')
table(dother$Feature)
colnames(decDNA)[1] <- "chromosome"
colnames(dother)[1] <- "chromosome"
# define bin size = 300kb
segment_length <- 300000  
segment_data_ec <- data.frame(
  chr = character(0),
  segment_start = numeric(0),
  segment_end = numeric(0))
for (i in seq_len(nrow(df))) {
  chr <- df$chromName[i]
  startPos <- df$chormStartPosFrom0[i]
  midPos <- df$chromlengthCumsum[i]
  segment_breaks <- seq(startPos, midPos, by = segment_length)
  for (j in seq_along(segment_breaks[-length(segment_breaks)])) {
    segment_start <- segment_breaks[j]
    segment_end <- segment_breaks[j + 1] - 1
    segment_data_ec <- rbind(segment_data_ec, c(chr, segment_start, segment_end))
  }
}
colnames(segment_data_ec) <- c("chr", "segment_start", "segment_end")
segment_data_ec$segment_start <- as.numeric(segment_data_ec$segment_start)
segment_data_ec$segment_end <- as.numeric(segment_data_ec$segment_end)
segment_data_ec <- na.omit(segment_data_ec)

segment_data_ec$chr <- as.numeric(gsub('chr','',segment_data_ec$chr))
segment_data_ec$cn <- 0
segment_data_ec$Freq <- 0
for (i in 1:nrow(segment_data_ec)) {
  seg_chr <- segment_data_ec$chr[i]
  seg_start <- segment_data_ec$segment_start[i]
  seg_end <- segment_data_ec$segment_end[i]
  # find ecDNA and ref overlap region
  overlapping_rows <- decDNA[(decDNA$chromosome == seg_chr & decDNA$StartPos <= seg_end & decDNA$EndPos >= seg_start), ]
  
  if (nrow(overlapping_rows) > 0) {
    max_cn <- max(overlapping_rows$cn)
    segment_data_ec$cn[i] <- max_cn
    unique_samples <- unique(overlapping_rows$Sample)
    segment_data_ec$Freq[i] <- length(unique_samples)
  }
}
# Freq = Freq number/total ecDNA+ patient
length(table(decDNA$CaseID)) #216 
segment_data_ec$Freq <- segment_data_ec$Freq/216
segment_data_ec$cn <- as.numeric(segment_data_ec$cn)
# Point Plot 
pdf('Figure_ecDNA_jitter.pdf',width=20,height = 5)
p1 <- ggplot(segment_data_ec, aes(x = segment_start, y = cn)) +
  geom_jitter(size = 0.8) +
  scale_fill_lancet(guide=guide_legend(reverse = T)) +
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=3)+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),labels = NULL)+
  scale_y_continuous(limits = c(0.00001,370),breaks = c(0,100,200,370),labels = c("0","100","200",">300"))+
  xlab('')+
  theme_classic()+
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank())
print(p1)
dev.off()

# Freq Plot 
pdf('Figure_ecDNA_Freq.pdf',width=30,height = 10)
p2 <- ggplot(segment_data_ec,aes(x=segment_start,y=Freq)) +
  geom_segment(aes(x = segment_start,xend=segment_end,y=0,yend=Freq),color="#d24569",position = "identity") + 
  scale_fill_lancet(guide=guide_legend(reverse = T))+
  geom_vline(data = df, mapping = aes(xintercept=chromlengthCumsum),linetype=3)+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),labels = NULL)+
  scale_y_continuous(limits = c(-0.01,0.2))+
  theme_classic()+
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank())+
  geom_path(colour='#d24569',linewidth=1)+
  xlab('')
print(p2)
dev.off()

# Other fSCNAs Plot
segment_length <- 300000  
segment_data_other <- data.frame(
  chr = character(0),
  segment_start = numeric(0),
  segment_end = numeric(0))

for (m in seq_len(nrow(df))) {
  chr <- df$chromName[m]
  startPos <- df$chormStartPosFrom0[m]
  midPos <- df$chromlengthCumsum[m]
  segment_breaks <- seq(startPos, midPos, by = segment_length)
  for (n in seq_along(segment_breaks[-length(segment_breaks)])) {
    segment_start <- segment_breaks[n]
    segment_end <- segment_breaks[n + 1] - 1
    segment_data_other <- rbind(segment_data_other, c(chr, segment_start, segment_end))
  }
}
colnames(segment_data_other) <- c("chr", "segment_start", "segment_end")
segment_data_other$segment_start <- as.numeric(segment_data_other$segment_start)
segment_data_other$segment_end <- as.numeric(segment_data_other$segment_end)
segment_data_other <- na.omit(segment_data_other)
segment_data_other$chr <- as.numeric(gsub('chr','',segment_data_other$chr))
segment_data_other$cn <- 0
segment_data_other$Freq <- 0
for (k in 1:nrow(segment_data_other)) {
  seg_chr <- segment_data_other$chr[k]
  seg_start <- segment_data_other$segment_start[k]
  seg_end <- segment_data_other$segment_end[k]
  overlapping_rows <- dother[(dother$chromosome == seg_chr & dother$StartPos <= seg_end & dother$EndPos >= seg_start), ]
  
  if (nrow(overlapping_rows) > 0) {
    max_cn <- max(overlapping_rows$cn)
    segment_data_other$cn[k] <- max_cn
    unique_samples <- unique(overlapping_rows$Sample)
    segment_data_other$Freq[k] <- length(unique_samples)
  }
}
length(table(dother$CaseID)) #324 
segment_data_other$Freq <- segment_data_other$Freq/324
segment_data_other$cn <- as.numeric(segment_data_other$cn)
# Point Plot
pdf('~/Desktop/ecDNA-revision-Figure/FigureS3a_other_jitter.pdf',width=20,height = 5)
p3 <- ggplot(segment_data_other, aes(x = segment_start, y = cn)) +
  geom_jitter(size = 0.8) +
  scale_fill_lancet(guide=guide_legend(reverse = T)) +
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=3)+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),labels = NULL)+
  scale_y_continuous(limits = c(0.00001,332),breaks = c(0,100,200,332),labels = c("0","100","200",">300"))+
  xlab('')+
  theme_classic()+
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank())
print(p3)
dev.off()


# show gene info
gene_pos <- read.table('protein_coding.hg38.position.txt',sep = '\t',header = F)
gene_pos$V1 <- gsub('chr','',gene_pos$V1)
gene_pos$V1 <- as.numeric(gene_pos$V1)
chromID <- gene_pos$V1
gene_pos$StartPos_gene <- gene_pos$V2 + df$chormStartPosFrom0[chromID]
gene_pos$EndPos_gene <- gene_pos$V3 + df$chormStartPosFrom0[chromID]
gene_pos <- gene_pos[,c(1,5,6,4)]
gene_pos$V1 <- paste0('chr',gene_pos$V1)
gene_pos <- na.omit(gene_pos)
segment_data_ec$chr <- paste0('chr',segment_data_ec$chr)
segment_data_gene <- bt.intersect(a=segment_data_ec,b=gene_pos,wa=T,wb=T)
segment_data_gene <- segment_data_gene[,c(1:3,9)] 
colnames(segment_data_gene) <- c('chr','segment_start','segment_end','gene_name')
show_gene <- c('CCND1','FGF4','FGF3','CTTN','E2F3','YWHAZ','SOX4','MDM2','PRARG',
               'YEATS4','RAF1','TBX3','UBR5','ERBB2')
segment_data_gene_show <- segment_data_gene %>%
  filter(segment_data_gene$gene_name%in%show_gene)

pdf('Figure_other_Freq.pdf',width=20,height = 5)
p4 <- ggplot(segment_data_other,aes(x=segment_start,y=Freq)) +
  geom_segment(aes(x = segment_start,xend=segment_end,y=0,yend=Freq),color="#2f6fb3",position = "identity") + 
  scale_fill_lancet(guide=guide_legend(reverse = T))+
  geom_vline(data = df, mapping = aes(xintercept=chromlengthCumsum),linetype=3)+
  geom_text(data = df, aes(x=chromMidelePosFrom0,y=-0.03,label=chromNum))+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),labels = NULL)+
  scale_y_continuous(limits = c(-0.03,0.2))+
  xlab('')+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank())+
  geom_path(colour='#2f6fb3',linewidth=1) +
  geom_text(data = segment_data_gene_show,
            aes(x = segment_start, y = 0.1, label = gene_name),
            vjust = 1.5, 
            hjust = 0.5,  
            size = 2)

print(p4)
dev.off()

library(grid)
library(gridExtra)
pdf('Figure.pdf',width = 25,height = 6)
p <- grid.arrange(p1,p2,p3,p4,ncol=1)
dev.off()
