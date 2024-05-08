library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg38)
ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
type_occurences <- mut_type_occurrences(vcfs, ref_genome)
p <- plot_spectrum(type_occurences)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
zero_cols <- colSums(mut_mat == 0) == nrow(mut_mat)
mut_mat <- mut_mat[, !zero_cols]
sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
cancer_signatures = as.matrix(cancer_signatures[,4:33])
hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")
cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
plot(hclust_cosmic)
cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)
cos_sim_samples_signatures <- na.omit(cos_sim_samples_signatures)
plot_cosine_heatmap(cos_sim_samples_signatures,
                    cluster_rows = F)
fit_res <-fit_to_signatures(mut_mat, cancer_signatures)
select <- which(rowSums(fit_res$contribution) > 0)
plot_contribution(fit_res$contribution[select,],cancer_signatures[,select],coord_flip = FALSE,mode = "absolute")
plot_contribution_heatmap(fit_res$contribution,cluster_samples = T,method = "complete")
res <- fit_res[["contribution"]]
res <- apply(res, 2, function(x){
  prop.table(x)})

