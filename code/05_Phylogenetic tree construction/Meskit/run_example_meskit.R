library(MesKit)
library(data.table)
library(dplyr)
maf <- readMaf(mafFile = 'maf.tsv',
                         clinicalFile  = 'cli.tsv',
                         refBuild = "hg38",
               nonSyn.vc = c("Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","Translation_Start_Site","Nonsense_Mutation","Nonstop_Mutation","In_Frame_Del","In_Frame_Ins","Missense_Mutation")) 
phyloTree_nonsilent_NJ <- getPhyloTree(maf, patient.id = "prx", method = "ML")
phylotree_nonsilent_NJ <- plotPhyloTree(phyloTree_nonsilent_NJ, use.tumorSampleLabel = TRUE)

