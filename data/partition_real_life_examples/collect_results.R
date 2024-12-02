
sobj_file <- '../../partition_real_life_examples/PBMC_SnapATAC_ndim15_sobj_SNN.RDS'
clustering_file <- '../../partition_real_life_examples/PBMC_SnapATAC_ndim15_seed0_r0.15.tsv'

sobj1 <- readRDS(sobj_file)
df_label <- read.table(clustering_file, header = T, sep="\t", comment.char = "")
sobj2 <- add_labels(sobj1, df_label, "barcode", "clusterings")

ground_truth <- sobj1$ground_truth
ground_truth <- factor(gsub("_", " ", ground_truth))
clusterings <- sobj2$ground_truth

df_ground_truth <- data.frame(ground_truth=ground_truth, barcode=names(clusterings))
write.table(df_ground_truth, "PBMC_ground_truth.tsv", sep='\t', quote=FALSE, row.names=FALSE)