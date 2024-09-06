dataset_ls <- c("Cell_line_mixing", "candidate1",  "candidate2", "Buenrostro_2018",  "Chen_2019",  "PBMC_multiomics") 
dataset_name <- c("Cell line mixing experiment", "human adult atlas subset1", "human adult atlas subset2",  "Buenrostro2018", "Chen2019", "10X PBMC multiomics")
k_ls <- c(10, 13, 10, 9, 13, 15)

# Put meta files of all datasets together
for(j in 1:length(dataset_ls)){
    dataset <- dataset_ls[j]
    k_optimal <- k_ls[j]
    if(j == 1){
        df <- read.table(file=paste0("/home/siluo/public/SiyuanLuo/projects/benchmark/outputs/", dataset, "/meta_info_for_evaluation.tsv"), sep="\t", header=TRUE)
        df_metrics <- read.table(file=paste0("/home/siluo/public/SiyuanLuo/projects/benchmark/outputs/", dataset, "/metrics_value.tsv"), sep="\t", header=TRUE)
        df$dataset <- dataset_name[j]
        df$dataset2 <- dataset_ls[j]
        df_metrics$dataset <- dataset_name[j]
        df$k_optimal <- k_optimal
        df_metrics$k_optimal <- k_optimal
        df_metrics$dataset2 <- dataset_ls[j]
    }else{
        df1 <- read.table(file=paste0("/home/siluo/public/SiyuanLuo/projects/benchmark/outputs/", dataset, "/meta_info_for_evaluation.tsv"), sep="\t", header=TRUE)
        df1$dataset <- dataset_name[j]
        df1$k_optimal <- k_optimal
        df1$dataset2 <- dataset_ls[j]
        
        df <- rbind(df, df1)
        df_metrics1 <- read.table(file=paste0("/home/siluo/public/SiyuanLuo/projects/benchmark/outputs/", dataset, "/metrics_value.tsv"), sep="\t", header=TRUE)
        df_metrics1$dataset <- dataset_name[j]
        df_metrics1$k_optimal <- k_optimal
        df_metrics1$dataset2 <- dataset_ls[j]
        df_metrics <- rbind(df_metrics, df_metrics1)
    }
}

df_seed0 <- df
df_metrics_seed0 <- df_metrics
df_seed0$seed <- 0
df_metrics_seed0$seed <- 0

# Put meta files of all datasets together
for(j in 1:length(dataset_ls)){
    dataset <- dataset_ls[j]
    k_optimal <- k_ls[j]
    if(j == 1){
        df <- read.table(file=paste0("/home/siluo/public/SiyuanLuo/projects/benchmark/outputs/", dataset, "/meta_info_for_evaluation_seed.tsv"), sep="\t", header=TRUE)
        df_metrics <- read.table(file=paste0("/home/siluo/public/SiyuanLuo/projects/benchmark/outputs/", dataset, "/metrics_value_seed.tsv"), sep="\t", header=TRUE)
        df$dataset <- dataset_name[j]
        df$dataset2 <- dataset_ls[j]
        df_metrics$dataset <- dataset_name[j]
        df$k_optimal <- k_optimal
        df_metrics$k_optimal <- k_optimal
        df_metrics$dataset2 <- dataset_ls[j]
    }else{
        df1 <- read.table(file=paste0("/home/siluo/public/SiyuanLuo/projects/benchmark/outputs/", dataset, "/meta_info_for_evaluation_seed.tsv"), sep="\t", header=TRUE)
        df1$dataset <- dataset_name[j]
        df1$k_optimal <- k_optimal
        df1$dataset2 <- dataset_ls[j]
        
        df <- rbind(df, df1)
        df_metrics1 <- read.table(file=paste0("/home/siluo/public/SiyuanLuo/projects/benchmark/outputs/", dataset, "/metrics_value_seed.tsv"), sep="\t", header=TRUE)
        df_metrics1$dataset <- dataset_name[j]
        df_metrics1$k_optimal <- k_optimal
        df_metrics1$dataset2 <- dataset_ls[j]
        df_metrics <- rbind(df_metrics, df_metrics1)
    }
}


df <- rbind(df, df_seed0)
df <- df[!duplicated(df), ]

df_metrics <- rbind(df_metrics, df_metrics_seed0)
df_metrics <- df_metrics[!duplicated(df_metrics), ]

df_metrics <- subset(df_metrics, select = -c(dataset2, dataset, clustering_file, rds_file, snn_file))
names(df_metrics)[names(df_metrics) == "method"] <- "tool"
names(df_metrics)[names(df_metrics) == "long_method"] <- "method"
names(df_metrics)[names(df_metrics) == "dataset_short"] <- "dataset"
head(df_metrics)

write.table(df_metrics, "results_all_metrics.tsv", sep='\t', quote=FALSE, row.names=FALSE)