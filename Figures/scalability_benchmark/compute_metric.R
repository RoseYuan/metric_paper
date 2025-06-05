suppressPackageStartupMessages({
library(optparse)
library(poem)
library(spARI)
})


option_list <- list(
    # parameters for preparing clustering
	make_option(c("-i", "--input"), type="character", default=NA, help="input rds file path."),
	make_option(c("-o", "--output"), type="character", default=NA, help="output folder path."),
	make_option(c("-m", "--metric"), type="character", default=NA, help="Metric name."),
    make_option(c("-s", "--size_id"), type="double", default=1, help="Which size to use. 1,2,3,4."),
    make_option(c("-l", "--level"), type="character", default=NA, help="Which level to calculate the metrics. 'dataset','class', or 'element'."),
    make_option(c("-t", "--true"), type="character", default=NA, help="Column name of true labels."),
    make_option(c("-q", "--pred"), type="character", default=NA, help="File path of predicted embeddings."),
    make_option(c("-p", "--pred_col"), type="character", default=NA, help="Column name of predicted labels."),
    make_option(c("-n", "--n_embed"), type="double", default=NA, help="Number of dimensions for the embedding."),
    make_option(c("-k", "--n_cluster"), type="double", default=NA, help="Number of clusters."),
    make_option(c("-a", "--args_id"), type="character", default=NA, help="Argument id."),
    make_option("--extra_args", type="character", default="", help="Extra arguments")
)
# -h should be preserved for --help!!!

opt <- parse_args(OptionParser(option_list=option_list))

df_ls <- readRDS(opt$input)
df <- df_ls[[opt$size_id]]
# parse the extra args string
extra_args_vec <- strsplit(opt$extra_args, " ")[[1]]
# convert to named list
parse_extra_args <- function(arg_vec) {
  args <- list()
  for (kv in arg_vec) {
    parts <- strsplit(kv, "=")[[1]]
    key <- parts[1]
    val <- parts[2]

    # Convert val to appropriate type
    if (!is.na(suppressWarnings(as.numeric(val)))) {
      val <- as.numeric(val)
    } else if (tolower(val) %in% c("true", "false")) {
      val <- as.logical(tolower(val))
    }

    args[[key]] <- val
  }
  return(args)
}
extra_args <- parse_extra_args(extra_args_vec)

if(file.exists(opt$pred)){
    embed <- read.table(opt$pred, header=TRUE, sep="\t")
    embed <- embed[embed$X %in% df$id, paste0("X",0:9)]
    if(opt$metric %in% c("AMSP", "NCE", "NP", "PWC", "cohesion", "adhesion")){
        res <- do.call(getGraphMetrics, c(list(x=embed, labels=df[, opt$true], metrics=opt$metric,level=opt$level), extra_args))
        # res <- getGraphMetrics(x=embed, labels=df[, opt$true], metrics=opt$metric,level=opt$level)
    }else if (opt$metric %in% c("cdbw", "compactness", "dbcv", "meanSW")){
        res <- do.call(getEmbeddingMetrics, c(list(x=embed, labels=df[, opt$true], metrics=opt$metric, level=opt$level), extra_args))
        # res <- getEmbeddingMetrics(x=embed, labels=df[, opt$true], metrics=opt$metric, level=opt$level)
        }else{
        stop("Invalid metric name. Please check the metric name.")
    }
}else{
    if (opt$metric %in% c( "ASPC", "SPC", "AWC", "AWH", "ARI","NCR")){
        res <- do.call(getPartitionMetrics, c(list(true=df[, opt$true], pred=df[, opt$pred_col], metrics=opt$metric,level=opt$level), extra_args))
        # res <- getPartitionMetrics(true=df[, opt$true], pred=df[, opt$pred_col], metrics=opt$metric,level=opt$level)
    }else if (opt$metric %in% c("SpatialARI", "SpatialRI", "SpatialSPC", "SpatialAccuracy", "nsARI", "nsAWC", "nsAWH", "nsSPC")){
        res <- do.call(getSpatialExternalMetrics, c(list(true=df[, opt$true], pred=df[, opt$pred_col], location=df[,c("x","y")], metrics=opt$metric,level=opt$level), extra_args))
        # res <- getSpatialExternalMetrics(true=df[, opt$true], pred=df[, opt$pred_col], location=df[,c("x","y")], metrics=opt$metric,level=opt$level)
    }else if(opt$metric == "SpatialARI_Yan"){
        res <- do.call(spARI, c(list(r_labels=df[, opt$true], c_labels=df[, opt$pred_col], coords=df[,c("x","y")]), extra_args))
    }else{
        stop("Invalid metric name. Please check the metric name.")
    }
}

outfile <- paste0(opt$output, "/", opt$metric, "_", opt$args_id, ".result")
write.table(res, file=outfile, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
