
####################
# Metrics
####################

adjusted_wallance_indices <-function(true=NULL, pred=NULL, contigency_res=NULL){  # cannot be calculated when there's singletons
    
    ## get pairs using C
    ## ensure that values of c1 and c2 are between 0 and n1
    if (is.null(contigency_res)){
      res <- aricode::sortPairs(true, pred)
      n <- length(true)
    } else{
      res <- contigency_res
      n <- sum(res$ni.)
    }
    spairs <- n*(n-1)/2 # N
    stot <- sum(choose(res$nij, 2), na.rm=TRUE) # T
    srow <- sum(choose(res$ni., 2), na.rm=TRUE) # P
    scol <- sum(choose(res$n.j, 2), na.rm=TRUE) # Q
    a <- stot
    b <- srow-stot
    c <- scol-stot
    d <- spairs+stot-srow-scol
    aw <- (a*d-b*c)/((a+b)*(b+d))
    av <- (a*d-b*c)/((a+c)*(c+d))
    ari <- 2*(a*d-b*c)/((a+b)*(b+d)+(a+c)*(c+d))
    
    awi <- list()
    avj <- list()
    for (i in sort(unique(res$pair_c1))){
      idx <- which(res$pair_c1 == i)
      term1 <- spairs * sum(choose(res$nij[idx], 2)) 
      term2 <- choose(sum(res$nij[idx]), 2) * scol
      term3 <- choose(sum(res$nij[idx]), 2) * (spairs - scol)
      awi[i+1] <- (term1 - term2) / term3
    }

    for (j in sort(unique(res$pair_c2))){
      idx <- which(res$pair_c2 == j)
      term1 <- spairs * sum(choose(res$nij[idx], 2)) 
      term2 <- choose(sum(res$nij[idx]), 2) * srow
      term3 <- choose(sum(res$nij[idx]), 2) * (spairs - srow)
      avj[j+1] <- (term1 - term2) / term3
    }
    
    # remove NA introduced by singletons
    aw2 <- mean(unlist(awi[!is.na(awi)]))
    av2 <- mean(unlist(avj[!is.na(avj)]))
    ari2 <- 2*aw2*av2/(aw2+av2)
    return(list("AW"=aw, "AV"=av, "ARI"=ari, "AW2"=aw2, "AV2"=av2, "ARI2"=ari2,"Awi"=awi, "Avj" = avj))
}

pwc <- function(graph, label, label_idx){
  suppressPackageStartupMessages({
    require(igraph)
  })
  V(graph)$label <- label
  label_ls <- unique(label)
  df_graph <- as_data_frame(graph, what = c("edges"))
  
  v_1 <- V(graph)[V(graph)$label == label_ls[label_idx]]
  v_2 <- V(graph)[V(graph)$label != label_ls[label_idx]]

  
  l1 <- length(v_1)
  
  j1 <- 0
  w1 <- 0
  for (i in v_1){
      df_tmp <- df_graph[df_graph$to == i | df_graph$from == i,] # edges of node i
      a <- sum(df_tmp$from %in% v_1 & df_tmp$to %in% v_1) # within cluster edges of node i
      b <- sum(df_tmp$from %in% v_1 & df_tmp$to %in% v_2) + sum(df_tmp$to %in% v_1 & df_tmp$from %in% v_2)# between cluster edges of node i
      c <- sum(df_tmp$weight[df_tmp$from %in% v_1 & df_tmp$to %in% v_1])
      d <- sum(df_tmp$weight[df_tmp$from %in% v_1 & df_tmp$to %in% v_2]) + sum(df_tmp$weight[df_tmp$to %in% v_1 & df_tmp$from %in% v_2])
      if(b>=a){j1 <- j1+1}
      if(d>=c){w1 <- w1+1}
  }
  return(list(j1=j1, j1_frac=j1/l1, w1=w1, w1_frac=w1/l1))
}

cal_distance <- function(x, metric="Euclidean") {
        if (metric == "Euclidean"){
        dist.matrix <- stats::dist(x)
    }else if(metric == "correlation"){
        dist.matrix <- amap::Dist(x, method = "correlation")
    }else if(metric == "manhattan"){
        dist.matrix <- stats::dist(x,method = "manhattan")
    }else if(metric == "cosine"){
        dist.matrix <- stylo::dist.cosine(x)
        }
}
# for a given distance matrix and given labels, calculate silhouette score, plot the silhouette distribution and average silhouette.
silhouette_result <- function(dist.matrix, labels, title="", my_col=NULL){
    suppressPackageStartupMessages({
        require(cluster)
        require(pals)
    })
    if (is.null(my_col)){
        my_col <- unlist(polychrome())[1:max(as.numeric(as.factor(labels)))]
    }

    sil = silhouette(as.numeric(as.factor(labels)), dist.matrix)
    p <- plot(sil, border=NA, col=my_col, main=title)
    avg_sil <- mean(sil[, 3])
    return(list("fig" = p, "avg" = avg_sil, "sil" = sil))
}

##################
# Visualization
##################
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(paletteer)
    library(ggpubr)
    library(RColorBrewer)
    library(scales)
    library(pals)
    library(igraph)
    library(dplyr)
})

my_col_paired <- brewer.pal(12, "Paired")
my_col_polychrome <- as.vector(unlist(polychrome()))
my_col_WRdDk <- c("#FFFFFF","#FFFF66","#FFCC33","#F57F17","#EE0000","#990000","#660000")
my_col_WPkRd <- c("#FFFFFF","#FFCCCC","#FF9999","#FF6666","#FF3333","#CC0033")
my_col_BlWPp <- c("#0066CC","#3399FF","#FFFFFF","#FF99FF","#FF00FF")


cross_table_plot <- function(ground_truth, clusterings, a=1.3, b=5.7, c=2, m=0, n=0.2){
    x <- unique(ground_truth)
    y <- as.factor(unique(clusterings))
    data <- expand.grid(X=x, Y=y)
    cross_count <- table(ground_truth, clusterings) 

    # cell count in the cross_count table
    data$Z1 <- apply(data, 1, function(x){cross_count[x[["X"]],as.numeric(x[["Y"]])]})
    # log transform Z1
    data$Z2 <- apply(data, 1, function(x){log(cross_count[x[["X"]],as.numeric(x[["Y"]])]+1)})
    # row normalize of Z1
    data <- data %>% group_by(X) %>% mutate(Z3 = 100*Z1/sum(Z1))

    top_cluster <- data %>% group_by(X) %>% top_n(1, Z1)

    unselected_Y <- setdiff(unique(data$Y), unique(top_cluster$Y))
    top_cluster <- top_cluster[order(top_cluster$X),]
    new_levels <- as.numeric(c(as.character(unique(top_cluster$Y)),unselected_Y))
    data$Y <- factor(data$Y, levels=new_levels)

    res <- adjusted_wallance_indices(ground_truth, clusterings)

    df_awi <- do.call(rbind, Map(data.frame, "awi"=res$Awi))
    df_awi$cell_type <- levels(ground_truth)
    df_count <- data.frame(table(ground_truth))
    df_count$Freq <- round(df_count$Freq / sum(df_count$Freq),3)
    colnames(df_count) <- c("cell_type", "frac")
    df_awi <- merge(df_awi, df_count, by=c("cell_type"))

    df_avj <- do.call(rbind, Map(data.frame, "avj"=res$Avj))
    df_avj$cell_type <- levels(clusterings)
    df_avj$cell_type <- factor(df_avj$cell_type, levels=new_levels)
    df_count <- data.frame(table(clusterings))
    df_count$Freq <- round(df_count$Freq / sum(df_count$Freq),3)
    colnames(df_count) <- c("cell_type", "frac")
    df_avj <- merge(df_avj, df_count, by=c("cell_type"))

    main <- ggplot(data, aes(Y, X, fill= Z3)) + 
    geom_tile(colour="black", size=0.1) + 
    scale_fill_gradientn(colours = my_col_WRdDk, breaks=seq(0,100,10), guide = guide_colourbar(title.position = "top")) +
    labs(x="ATAC cluster", y="True class", fill = "Cells per true class %") + 
    theme(legend.direction = "horizontal", 
        legend.position = "bottom", 
        legend.key.width= unit(a, 'cm'),
        legend.text=element_text(size=16, face="bold"),
        legend.title=element_text(size=16, face="bold")) + 
    theme(plot.margin = unit(c(0, 0, 0, 1), "cm"),
        panel.background = element_blank(),
        axis.title.x = element_text(size = 16, margin = margin(5,0,0,0), face="bold"), 
        axis.text.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, margin = margin(0,0,0,0), face="bold"), 
        axis.text.y = element_text(size = 16, face="bold"))

    bp.x <- ggplot(data = df_avj, aes(x = cell_type, y = avj, fill=frac)) + 
    geom_bar(stat = "identity", colour = "black", size=0.1) + 
    scale_fill_gradientn(colours = my_col_WPkRd, guide = guide_colourbar(), limits = c(m, n)) +
    theme(
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size=14, face="bold"), 
        axis.title.y = element_text(size = 16, margin = margin(10,5,0,0), face="bold"),
        legend.position = "none",
        # legend.position = "top",
        panel.background = element_blank()) + 
    labs(y = "AWH", fill="Fraction")+ theme(plot.margin = unit(c(1, 0, 0, b), "cm"))

    bp.y <- ggplot(data = df_awi, aes(x = cell_type, y = awi, fill=frac)) +
    geom_bar(stat = "identity", colour = "black", size=0.1) + coord_flip() + # theme_ipsum() + #theme_gray() +
    scale_fill_gradientn(colours = my_col_WPkRd, guide = guide_colourbar(), limits = c(m, n)) +
    theme(axis.title.x = element_text(size = 16, margin = margin(5,5,0,0), face="bold"), 
            axis.text.x = element_text(size = 7, face="bold"),
            axis.text.y = element_blank(), 
            axis.title.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            legend.text=element_text(size=12, face="bold"),
            legend.title=element_text(size=14, face="bold"),
            # legend.position="none",
            panel.background = element_blank()) + 
    labs(y="AWC", fill="Fraction")+ theme(plot.margin = unit(c(0, 0, c, 0), "cm"))

    df_hm <- data.frame(cols = numeric(0), value = numeric(0))

    gg_empty <- df_hm %>% 
    ggplot(aes(x = cols, y = value)) +
    geom_blank() +
    theme(axis.text = element_blank(),
            axis.title = element_blank(),
            line = element_blank(),
            panel.background = element_blank()) + 
            geom_text(size=7,  fontface = "bold", aes(label = paste0("AWC = ", round(res$AW,3), "\n", "AWH = ", round(res$AV,3), "\n", "ARI = ", round(res$ARI,3))), x = 0.5, y = 0.5)

    plot <- ggarrange(
    bp.x, gg_empty, main, bp.y,
    nrow = 2, ncol = 2, widths = c(3, 1.2), heights = c(1, 3)
    )
    return(plot)
}
##################
# Seurat object related
##################

# function to add ground truth to a Signac object, and check the cell id
# df_label: a dataframe containing a barcode column and a label column
add_labels <- function(sobj, df_label, barcode_col, label_col){
    suppressPackageStartupMessages({
        require(Signac)
        require(Seurat)
    })
    cell_id <- Cells(sobj)
    rownames(df_label) <- df_label[, barcode_col]
    cell_labels <- df_label[cell_id, label_col]
    Idents(sobj) <- factor(cell_labels)
    sobj$ground_truth <- factor(cell_labels)
    return(sobj)
}

# function to add an embedding object to a Signac object
# embed: a table where rows are cells and columns are latent dimensions. row names should be the cell id
add_embedding <- function(sobj, embed, embed_name="learned_embedding", max_dim=NULL){
    suppressPackageStartupMessages({
        require(Signac)
        require(Seurat)
    })
    cell_id <- Cells(sobj)

    if(!is.null(max_dim)){
        embed <- embed[,1:max_dim]
    }

    cell_id <- cell_id[toupper(cell_id) %in% toupper(rownames(embed))]
    # take only the intersection of cells
    sobj <- subset(x=sobj, cells=cell_id)
    
    g <- rep(seq_along(cell_id), sapply(cell_id, length))
    embed_id <- g[match(toupper(cell_id), toupper(rownames(embed)))]
    embed <- embed[embed_id,]
    if (all(toupper(rownames(embed)) == toupper(Cells(sobj)))){
        rownames(embed) <- Cells(sobj)
    } else(stop("All cells in the embedding being added must match the cells in the object!"))

    colnames(embed) <- NULL

    sobj@reductions[[embed_name]] <- CreateDimReducObject(embeddings = as.matrix(embed), key = "LSI_", assay = DefaultAssay(sobj))
    return(sobj)
}

graph_from_sobj <- function(sobj, embedding_name="learned_embedding", n_neighbors=20){
    sobj <- Seurat::FindNeighbors(object = sobj, 
                        reduction = embedding_name, 
                        graph.name = c(paste0("knn_k", n_neighbors), paste0("snn_k", n_neighbors)),
                        dims=1:15,
                        k.param=n_neighbors)
    gm <- sobj@graphs[[paste0("snn_k", n_neighbors)]]
    attributes(gm)$class <- "dgCMatrix"
    g <- igraph::graph_from_adjacency_matrix(adjmatrix = gm, add.colnames = TRUE, mode = "directed", weighted = TRUE)
    return(g)
}