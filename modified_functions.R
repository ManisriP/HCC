#' Helper to handle SummarizedExperiment or expression data frame as input
#'
#' @param exp Expression data as a data frame or a SummarizedExperiment object
#'
#' @return If exp is a SummarizedExperiment object, it will return `assay(exp)`.
#' Otherwhise, it will simply return exp as it is.
#' @noRd
#' @importFrom SummarizedExperiment assay
handleSE <- function(exp) {
  if(is(exp, "SummarizedExperiment")) {
    fexp <- SummarizedExperiment::assay(exp)
  } else {
    fexp <- exp
  }
  return(fexp)
}

#' Replace content of a SummarizedExperiment object based on filtered expression data frame
#'
#' @param exp Expression data frame with genes IDs in row names and samples in column names.
#' @param SE Original SummarizedExperiment object.
#' @return A SummarizedExperiment object
#' @noRd
#' @importFrom SummarizedExperiment SummarizedExperiment colData assay
exp2SE <- function(exp, SE) {
  
  # Get overlaps between samples in both sets
  SE_cols <- rownames(SummarizedExperiment::colData(SE))
  filt_exp_cols <- colnames(exp)
  overlap <- which(SE_cols %in% filt_exp_cols)
  
  # Modify original SE based on filtered expression data frame
  SE_final <- SummarizedExperiment::SummarizedExperiment(
    assays = exp,
    colData = SummarizedExperiment::colData(SE)[overlap, , drop=FALSE]
  )
  
  return(SE_final)
}


#' Generate custom color palette
#'
#' @param pal Numeric specifying palette number, from 1 to 3.
#'
#' @return Character vector of custom color palette with 20 colors
#' @noRd
custom_palette <- function(pal = 1) {
  pal1 <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF",
            "#9467BDFF", "#8C564BFF", "#E377C2FF", "#7F7F7FFF",
            "#BCBD22FF", "#17BECFFF", "#AEC7E8FF", "#FFBB78FF",
            "#98DF8AFF", "#FF9896FF", "#C5B0D5FF", "#C49C94FF",
            "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF")
  
  pal2 <- c("#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF",
            "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF",
            "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF",
            "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", "#C6DBEFFF",
            "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF")
  
  pal3 <- c("#393B79FF", "#637939FF", "#8C6D31FF", "#843C39FF",
            "#7B4173FF", "#5254A3FF", "#8CA252FF", "#BD9E39FF",
            "#AD494AFF", "#A55194FF", "#6B6ECFFF", "#B5CF6BFF",
            "#E7BA52FF", "#D6616BFF", "#CE6DBDFF", "#9C9EDEFF",
            "#CEDB9CFF", "#E7CB94FF", "#E7969CFF", "#DE9ED6FF")
  
  l <- list(pal1, pal2, pal3)
  l_final <- l[[pal]]
  return(l_final)
}


#' Wrapper to handle sample color annotation in heatmap
#'
#' @param col_metadata Sample annotation.
#' @param fexp Gene expression data frame
#'
#' @return List containing processed col_metadata, fexp and annotation_color.
#' @noRd
sample_cols_heatmap <- function(col_metadata, fexp) {
  col_names <- c("Sample group 1", "Sample group 2")
  if(is.data.frame(col_metadata)) {
    colnames(col_metadata) <- col_names[seq_along(col_metadata)]
    col_metadata <- col_metadata[order(col_metadata[, 1]), , drop=FALSE]
    fexp <- fexp[, rownames(col_metadata)]
    if(ncol(col_metadata) == 1) {
      colors1 <- custom_palette(1)[seq_along(unique(col_metadata[,1]))]
      annotation_color <- list(`Sample group 1` = colors1)
      names(annotation_color$`Sample group 1`) <- unique(col_metadata[,1])
    } else if(ncol(col_metadata) == 2) {
      colors1 <- custom_palette(1)[seq_along(unique(col_metadata[,1]))]
      colors2 <- custom_palette(2)[seq_along(unique(col_metadata[,2]))]
      annotation_color <- list(`Sample group 1` = colors1,
                               `Sample group 2` = colors2)
      names(annotation_color$`Sample group 1`) <- unique(col_metadata[,1])
      names(annotation_color$`Sample group 2`) <- unique(col_metadata[,2])
    } else {
      stop("Maximum number of columns for col_metadata is 2.")
    }
  }
  
  if(!exists("annotation_color")) {
    annotation_color <- NA
  }
  results <- list(col_metadata = col_metadata,
                  fexp = fexp,
                  annotation_colors = annotation_color)
  return(results)
  
}


#' Wrapper to handle gene color annotation in heatmap
#'
#' @param row_metadata Gene annotation.
#' @param fexp Gene expression data frame.
#' @param annotation_color Object returned by \code{sample_cols_heatmap}.
#'
#' @return List containing processed row_metadata, fexp and annotation_color.
#' @noRd
gene_cols_heatmap <- function(row_metadata, fexp, annotation_color) {
  if(is.data.frame(row_metadata)) {
    colnames(row_metadata) <- "Gene annotation"
    row_metadata <- row_metadata[order(row_metadata[, 1]), , drop=FALSE]
    fexp <- fexp[rownames(row_metadata), ]
    if(!is.list(annotation_color)) {
      annotation_color <- list()
    }
    annotation_color$`Gene annotation` <- custom_palette(3)[
      seq_along(unique(row_metadata[,1]))
    ]
    names(annotation_color$`Gene annotation`) <- unique(row_metadata[,1])
  }
  
  results <- list(row_metadata = row_metadata,
                  fexp = fexp,
                  annotation_color = annotation_color)
  return(results)
}


#' Set theme for expression profile plot
#'
#' @return Custom theme for the function \code{plot_expression_profile}.
#' @noRd
#' @importFrom ggplot2 theme element_text element_blank element_rect
theme_exp_profile <- function() {
  theme <- ggplot2::theme(
    plot.title = ggplot2::element_text(lineheight=0.8,
                                       face='bold',
                                       colour='black',
                                       size=13,
                                       hjust=0.5),
    axis.title = ggplot2::element_text(size=11),
    axis.text.y = ggplot2::element_text(angle=0,
                                        vjust=0.5,
                                        size=8),
    axis.text.x = ggplot2::element_text(angle=90,
                                        vjust=0.5,
                                        size=6),
    panel.grid = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 8),
    legend.background = ggplot2::element_rect(fill='gray90',
                                              size=0.5,
                                              linetype='dotted'),
    legend.position='bottom'
  )
  return(theme)
}


#' Wrapper to calculate correlation matrix and adjacency matrices
#'
#' @param cor_method Correlation method.
#' @param norm.exp Gene expression data frame.
#' @param SFTpower SFT power calculated with \code{SFT_fit}.
#' @param net_type Network type. One of 'signed', 'signed hybrid' or 'unsigned'.
#'
#' @return List containing the correlation matrix and the adjacency matrix.
#' @noRd
#' @importFrom WGCNA adjacency.fromSimilarity bicor
calculate_cor_adj <- function(cor_method, norm.exp, SFTpower,
                              net_type) {
  
  if(cor_method == "pearson") {
    cor_matrix <- cor(t(norm.exp), method = "pearson")
    adj_matrix <- WGCNA::adjacency.fromSimilarity(cor_matrix,
                                                  power = SFTpower,
                                                  type=net_type)
  } else if(cor_method == "spearman") {
    cor_matrix <- cor(t(norm.exp), use="p", method = "spearman")
    adj_matrix <- WGCNA::adjacency.fromSimilarity(cor_matrix,
                                                  power=SFTpower,
                                                  type=net_type)
  } else if (cor_method == "biweight") {
    cor_matrix <- WGCNA::bicor(t(norm.exp), maxPOutliers = 0.1)
    adj_matrix <- WGCNA::adjacency.fromSimilarity(cor_matrix,
                                                  power=SFTpower,
                                                  type=net_type)
  } else {
    stop("Please, specify a correlation method. One of 'spearman', 'pearson' or 'biweight'.")
  }
  results <- list(cor = cor_matrix,
                  adj = adj_matrix)
  return(results)
}


#' Wrapper to assign TOM type
#'
#' @param net_type Network type. One of 'signed', 'signed hybrid' or 'unsigned'.
#'
#' @return Character of TOM type
#' @noRd
get_TOMtype <- function(net_type) {
  if(net_type == "signed hybrid") {
    TOMType <- "signed"
  } else if(net_type == "signed") {
    TOMType <- "signed Nowick"
  } else {
    TOMType <- "unsigned"
  }
  return(TOMType)
}


#' Wrapper to handle variable type for trait object
#'
#' @param metadata A data frame containing sample names in row names and
#' sample annotation in the first column.
#' @param continuous_trait Logical indicating if trait is a continuous variable.
#' Default is FALSE.
#'
#' @return Processed trait object.
#' @noRd
handle_trait_type <- function(metadata, continuous_trait = FALSE) {
  if(!continuous_trait) {
    sampleinfo <- cbind(Samples=rownames(metadata), metadata)
    tmpdir <- tempdir()
    tmpfile <- tempfile(tmpdir = tmpdir, fileext = "traitmatrix.txt")
    tablesamples <- table(sampleinfo)
    write.table(tablesamples, file = tmpfile,
                quote = FALSE, sep="\t", row.names=TRUE)
    trait <- read.csv(tmpfile, header=TRUE,
                      sep="\t", row.names=1, stringsAsFactors = FALSE)
    unlink(tmpfile)
  } else {
    trait <- metadata
  }
  return(trait)
}


#' Transform a correlation matrix to an edge list
#'
#' @param matrix Symmetrical correlation matrix.
#'
#' @return A 2-column data frame containing node 1, node 2 and edge weight.
#' @export
#' @rdname cormat_to_edgelist
#' @examples
#' data(filt.se)
#' cor_mat <- cor(t(SummarizedExperiment::assay(filt.se)))
#' edgelist <- cormat_to_edgelist(cor_mat)
cormat_to_edgelist <- function(matrix) {
  edgelist <- matrix
  edgelist[lower.tri(edgelist, diag=TRUE)] <- NA
  edgelist <- na.omit(data.frame(as.table(edgelist), stringsAsFactors=FALSE))
  colnames(edgelist) <- c("Node1", "Node2", "Weight")
  edgelist$Node1 <- as.character(edgelist$Node1)
  edgelist$Node2 <- as.character(edgelist$Node2)
  edgelist$Weight <- as.numeric(edgelist$Weight)
  return(edgelist)
}


#' Check scale-free topology fit for a given network
#'
#' @param edgelist Edge list as a data frame containing node 1,
#' node 2 and edge weight.
#' @param net_type Type of biological network. One of "gcn", "grn", or "ppi".
#' Default: gcn.
#'
#' @return A list with SFT fit statistics and a message indicating if
#' the network is scale-free.
#' @rdname check_SFT
#' @export
#' @importFrom igraph graph_from_data_frame as_adjacency_matrix fit_power_law
#' @examples
#' set.seed(1)
#' exp <- t(matrix(rnorm(10000), ncol=1000, nrow=200))
#' rownames(exp) <- paste0("Gene", 1:nrow(exp))
#' colnames(exp) <- paste0("Sample", 1:ncol(exp))
#' cormat <- cor(t(exp))
#' edges <- cormat_to_edgelist(cormat)
#' edges <- edges[abs(edges$Weight) > 0.10, ]
#' check_SFT(edges)
check_SFT <- function(edgelist, net_type = "gcn") {
  
  # Calculate degree of the resulting graph
  if(net_type == "gcn") {
    graph <- igraph::graph_from_data_frame(edgelist, directed=FALSE)
    adj <- igraph::as_adjacency_matrix(graph, sparse = FALSE)
    diag(adj) <- 0
    degree <- apply(adj, 1, sum, na.rm=TRUE)
  } else if(net_type == "grn") {
    graph <- igraph::graph_from_data_frame(edgelist, directed=TRUE)
    degree <- igraph::degree(graph, mode = "out")
  } else if(net_type == "ppi") {
    graph <- igraph::graph_from_data_frame(edgelist, directed=FALSE)
    degree <- igraph::degree(graph, mode)
  } else {
    stop("Invalid network type. Please, input one of 'gcn', 'grn', or 'ppi'.")
  }
  
  # Test for scale-free topology fit
  test <- igraph::fit_power_law(degree)
  if(test$KS.p < 0.05) {
    message("At the 95% confidence level for the Kolmogorov-Smirnov statistic, your graph does not fit the scale-free topology. P-value:", test$KS.p)
  } else {
    message("Your graph fits the scale-free topology. P-value:", test$KS.p)
  }
  
  return(test)
}


#' Helper to handle list of SummarizedExperiment objects as input
#'
#' @param exp List of data frames or SummarizedExperiment objects.
#'
#' @return If exp is a list of SummarizedExperiment objects,
#' it will return a list of data frames. Otherwise, it will simply return exp as it is.
#' @noRd
#' @importFrom SummarizedExperiment assay
handleSElist <- function(exp) {
  if(is(exp[[1]], "SummarizedExperiment")) {
    list <- lapply(exp, handleSE)
  } else {
    list <- exp
  }
  return(list)
}


#' Helper function to handle metadata input for consensus modules identification
#'
#' @param exp List of data frames or SummarizedExperiment objects.
#' @param metadata A data frame containing sample names in row names and sample annotation in the first column.
#' @return Data frame of metadata for all expression sets.
#' @noRd
#' @importFrom SummarizedExperiment colData
handle_metadata <- function(exp, metadata) {
  if(is(exp[[1]], "SummarizedExperiment")) {
    metadata <- Reduce(rbind, lapply(exp, function(x) {
      meta <- as.data.frame(SummarizedExperiment::colData(x))
      return(meta)
    }))
  } else {
    metadata <- metadata
  }
  return(metadata)
}



#' Convert p-values in matrix to symbols
#'
#' @param matrix Matrix of p-values.
#'
#' @return Matrix of symbols.
#' @noRd
pval2symbol <- function(matrix) {
  modtraitsymbol <- matrix
  modtraitsymbol[modtraitsymbol < 0.001] <- "***"
  modtraitsymbol[modtraitsymbol >= 0.001 & modtraitsymbol < 0.01] <- "**"
  modtraitsymbol[modtraitsymbol >= 0.01 & modtraitsymbol < 0.05] <- "*"
  modtraitsymbol[modtraitsymbol >= 0.05] <- ""
  return(modtraitsymbol)
}





######### modified exp2gcn function to tune minimum module size
exp2gcn_1 = function (exp, net_type = "signed", module_merging_threshold = 0.8,include_TOM = F, 
                      SFTpower = NULL, cor_method = "spearman", verbose = FALSE, min_module_size = 100) 
{
  params <- list(net_type = net_type, module_merging_threshold = module_merging_threshold, 
                 SFTpower = SFTpower, cor_method = cor_method)
  norm.exp <- handleSE(exp)
  if (is.null(SFTpower)) {
    stop("Please, specify the SFT power.")
  }
  if (verbose) {
    message("Calculating adjacency matrix...")
  }
  cor_matrix <- calculate_cor_adj(cor_method, norm.exp, SFTpower, 
                                  net_type)$cor
  adj_matrix <- calculate_cor_adj(cor_method, norm.exp, SFTpower, 
                                  net_type)$adj
  gene_ids <- rownames(adj_matrix)
  adj_matrix <- matrix(adj_matrix, nrow = nrow(adj_matrix))
  rownames(adj_matrix) <- gene_ids
  colnames(adj_matrix) <- gene_ids
  if (verbose) {
    message("Calculating topological overlap matrix (TOM)...")
  }
  tomtype <- get_TOMtype(net_type)
  TOM <- WGCNA::TOMsimilarity(adj_matrix, TOMType = tomtype)
  dissTOM <- 1 - TOM
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  if (verbose) {
    message("Detecting coexpression modules...")
  }
  old.module_labels <- dynamicTreeCut::cutreeDynamicTree(dendro = geneTree, 
                                                         minModuleSize = min_module_size, deepSplit = TRUE)
  nmod <- length(unique(old.module_labels))
  palette <- rev(WGCNA::standardColors(nmod))
  old.module_colors <- WGCNA::labels2colors(old.module_labels, 
                                            colorSeq = palette)
  if (verbose) {
    message("Calculating module eigengenes (MEs)...")
  }
  old.MElist <- WGCNA::moduleEigengenes(t(norm.exp), colors = old.module_colors, 
                                        softPower = SFTpower)
  old.MEs <- old.MElist$eigengenes
  MEDiss1 <- 1 - cor(old.MEs)
  old.METree <- hclust(as.dist(MEDiss1), method = "average")
  MEDissThreshold <- 1 - module_merging_threshold
  if (verbose) {
    message("Merging similar modules...")
  }
  if (cor_method == "pearson") {
    merge1 <- WGCNA::mergeCloseModules(t(norm.exp), old.module_colors, 
                                       cutHeight = MEDissThreshold, verbose = 0, colorSeq = palette)
  }
  else if (cor_method == "spearman") {
    merge1 <- WGCNA::mergeCloseModules(t(norm.exp), old.module_colors, 
                                       cutHeight = MEDissThreshold, verbose = 0, corOptions = list(use = "p", 
                                                                                                   method = "spearman"), colorSeq = palette)
  }
  else if (cor_method == "biweight") {
    merge1 <- WGCNA::mergeCloseModules(t(norm.exp), old.module_colors, 
                                       cutHeight = MEDissThreshold, verbose = 0, corFnc = bicor, 
                                       colorSeq = palette)
  }
  else {
    stop("Please, specify a correlation method. One of 'spearman', 'pearson' or 'biweight'.")
  }
  new.module_colors <- merge1$colors
  new.MEs <- merge1$newMEs
  new.METree <- hclust(as.dist(1 - cor(new.MEs)), method = "average")
  genes_and_modules <- as.data.frame(cbind(gene_ids, new.module_colors), 
                                     stringsAsFactors = FALSE)
  colnames(genes_and_modules) <- c("Genes", "Modules")
  if (verbose) {
    message("Calculating intramodular connectivity...")
  }
  kwithin <- WGCNA::intramodularConnectivity(adj_matrix, new.module_colors)
  if (include_TOM == T){
  result.list <- list(tom_matrix = TOM,adjacency_matrix = adj_matrix, MEs = new.MEs, 
                      genes_and_modules = genes_and_modules, kIN = kwithin, 
                      moduleColors = new.module_colors, correlation_matrix = cor_matrix, 
                      params = params, dendro_plot_objects = list(tree = geneTree, 
                                                                  unmerged = old.module_colors))
  }else{
    result.list <- list(adjacency_matrix = adj_matrix, MEs = new.MEs, 
                        genes_and_modules = genes_and_modules, kIN = kwithin, 
                        moduleColors = new.module_colors, correlation_matrix = cor_matrix, 
                        params = params, dendro_plot_objects = list(tree = geneTree, 
                                                                    unmerged = old.module_colors))
  }
  return(result.list)
}


consensus_modules_1=function (exp_list, metadata, power, cor_method = "spearman", 
                              net_type = "signed hybrid", module_merging_threshold = 0.8, min_cluster_size = 100,
                              verbose = FALSE) 
{
  metadata <- handle_metadata(exp_list, metadata)
  exp_list <- handleSElist(exp_list)
  nSets <- length(exp_list)
  multiExp <- lapply(exp_list, function(x) {
    element <- list(data = as.data.frame(t(x)))
    return(element)
  })
  expSize <- WGCNA::checkSets(multiExp)
  nGenes <- expSize$nGenes
  nSamples <- expSize$nSamples
  sampleinfo <- lapply(seq_along(multiExp), function(x) {
    sinfo <- metadata[rownames(multiExp[[x]]$data), , drop = FALSE]
    return(sinfo)
  })
  if (verbose) {
    message("Calculating adjacency matrix...")
  }
  adj <- lapply(seq_len(nSets), function(x) {
    if (cor_method == "pearson") {
      adjacencies <- WGCNA::adjacency(multiExp[[x]]$data, 
                                      power = power[x], type = net_type)
    }
    else if (cor_method == "spearman") {
      adjacencies <- WGCNA::adjacency(multiExp[[x]]$data, 
                                      power = power[x], type = net_type, corOptions = list(use = "p", 
                                                                                           method = "spearman"))
    }
    else if (cor_method == "biweight") {
      adjacencies <- WGCNA::adjacency(multiExp[[x]]$data, 
                                      power = power[x], type = net_type, corFnc = bicor)
    }
    else {
      stop("Please, specify a correlation method. One of 'spearman', 'pearson' or 'biweight'.")
    }
    return(adjacencies)
  })
  if (verbose) {
    message("Calculating topological overlap matrix (TOM)...")
  }
  TOM <- lapply(adj, function(x) {
    if (net_type == "signed hybrid") {
      tom <- WGCNA::TOMsimilarity(x, TOMType = "signed")
    }
    else if (net_type == "signed") {
      tom <- WGCNA::TOMsimilarity(x, TOMType = "signed Nowick")
    }
    else {
      tom <- WGCNA::TOMsimilarity(x, TOMType = "unsigned")
    }
    return(tom)
  })
  scaleP <- 0.95
  nSamples <- as.integer(1/(1 - scaleP) * 1000)
  scaleSample <- sample(nGenes * (nGenes - 1)/2, size = nSamples)
  scaledTOM <- lapply(seq_len(nSets), function(x) {
    TOMScalingSamples <- as.dist(TOM[[x]])[scaleSample]
    scaleQuant <- quantile(TOMScalingSamples, probs = scaleP, 
                           type = 8)
    tom_scaled <- TOM[[x]]
    if (x > 1) {
      scalePowers <- log(scaleQuant[1])/log(scaleQuant)
      tom_scaled <- tom_scaled^scalePowers
    }
    return(tom_scaled)
  })
  if (nSets <= 3) {
    consensusTOM <- do.call(pmin, lapply(seq_along(TOM), 
                                         function(i) TOM[[i]]))
  }
  else {
    consensusTOM <- do.call(function(x) WGCNA::pquantile(x, 
                                                         prob = 0.25), lapply(seq_along(TOM), function(i) TOM[[i]]))
  }
  consTree <- hclust(as.dist(1 - consensusTOM), method = "average")
  unmergedLabels <- dynamicTreeCut::cutreeDynamic(dendro = consTree, 
                                                  distM = 1 - consensusTOM, deepSplit = 2, cutHeight = 0.995, 
                                                  minClusterSize = min_cluster_size, pamRespectsDendro = FALSE)
  nmod <- length(unique(unmergedLabels))
  palette <- rev(WGCNA::standardColors(nmod))
  unmergedColors <- WGCNA::labels2colors(unmergedLabels, colorSeq = palette)
  unmergedMEs <- WGCNA::multiSetMEs(multiExp, colors = NULL, 
                                    universalColors = unmergedColors)
  consMEDiss <- WGCNA::consensusMEDissimilarity(unmergedMEs)
  consMETree <- hclust(as.dist(consMEDiss), method = "average")
  merging_threshold <- 1 - module_merging_threshold
  merge <- WGCNA::mergeCloseModules(multiExp, unmergedColors, 
                                    cutHeight = merging_threshold, verbose = 0)
  moduleLabels <- merge$colors
  moduleColors <- WGCNA::labels2colors(moduleLabels, colorSeq = palette)
  genes_cmod <- data.frame(Genes = rownames(adj[[1]]), Cons_modules = moduleColors)
  consMEs <- merge$newMEs
  result_list <- list(consModules = moduleColors, consMEs = consMEs, 
                      exprSize = expSize, sampleInfo = sampleinfo, genes_cmodules = genes_cmod, 
                      dendro_plot_objects = list(tree = consTree, unmerged = unmergedColors))
  return(result_list)
}


SFT_fit_1 = function (exp, net_type = "signed", rsquared = 0.8, cor_method = "spearman") 
{
  exp <- handleSE(exp)
  texp <- t(exp)
  if (cor_method == "pearson") {
    sft <- WGCNA::pickSoftThreshold(texp, networkType = net_type, 
                                    powerVector = 3:20, RsquaredCut = rsquared)
  }
  else if (cor_method == "biweight") {
    sft <- WGCNA::pickSoftThreshold(texp, networkType = net_type, 
                                    powerVector = 3:20, RsquaredCut = rsquared, corFnc = bicor, 
                                    corOptions = list(use = "p", maxPOutliers = 0.05))
  }
  else if (cor_method == "spearman") {
    sft <- WGCNA::pickSoftThreshold(texp, networkType = net_type, 
                                    powerVector = 3:20, RsquaredCut = rsquared, corOptions = list(use = "p", 
                                                                                                  method = "spearman"))
  }
  else {
    stop("Please, specify a correlation method (one of 'spearman', 'pearson' or 'biweight').")
  }
  wgcna_power <- sft$powerEstimate
  if (is.na(wgcna_power)) {
    message("No power reached R-squared cut-off, now choosing max R-squared based power")
    wgcna_power <- sft$fitIndices$Power[which(sft$fitIndices$SFT.R.sq == 
                                                max(sft$fitIndices$SFT.R.sq))]
  }
  sft_df <- data.frame(power = sft$fitIndices$Power, fit = -sign(sft$fitIndices$slope) * 
                         sft$fitIndices$SFT.R.sq, meank = sft$fitIndices$mean.k.)
  p1 <- ggpubr::ggscatter(sft_df, x = "power", y = "fit", 
                           ylim = c(0, 1), 
                          label = "power", size = 1, font.label = 10, color = "gray10", 
                          font.tickslab = 10, ytickslab.rt = 90) + geom_hline(yintercept = rsquared, 
                                                                              color = "brown3") +
    theme(plot.title = element_text(hjust = 0.5)) + xlab("Soft threshold (power)") +
    ylab(expression(paste("Scale-free topology fit - ", 
                          R^{2}) ) ) + ggtitle("Scale independence")
    
  p2 <- ggpubr::ggscatter(sft_df, x = "power", y = "meank", xlab = "Soft threshold (power)", 
                          ylab = "Mean connectivity (k)", title = "Mean connectivity", 
                          label = "power", size = 1, font.label = 10, color = "gray10", 
                          font.tickslab = 10, ytickslab.rt = 90) + 
    theme(plot.title = element_text(hjust = 0.5))
  sft_plot <- ggpubr::ggarrange(p1, p2)
  result <- list(power = as.numeric(wgcna_power), plot = sft_plot, dtf = sft_df, plot1=p1, plot2=p2)
  return(result)
}



##### SystemPiper vennPlot modified function
vennplot_1 = function (x, mymain = "Venn Diagram", mysub = "default", setlabels = "default", 
  yoffset = seq(0, 10, by = 0.34), ccol = rep(1, 31), colmode = 1, 
  lcol = c("#FF0000", "#008B00", "#0000FF", "#FF00FF", "#CD8500"), 
  lines = c("#FF0000", "#008B00", "#0000FF", "#FF00FF", "#CD8500"), 
  mylwd = 2, diacol = 1, type = "ellipse", ccex = 1, lcex = 1, 
  sepsplit = "_", ...) 
{
  if (!any(c(class(x) == "VENNset", class(x) == "list", is.numeric(x)))) {
    stop("x needs to be one of: VENNset, list of VENNsets, named numeric vector, or list of named numeric vectors.")
  }
  if (class(x) == "list") {
    if (length(unique(sapply(x, length))) != 1) 
      stop("List components need to have identical length.")
  }
  if (class(x) == "VENNset") {
    counts <- list(sapply(vennlist(x), length))
    myclass <- "VENNset"
  }
  else if (class(x) == "list" & all(sapply(x, class) == "VENNset")) {
    counts <- lapply(x, function(y) sapply(vennlist(y), 
      length))
    myclass <- "VENNset"
  }
  else if (is.numeric(x) & is.list(x) == FALSE) {
    counts <- list(x)
    myclass <- "numeric"
  }
  else if (class(x) == "list" & all(sapply(x, is.numeric))) {
    counts <- x
    myclass <- "numeric"
  }
  else {
    stop("x needs to be one of: VENNset, list of VENNsets, named numeric vector, or list of named numeric vectors.")
  }
  if (!length(counts[[1]]) %in% c(3, 7, 15, 31)) 
    stop("Only 2-5 way venn comparisons are supported.")
  grepLabel <- function(label, x = names(counts[[1]])) {
    x <- strsplit(x, sepsplit)
    as.numeric(which(sapply(x, function(y) any(y == label))))
  }
  if (length(counts[[1]]) == 3) {
    if (mysub == "default") {
      if (myclass == "numeric") {
        n <- names(counts[[1]])[1:2]
        if (!all(rowSums(sapply(n, function(x) sapply(n, 
          function(y) grepl(y, x)))) == 1)) {
          sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, 
            names(counts[[1]]))]))
          if (sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, 
            names(counts[[1]][-c(1:length(n))])))) {
            sample_counts <- rep("?", length(n))
            warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by the character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")
          }
        }
        else {
          sample_counts <- sapply(n, function(x) sum(counts[[1]][grep(x, 
            names(counts[[1]]))]))
        }
        mysub <- paste(paste("Unique objects: All =", 
          sum(counts[[1]])), paste("; S1 =", sample_counts[1]), 
          paste("; S2 =", sample_counts[2]), sep = "")
      }
      else if (myclass == "VENNset") {
        if (class(x) == "list") 
          x <- x[[1]]
        sample_counts <- sapply(setlist(x), function(y) unique(length(y)))
        mysub <- paste(paste("Unique objects: All =", 
          length(unique(unlist(setlist(x))))), paste("; S1 =", 
          sample_counts[1]), paste("; S2 =", sample_counts[2]), 
          sep = "")
      }
      else {
        mysub <- mysub
      }
    }
    graphics::symbols(x = c(4, 6), y = c(4, 4), circles = c(2, 
      2), xlim = c(0, 11), ylim = c(0, 11), inches = FALSE, 
      main = mymain, sub = mysub, lwd = mylwd, xlab = "", 
      ylab = "",xaxt = "n", 
            yaxt = "n", bty = "n",  fg = lines, 
      ...)
    for (i in seq(along = counts)) {
      olDF <- data.frame(x = c(3.1, 7, 5), y = c(4, 4, 
        4), counts = counts[[i]])
      if (colmode == 1) {
        graphics::text(olDF$x, olDF$y + yoffset[i], 
          olDF$counts, col = ccol, cex = ccex, ...)
      }
      if (colmode == 2) {
        graphics::text(olDF$x, olDF$y + yoffset[i], 
          olDF$counts, col = ccol[[i]], cex = ccex[i], 
          ...)
      }
    }
    if (length(setlabels) == 1 & setlabels[1] == "default") {
      setlabels <- names(counts[[1]][1:2])
    }
    else {
      setlabels <- setlabels
    }
    graphics::text(c(2, 8), c(8.8, 8.8), labels = setlabels, 
      col = lcol, cex = lcex, ...)
  }
  if (length(counts[[1]]) == 7) {
    if (mysub == "default") {
      if (myclass == "numeric") {
        n <- names(counts[[1]])[1:3]
        if (!all(rowSums(sapply(n, function(x) sapply(n, 
          function(y) grepl(y, x)))) == 1)) {
          sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, 
            names(counts[[1]]))]))
          if (sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, 
            names(counts[[1]][-c(1:length(n))])))) {
            sample_counts <- rep("?", length(n))
            warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by the character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")
          }
        }
        else {
          sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, 
            names(counts[[1]]))]))
        }
        mysub <- paste(paste("Unique objects: All =", 
          sum(counts[[1]])), paste("; S1 =", sample_counts[1]), 
          paste("; S2 =", sample_counts[2]), paste("; S3 =", 
            sample_counts[3]), sep = "")
      }
      else if (myclass == "VENNset") {
        if (class(x) == "list") 
          x <- x[[1]]
        sample_counts <- sapply(setlist(x), function(y) unique(length(y)))
        mysub <- paste(paste("Unique objects: All =", 
          length(unique(unlist(setlist(x))))), paste("; S1 =", 
          sample_counts[1]), paste("; S2 =", sample_counts[2]), 
          paste("; S3 =", sample_counts[3]), sep = "")
      }
      else {
        mysub <- mysub
      }
    }
    graphics::symbols(x = c(4, 6, 5), y = c(6, 6, 4), circles = c(2, 
      2, 2), xlim = c(0, 10), ylim = c(0, 10), inches = FALSE, 
      main = mymain, sub = mysub, lwd = mylwd, xlab = "", 
      ylab = "", xaxt = "n", yaxt = "n", bty = "n", fg = lines, 
      ...)
    for (i in seq(along = counts)) {
      olDF <- data.frame(x = c(3, 7, 5, 5, 3.8, 6.3, 5), 
        y = c(6.5, 6.5, 3, 7, 4.6, 4.6, 5.3), counts = counts[[i]])
      if (colmode == 1) {
        graphics::text(olDF$x, olDF$y + yoffset[i], 
          olDF$counts, col = ccol, cex = ccex, ...)
      }
      if (colmode == 2) {
        graphics::text(olDF$x, olDF$y + yoffset[i], 
          olDF$counts, col = ccol[[i]], cex = ccex[i], 
          ...)
      }
    }
    if (length(setlabels) == 1 & setlabels[1] == "default") {
      setlabels <- names(counts[[1]][1:3])
    }
    else {
      setlabels <- setlabels
    }
    graphics::text(c(2, 8, 6), c(8.8, 8.8, 1.1), labels = setlabels, 
      col = lcol, cex = lcex, ...)
  }
  if (length(counts[[1]]) == 15 & type == "ellipse") {
    if (mysub == "default") {
      if (myclass == "numeric") {
        n <- names(counts[[1]])[1:4]
        if (!all(rowSums(sapply(n, function(x) sapply(n, 
          function(y) grepl(y, x)))) == 1)) {
          sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, 
            names(counts[[1]]))]))
          if (sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, 
            names(counts[[1]][-c(1:length(n))])))) {
            sample_counts <- rep("?", length(n))
            warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by the character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")
          }
        }
        else {
          sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, 
            names(counts[[1]]))]))
        }
        mysub <- paste(paste("Unique objects: All =", 
          sum(counts[[1]])), paste("; S1 =", sample_counts[1]), 
          paste("; S2 =", sample_counts[2]), paste("; S3 =", 
            sample_counts[3]), paste("; S4 =", sample_counts[4]), 
          sep = "")
      }
      else if (myclass == "VENNset") {
        if (class(x) == "list") 
          x <- x[[1]]
        sample_counts <- sapply(setlist(x), function(y) unique(length(y)))
        mysub <- paste(paste("Unique objects: All =", 
          length(unique(unlist(setlist(x))))), paste("; S1 =", 
          sample_counts[1]), paste("; S2 =", sample_counts[2]), 
          paste("; S3 =", sample_counts[3]), paste("; S4 =", 
            sample_counts[4]), sep = "")
      }
      else {
        mysub <- mysub
      }
    }
    plotellipse <- function(center = c(1, 1), radius = c(1, 
      2), rotate = 1, segments = 360, xlab = "", ylab = "", 
      ...) {
      angles <- (0:segments) * 2 * pi/segments
      rotate <- rotate * pi/180
      ellipse <- cbind(radius[1] * cos(angles), radius[2] * 
        sin(angles))
      ellipse <- cbind(ellipse[, 1] * cos(rotate) + ellipse[, 
        2] * sin(rotate), ellipse[, 2] * cos(rotate) - 
        ellipse[, 1] * sin(rotate))
      ellipse <- cbind(center[1] + ellipse[, 1], center[2] + 
        ellipse[, 2])
      plot(ellipse, type = "l", xlim = c(0, 10), ylim = c(0, 
        10), xlab = "", ylab = "", ...)
    }
    ellipseVenn <- function(...) {
      graphics::split.screen(c(1, 1))
      plotellipse(center = c(3.5, 3.6), radius = c(2, 
        4), rotate = -35, segments = 360, xlab = "", 
        ylab = "", col = lines[1], axes = FALSE, main = mymain, 
        sub = mysub, lwd = mylwd, ...)
      graphics::screen(1, new = FALSE)
      plotellipse(center = c(4.7, 4.4), radius = c(2, 
        4), rotate = -35, segments = 360, xlab = "", 
        ylab = "", col = lines[2], axes = FALSE, lwd = mylwd, 
        ...)
      graphics::screen(1, new = FALSE)
      plotellipse(center = c(5.3, 4.4), radius = c(2, 
        4), rotate = 35, segments = 360, xlab = "", 
        ylab = "", col = lines[3], axes = FALSE, lwd = mylwd, 
        ...)
      graphics::screen(1, new = FALSE)
      plotellipse(center = c(6.5, 3.6), radius = c(2, 
        4), rotate = 35, segments = 360, xlab = "", 
        ylab = "", col = lines[4], axes = FALSE, lwd = mylwd, 
        ...)
      for (i in seq(along = counts)) {
        olDF <- data.frame(x = c(1.5, 3.5, 6.5, 8.5, 
          2.9, 3.1, 5, 5, 6.9, 7.1, 3.6, 5.8, 4.2, 6.4, 
          5), y = c(4.8, 7.2, 7.2, 4.8, 5.9, 2.2, 0.7, 
          6, 2.2, 5.9, 4, 1.4, 1.4, 4, 2.8), counts = counts[[i]])
        if (colmode == 1) {
          graphics::text(olDF$x, olDF$y + yoffset[i], 
            olDF$counts, col = ccol, cex = ccex, ...)
        }
        if (colmode == 2) {
          graphics::text(olDF$x, olDF$y + yoffset[i], 
            olDF$counts, col = ccol[[i]], cex = ccex[i], 
            ...)
        }
      }
      if (length(setlabels) == 1 & setlabels[1] == "default") {
        setlabels <- names(counts[[1]][1:4])
      }
      else {
        setlabels <- setlabels
      }
      graphics::text(c(0.4, 2.8, 7.5, 9.4), c(7.3, 8.3, 
        8.3, 7.3), labels = setlabels, col = lcol, cex = lcex, 
        ...)
      graphics::close.screen(all.screens = TRUE)
    }
    ellipseVenn(...)
  }
  if (length(counts[[1]]) == 15 & type == "circle") {
    if (mysub == "default") {
      if (myclass == "numeric") {
        n <- names(counts[[1]])[1:4]
        if (!all(rowSums(sapply(n, function(x) sapply(n, 
          function(y) grepl(y, x)))) == 1)) {
          sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, 
            names(counts[[1]]))]))
          if (sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, 
            names(counts[[1]][-c(1:length(n))])))) {
            sample_counts <- rep("?", length(n))
            warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by the character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")
          }
        }
        else {
          sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, 
            names(counts[[1]]))]))
        }
        mysub <- paste(paste("Unique objects: All =", 
          sum(counts[[1]])), paste("; S1 =", sample_counts[1]), 
          paste("; S2 =", sample_counts[2]), paste("; S3 =", 
            sample_counts[3]), paste("; S4 =", sample_counts[4]), 
          sep = "")
      }
      else if (myclass == "VENNset") {
        if (class(x) == "list") 
          x <- x[[1]]
        sample_counts <- sapply(setlist(x), function(y) unique(length(y)))
        mysub <- paste(paste("Unique objects: All =", 
          length(unique(unlist(setlist(x))))), paste("; S1 =", 
          sample_counts[1]), paste("; S2 =", sample_counts[2]), 
          paste("; S3 =", sample_counts[3]), paste("; S4 =", 
            sample_counts[4]), sep = "")
      }
      else {
        mysub <- mysub
      }
    }
    graphics::symbols(x = c(4, 5.5, 4, 5.5), y = c(6, 6, 
      4.5, 4.5), circles = c(2, 2, 2, 2), xlim = c(0, 
      10), ylim = c(0, 10), inches = FALSE, main = mymain, 
      sub = mysub, lwd = mylwd, xlab = "", ylab = "", 
      xaxt = "n", yaxt = "n", bty = "n", fg = lines, ...)
    for (i in seq(along = counts)) {
      olDF <- data.frame(x = c(3, 6.5, 3, 6.5, 4.8, 3, 
        4.8, 4.8, 6.5, 4.8, 3.9, 5.7, 3.9, 5.7, 4.8), 
        y = c(7.2, 7.2, 3.2, 3.2, 7.2, 5.2, 0.4, 0.4, 
          5.2, 3.2, 6.3, 6.3, 4.2, 4.2, 5.2), counts = counts[[i]])
      if (colmode == 1) {
        graphics::text(olDF$x[-c(7, 8)], olDF$y[-c(7, 
          8)] + yoffset[i], olDF$counts[-c(7, 8)], col = ccol, 
          cex = ccex, ...)
      }
      if (colmode == 2) {
        graphics::text(olDF$x[-c(7, 8)], olDF$y[-c(7, 
          8)] + yoffset[i], olDF$counts[-c(7, 8)], col = ccol[[i]], 
          cex = ccex[i], ...)
      }
      graphics::text(c(4.8), c(0.8) + yoffset[i], paste("Only in ", 
        names(counts[[1]][1]), " & ", names(counts[[1]][4]), 
        ": ", olDF$counts[7], "; Only in ", names(counts[[1]][2]), 
        " & ", names(counts[[1]][3]), ": ", olDF$counts[8], 
        sep = ""), col = diacol, cex = ccex, ...)
    }
    if (length(setlabels) == 1 & setlabels[1] == "default") {
      setlabels <- names(counts[[1]][1:4])
    }
    else {
      setlabels <- setlabels
    }
    graphics::text(c(2, 7.5, 2, 7.5), c(8.3, 8.3, 2, 2), 
      labels = setlabels, col = lcol, cex = lcex, ...)
  }
  if (length(counts[[1]]) == 31) {
    if (mysub == "default") {
      if (myclass == "numeric") {
        n <- names(counts[[1]])[1:5]
        if (!all(rowSums(sapply(n, function(x) sapply(n, 
          function(y) grepl(y, x)))) == 1)) {
          sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, 
            names(counts[[1]]))]))
          if (sum(grepl(sepsplit, n)) > 0 | !all(grepl(sepsplit, 
            names(counts[[1]][-c(1:length(n))])))) {
            sample_counts <- rep("?", length(n))
            warning("Set labels are substrings of one another. To fix this, the set labels need to be separated by the character provided under \"sepsplit\", but the individual names cannot contain this character themselves.")
          }
        }
        else {
          sample_counts <- sapply(n, function(x) sum(counts[[1]][grepLabel(x, 
            names(counts[[1]]))]))
        }
        mysub <- paste(paste("Unique objects: All =", 
          sum(counts[[1]])), paste("; S1 =", sample_counts[1]), 
          paste("; S2 =", sample_counts[2]), paste("; S3 =", 
            sample_counts[3]), paste("; S4 =", sample_counts[4]), 
          paste("; S5 =", sample_counts[5]), sep = "")
      }
      else if (myclass == "VENNset") {
        if (class(x) == "list") 
          x <- x[[1]]
        sample_counts <- sapply(setlist(x), function(y) unique(length(y)))
        mysub <- paste(paste("Unique objects: All =", 
          length(unique(unlist(setlist(x))))), paste("; S1 =", 
          sample_counts[1]), paste("; S2 =", sample_counts[2]), 
          paste("; S3 =", sample_counts[3]), paste("; S4 =", 
            sample_counts[4]), paste("; S5 =", sample_counts[5]), 
          sep = "")
      }
      else {
        mysub <- mysub
      }
    }
    plotellipse <- function(center = c(1, 1), radius = c(1, 
      2), rotate = 1, segments = 360, xlab = "", ylab = "", 
      ...) {
      angles <- (0:segments) * 2 * pi/segments
      rotate <- rotate * pi/180
      ellipse <- cbind(radius[1] * cos(angles), radius[2] * 
        sin(angles))
      ellipse <- cbind(ellipse[, 1] * cos(rotate) + ellipse[, 
        2] * sin(rotate), ellipse[, 2] * cos(rotate) - 
        ellipse[, 1] * sin(rotate))
      ellipse <- cbind(center[1] + ellipse[, 1], center[2] + 
        ellipse[, 2])
      plot(ellipse, type = "l", xlim = c(0, 10), ylim = c(0, 
        10), xlab = "", ylab = "", ...)
    }
    ellipseVenn <- function(...) {
      graphics::split.screen(c(1, 1))
      graphics::screen(1, new = FALSE)
      plotellipse(center = c(4.83, 6.2), radius = c(1.43, 
        4.11), rotate = 0, segments = 360, xlab = "", 
        ylab = "", col = lines[1], axes = FALSE, main = mymain, 
        sub = mysub, lwd = mylwd, ...)
      graphics::screen(1, new = FALSE)
      plotellipse(center = c(6.25, 5.4), radius = c(1.7, 
        3.6), rotate = 66, segments = 360, xlab = "", 
        ylab = "", col = lines[2], axes = FALSE, lwd = mylwd, 
        ...)
      graphics::screen(1, new = FALSE)
      plotellipse(center = c(6.1, 3.5), radius = c(1.55, 
        3.9), rotate = 150, segments = 360, xlab = "", 
        ylab = "", col = lines[3], axes = FALSE, lwd = mylwd, 
        ...)
      graphics::screen(1, new = FALSE)
      plotellipse(center = c(4.48, 3.15), radius = c(1.55, 
        3.92), rotate = 210, segments = 360, xlab = "", 
        ylab = "", col = lines[4], axes = FALSE, lwd = mylwd, 
        ...)
      graphics::screen(1, new = FALSE)
      plotellipse(center = c(3.7, 4.8), radius = c(1.7, 
        3.6), rotate = 293.5, segments = 360, xlab = "", 
        ylab = "", col = lines[5], axes = FALSE, lwd = mylwd, 
        ...)
      for (i in seq(along = counts)) {
        olDF <- data.frame(x = c(4.85, 8, 7.1, 3.5, 
          2, 5.9, 4.4, 4.6, 3.6, 7.1, 6.5, 3.2, 5.4, 
          6.65, 3.4, 5, 6.02, 3.6, 5.2, 4.03, 4.2, 6.45, 
          6.8, 3.39, 6.03, 5.74, 4.15, 3.95, 5.2, 6.4, 
          5.1), y = c(8.3, 6.2, 1.9, 1.6, 5.4, 6.85, 
          6.6, 2.45, 6.4, 4.3, 6, 4.6, 2.1, 3.4, 3.25, 
          6.43, 6.38, 5.1, 2.49, 6.25, 3.08, 5.3, 4, 
          3.8, 3.2, 5.95, 5.75, 3.75, 3, 4.5, 4.6), 
          counts = counts[[i]])
        if (colmode == 1) {
          graphics::text(olDF$x, olDF$y + yoffset[i], 
            olDF$counts, col = ccol, cex = ccex, ...)
        }
        if (colmode == 2) {
          graphics::text(olDF$x, olDF$y + yoffset[i], 
            olDF$counts, col = ccol[[i]], cex = ccex[i], 
            ...)
        }
      }
      if (length(setlabels) == 1 & setlabels[1] == "default") {
        setlabels <- names(counts[[1]][1:5])
      }
      else {
        setlabels <- setlabels
      }
      graphics::text(c(5.7, 7.9, 8.5, 4.2, 0.8), c(9.9, 
        7.9, 1.9, 0, 7.3), adj = c(0, 0.5), labels = setlabels, 
        col = lcol, cex = lcex, ...)
      graphics::close.screen(all.screens = TRUE)
    }
    ellipseVenn(...)
  }
}



m_heatmap = function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, lwid = c(0.05,4.00), lhei = c(0.25, 4.00),
                      distfun = dist, hclustfun = hclust, reorderfun = function(d, 
                                                                                w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv, 
                                                                                                                                           "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE, 
                      margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 + 
                        1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
                      labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
                      verbose = getOption("verbose"), ...) 
{
  scale <- if (symm && missing(scale)) 
    "none"
  else match.arg(scale)
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
    stop("'x' must be a numeric matrix")
  nr <- di[1L]
  nc <- di[2L]
  if (nr <= 1 || nc <= 1) 
    stop("'x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2L) 
    stop("'margins' must be a numeric vector of length 2")
  doRdend <- !identical(Rowv, NA)
  doCdend <- !identical(Colv, NA)
  if (!doRdend && identical(Colv, "Rowv")) 
    doCdend <- FALSE
  if (is.null(Rowv)) 
    Rowv <- rowMeans(x, na.rm = na.rm)
  if (is.null(Colv)) 
    Colv <- colMeans(x, na.rm = na.rm)
  if (doRdend) {
    if (inherits(Rowv, "dendrogram")) 
      ddr <- Rowv
    else {
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      if (!is.logical(Rowv) || Rowv) 
        ddr <- reorderfun(ddr, Rowv)
    }
    if (nr != length(rowInd <- order.dendrogram(ddr))) 
      stop("row dendrogram ordering gave index of wrong length")
  }
  else rowInd <- 1L:nr
  if (doCdend) {
    if (inherits(Colv, "dendrogram")) 
      ddc <- Colv
    else if (identical(Colv, "Rowv")) {
      if (nr != nc) 
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      ddc <- ddr
    }
    else {
      hcc <- hclustfun(distfun(if (symm) 
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      if (!is.logical(Colv) || Colv) 
        ddc <- reorderfun(ddc, Colv)
    }
    if (nc != length(colInd <- order.dendrogram(ddc))) 
      stop("column dendrogram ordering gave index of wrong length")
  }
  else colInd <- 1L:nc
  x <- x[rowInd, colInd]
  labRow <- labRow[rowInd] %||% rownames(x) %||% (1L:nr)[rowInd]
  labCol <- labCol[colInd] %||% colnames(x) %||% (1L:nc)[colInd]
  if (scale == "row") {
    x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 1L, sd, na.rm = na.rm)
    x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
  }
  else if (scale == "column") {
    x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 2L, sd, na.rm = na.rm)
    x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
  }
  lmat <- rbind(c(NA, 3), 2:1)
  lwid <- c(if (doRdend) 1 else lwid)
  lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else lhei)
  if (!missing(ColSideColors)) {
    if (!is.character(ColSideColors) || length(ColSideColors) != 
        nc) 
      stop("'ColSideColors' must be a character vector of length ncol(x)")
    lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
    lhei <- c(lhei[1L], 0.2, lhei[2L])
  }
  if (!missing(RowSideColors)) {
    if (!is.character(RowSideColors) || length(RowSideColors) != 
        nr) 
      stop("'RowSideColors' must be a character vector of length nrow(x)")
    lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                                         1), 1), lmat[, 2] + 1)
    lwid <- c(lwid[1L], 0.2, lwid[2L])
  }
  lmat[is.na(lmat)] <- 0
  if (verbose) {
    cat("layout: widths = ", lwid, ", heights = ", lhei, 
        "; lmat=\n")
    print(lmat)
  }
  dev.hold()
  on.exit(dev.flush())
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1L], 0, 0, 0.5))
    image(rbind(if (revC) 
      nr:1L
      else 1L:nr), col = RowSideColors[rowInd], axes = FALSE)
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2L]))
    image(cbind(1L:nc), col = ColSideColors[colInd], axes = FALSE)
  }
  par(mar = c(margins[1L], 0, 0, margins[2L]))
  if (!symm || scale != "none") 
    x <- t(x)
  if (revC) {
    iy <- nr:1
    if (doRdend) 
      ddr <- rev(ddr)
    x <- x[, iy]
  }
  else iy <- 1L:nr
  image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
          c(0, nr), axes = FALSE, xlab = "", ylab = "", ...)
  axis(1, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
       cex.axis = cexCol)
  if (!is.null(xlab)) 
    mtext(xlab, side = 1, line = margins[1L] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
       cex.axis = cexRow)
  if (!is.null(ylab)) 
    mtext(ylab, side = 4, line = margins[2L] - 1.25)
  if (!missing(add.expr)) 
    eval.parent(substitute(add.expr))
  par(mar = c(margins[1L], 0, 0, 0))
  if (doRdend) 
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  else frame()
  par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))
  if (doCdend) 
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  else if (!is.null(main)) 
    frame()
  if (!is.null(main)) {
    par(xpd = NA)
    title(main, cex.main = 1.5 * op[["cex.main"]])
  }
  invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
                                                              doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}
