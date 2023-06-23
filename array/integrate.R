# Load libraries
library(CEMiTool)
library(rmarkdown)
library(ggplot2)
library(dplyr)
### For the most basic usage of `CEMiTool` only a `data.frame` containing expression data
### with gene symbols in the rows and sample names in the columns is needed
#data(expr0)
#data(sample_annot)
#head(sample_annot)

#' @importFrom gRbase combnPrim
#' @import igraph
#' @GeneOverlap
#'
NULL

#' Generates communities from edgelist 
#'
#' Returns communities from edgelist created by cemoverlap function.
#'
#' @param mod_intersection_df Module intersection dataframe obtained from cemoverlap function
#' @param presence_as_weights Logical. Should sums of node pair occurence be considered edge weights?
#' @param smallest_community Minimal number of genes in community (default:15)
#'
#' @details Function takes edgelist as inputs and generates communities using functions 
#'    provided in igraph package (default:'cluster_fast_greedy:
#'
#' @export 
overlap_community <- function(mod_intersection_df, presence_as_weights = FALSE, smallest_community = 15, 
                              method = 'cluster_fast_greedy'){
  # Available methods:
  # cluster_edge_betweenness // cluster_fast_greedy // cluster_label_prop // cluster_leading_eigen
  # cluster_louvain // cluster_optimal // cluster_spinglass // cluster_walktrap
  edgemat <- as.matrix(mod_intersection_df[,c('gene1', 'gene2')])
  edgegraph <- graph_from_edgelist(edgemat, directed = FALSE)
  if(presence_as_weights){
    weight <- mod_intersection_df[,'proportion']
    edgegraph$weight <- weight
  }
  commfunc <- get(method)
  comm <- commfunc(edgegraph)
  comm <- communities(comm)
  comm <- as.list(comm)
  len_vec <- sapply(comm, length)
  # comm <- comm[order(len_vec, decreasing = T)]
  names(comm) <- paste0('CM', seq_along(1:length(comm)))
  names(comm) <- ifelse(len_vec >= smallest_community, 
                        names(comm), 
                        paste0(names(comm), '.SMALL'))
  # comm <- comm[sapply(comm, length) >= smallest_community]
  return(comm)
}

#' Enriches communities 
#'
#' Returns enrichment of communities from edgelist created by cemoverlap function.
#'
#' @param community_list Community list obtained from overlap community function
#' @param list_of_cem List of Cemitool objects
#' @param list_of_cem_names Character vector to name list of Cemitool objects
#' @param run_fgsea Logical. Should I run fgsea?
#' @param comp_group Which group will be used as base for comparison. If 'none', then all combinations of comparisons will be made
#' @param subject_col Column containing subject information in sample annotation slot of CEMiTool objects 
#'
#' @export 
enrich_mods <- function(community_list, list_of_cem, list_of_cem_names, run_fgsea = F){
  # NOTE: ASSUMPTION of current approach: Relevant modules for a comparison in a study will have a high proportion of differentially
  #         regulated genes in a certain direction. Base assumption is that NON relevant modules will be centered at zero
  # NOTE: Results can be improved in a myriad of ways. We can use information of GMTs to prioritize
  #       modules based on functionality. We can use information inside module slot to prioritize genes with
  #       bigger centrality in individual results. We can cross-reference results with a Fisher comparison approach.
  # NOTE: Think if it is possible to apply approaches such as machine learning or single sample GSEA for stratification.
  # NOTE: Think if it is possible to automatically select p.value or FDR thresholds
  if(!missing(list_of_cem_names)){
    names(list_of_cem) <- list_of_cem_names 
  }else{
    names(list_of_cem) <- paste0('cem_', seq_along(1:length(list_of_cem)))
  }
  # Setting
  mod_exp <- Map(function(cem, cemname){
    # NOTE: Add a subject_column or blocking_column slot to cemitool object 
    # Group and subject are hardcoded
    contmat <- makeContMatrix(samp_annot = cem@sample_annotation, class_column = cem@class_column,
                              comp_group = 'D0', exprs = cem@expression, subject_col = 'subject') 
    toptables <- do.call(makeLimmaComp, contmat)
    iter_comp <- Map(function(comp, compname){
      if(run_fgsea){
        ranks <- rank(comp$logFC)
        names(ranks) <- comp$gene
        fgseaRes <- fgsea::fgsea(pathways = community_list, stats = ranks,
                                 minSize = 15, maxSize = 500, nperm = 1000)
        fgseaRes <- as.data.frame(fgseaRes)
        fgseaRes <- fgseaRes[order(fgseaRes$padj),]
      }else{
        fgseaRes <- data.frame()
      }
      percentage_in_comparison <- Map(function(mod, modname){
        comp <- subset(comp, gene %in% mod & P.Value < 0.05)
        per_up <- sum(comp$logFC > 0)/length(mod)
        per_do <- sum(comp$logFC < 0)/length(mod)
        # Fernando's suggestion
        fcs_pert <- (median(comp$logFC)/sd(comp$logFC)) * median(-log10(comp$P.Value))
        matched_gsea <- match(modname, fgseaRes$pathway)
        if(!is.na(matched_gsea)){
          gsea_vec <- as.numeric(fgseaRes[matched_gsea, 2:7])
        }else{ 
          gsea_vec <- rep(NA, 6)
        }
        scaled <- per_up - per_do
        scaled_percentage <- c(per_up, per_do, scaled, fcs_pert, compname, modname, cemname, paste0(compname, '_in_', cemname), gsea_vec)
        return(scaled_percentage)
      }, community_list, names(community_list))
      percentage_in_comparison <- as.data.frame(do.call(rbind, percentage_in_comparison))
      colnames(percentage_in_comparison) <- c('percentage_up', 'percentage_down', 'scaled_percentage', 'score_mod','comparison', 'community', 
                                              'cem_obj', 'name_comp', 'pval', 'padj', 'ES', 'NES', 'nMoreExtreme', 'size')
      return(percentage_in_comparison)
    }, toptables, names(toptables))
    p_incomp_df <- do.call(rbind, iter_comp)
    return(p_incomp_df)
  }, list_of_cem, names(list_of_cem))
  mod_exp_df <- as.data.frame(do.call(rbind, mod_exp))
  ncolumns <- c('percentage_up', 'percentage_down', 'scaled_percentage', 'score_mod','pval', 'padj', 'ES', 'NES', 'nMoreExtreme', 'size')
  for(i in ncolumns){
    if(i %in% colnames(mod_exp_df)){
      mod_exp_df[[i]] <- as.numeric(mod_exp_df[[i]])  
    }
  }
  mod_exp_df <- Filter(function(x) !all(is.na(x)), mod_exp_df)
  rownames(mod_exp_df) <- NULL
  return(mod_exp_df)
}

#' Integrates CEMiTool analyses
#'
#' Returns the occurrence of edges between different analyses
#'
#' @param ... Objects of class \code{CEMiTool}, data.frames or character string containing
#' the path to a file with genes and modules.
#' @param analyses_list List of objects of class \code{CEMiTool}, data.frames or character
#' strings containing the path to files with genes and modules.
#' @param fraction The fraction of objects in which an edge pair must be present to 
#' be selected (default = 1, accepts values from 0-1)
#'
#' @return Object of class \code{data.frame} containing edgelist describing common 
#' edges between the networks defined in module slots from the input objects
#'
#' @details The method assumes that all genes inside each module are connected to
#' every other gene from the same module. 
#'
#' @examples
#' \dontrun{
#' # Run cemitool twice on expr dataset. In each time, one sample will be removed
#' data(expr0)
#' set.seed(10)
#' dset1 <- expr0[,-sample(1:ncol(expr0), 1)]
#' dset2 <- expr0[,-sample(1:ncol(expr0), 1)]
#' cem1 <- cemitool(dset1) 
#' cem2 <- cemitool(dset2) 
#' cemoverlap_df <- cem_overlap(cem1, cem2)
#' # Can also be run with a list: cemoverlap_df <- cemoverlap(list(cem1, cem2))
#' @export
cem_overlap <- function(..., analyses_list = NULL, fraction = 0, 
                        desired_table = 'adjacency', mc.cores = 1){
  # Desired_table must be one of the three: 
  #   spearman = Return spearman rho
  #   pearson = Return pearson rho
  #   b_correlations = Return adjacency list defined in cemitool object
  #   adjacency = Return discretized edges
  analyses <- c(list(...), analyses_list)
  if(is.null(names(analyses))){
    names(analyses) <- paste0('cem', seq_along(analyses)) 
  }
  analyses <- Filter(Negate(is.null), analyses)
  edgelist <- lapply(seq_along(analyses), function(index){
    cem_obj <- analyses[[index]]
    cem <- module_genes(cem_obj)
    cem_name <- names(analyses[index])
    # Splits by module. Removes  
    mods <- split(cem[,'genes'], cem[,'modules'])
    mods_log <- sapply(mods, length) < 2
    mods <- mods[!mods_log]
    mods <- mods[names(mods) != 'Not.Correlated']
    # combines all genes inside each module
    per_mod <- lapply(mods, function(mod) {
      # Compute pairwise correlations between genes in modules
      if(desired_table %in% c('pearson', 'spearman')){
        mod_cor <- data.table::melt(cor(t(cem_obj@expression[mod,]), method = desired_table))
        mod_gg <- t(apply(mod_cor[,c(1,2)], 1, sort))
        mod_outp <- data.table(mod_gg, mod_cor[,3,drop = T])
        colnames(mod_outp) <- c('gene1', 'gene2', cem_name)
        # Filter out nodes connecting to itself
        mod_outp <- subset(mod_outp, gene1 != gene2)
        # Filter out duplicated nodes
        mod_outp <- subset(mod_outp, !duplicated(mod_outp))
        return(mod_outp)
      }else if(desired_table == 'adjacency'){
        if(mc.cores == 1){
          mod_outp <- do.call(rbind, lapply(gRbase::combnPrim(mod, 2,simplify = F), sort))
        } else{
          mod_outp <- do.call(rbind, 
                              parallel::mclapply(gRbase::combnPrim(mod, 2,simplify = F), 
                                                 sort, mc.cores = mc.cores))
        }
        mod_outp <- data.frame(mod_outp, TRUE)
        colnames(mod_outp) <- c('gene1', 'gene2', cem_name)
        rownames(mod_outp) <- NULL
        return(mod_outp)
      }else if(desired_table == 'b_correlations'){
        adj_mod <- cem_obj@adjacency[mod,mod]
        adj_mod <- melt(adj_mod)
        colnames(adj_mod) <- c('gene1', 'gene2', cem_name)
        # Filter out nodes connecting to itself
        adj_mod <- subset(adj_mod, gene1 != gene2)
        # Filter out duplicated nodes
        adj_mod <- subset(adj_mod, !duplicated(adj_mod))
        return(adj_mod)
      }
    })
    edges <- do.call(rbind, c(per_mod, make.row.names = F))
    return(edges)
  })
  # Merges all studies
  out <- Reduce(function(...){merge(..., by=c('gene1', 'gene2'), all=TRUE)}, edgelist)
  ######setDF(out)
  
  
  # Sum of studies containing pair and order dataframe by sum of occurrences 
  study_names <- colnames(out)[!colnames(out) %in% c('gene1', 'gene2')] 
  if(desired_table %in% c('spearman', 'pearson', 'b_correlations')){
    presentin <- apply(out[,study_names], 1, function(x) length(x) - sum(is.na(x)))
    cor_median <- apply(out[,study_names], 1, function(x) median(x, na.rm = T)) 
    cor_mean <- apply(out[,study_names], 1, function(x) mean(x, na.rm = T)) 
    cor_sd <- apply(out[,study_names], 1, function(x) sd(x, na.rm = T)) 
    out <- out%>%
      mutate(edgeCount = presentin, proportion = edgeCount/length(study_names), 
             edgeCorMedian = cor_median, edgeSd = cor_sd, edgeCorMean = cor_mean)%>%
      # Keep only edges present in at least the proportion of cemitool objects specified in 'fraction' variable
      filter(proportion >= fraction)%>%
      arrange(desc(proportion), desc(edgeCorMedian))%>%
      as.data.frame()
    return(out)
  } else{
    presentin <- apply(out[,study_names], 1, function(x) length(x) - sum(is.na(x)))
    out <- out%>%
      mutate(edgeCount = presentin, proportion = edgeCount/length(study_names))%>%
      # Keep only edges present in at least the proportion of cemitool objects specified in 'fraction' variable
      filter(proportion >= fraction)%>%
      arrange(desc(proportion))%>%
      as.data.frame()
    return(out)
  }
}

# Generate Fisher correspondance between modules
# TODO: Document function
stat_overlap_mods <- function(..., analyses_list, p_thresh = 1, fdr_thresh = 1, jac_thresh = 0,
                              gsea_metric = 'nes'){
  # Analyses_list
  # cutoff for fisher pvalue, fdr or jaccard index (no cutoff as default) 
  # gsea_metric = 'nes'. Metric can als be 'pval' or 'es'
  analyses <- c(list(...), analyses_list)  
  if(is.null(names(analyses))){
    names(analyses) <- paste0('cem',seq_along(analyses)) 
  }
  stopifnot(sum(duplicated(names(analyses))) == 0)
  new_mods_per_cem <- lapply(analyses, function(cem){
    tmpmod <- subset(cem@module, modules != 'Not.Correlated')
    spmod <- split(tmpmod$genes, tmpmod$modules)
    spmod
  })
  new_mods <- unlist(new_mods_per_cem, recursive = F)
  universe <- length(unique(Reduce(union, new_mods)))
  first_lev <- Map(function(fmod_name, firstlev){
    second_lev <- Map(function(smod_name, seclev){
      gene_ovlp <- GeneOverlap::newGeneOverlap(firstlev, seclev, genome.size = universe)
      gene_ovlp <- GeneOverlap::testGeneOverlap(gene_ovlp)
      jac <- slot(gene_ovlp, 'Jaccard')
      fis <- slot(gene_ovlp, 'pval')
      # fmod_len <- paste0(fmod_name, '//', length(firstlev))
      # smod_len <- paste0(smod_name, '//', length(seclev))
      mod_ord <- sort(c(fmod_name, smod_name))
      df_row <- c(mod_ord, jac, fis)
      return(df_row)
    }, names(new_mods), new_mods)
    do.call(rbind, second_lev)
  }, names(new_mods), new_mods)
  df_output <- do.call(rbind, first_lev)
  df_output <- as.data.frame(df_output)
  rownames(df_output) <- NULL
  colnames(df_output) <- c('first_mod', 'second_mod', 'Jaccard', 'Fisherp')
  # remove interactions from the same modules, duplicates and filter out significant interactions
  df_output <- subset(df_output, first_mod != second_mod)
  remove_pair <- duplicated(paste(df_output$first_mod, df_output$second_mod, sep = '.'))
  df_output <- df_output[!remove_pair,]
  df_output$Jaccard <- as.numeric(df_output$Jaccard)
  df_output$Fisherp <- as.numeric(df_output$Fisherp)
  df_output$fdr <- p.adjust(df_output$Fisherp)
  df_output$logp <- -log10(df_output$Fisherp) 
  df_output$logfdr <- -log10(df_output$fdr) 
  # df_output <- df_output[order(df_output$Fisherp),]
  df_output <- df_output[order(df_output$Jaccard, decreasing = T),]
  df_output <- subset(df_output, Jaccard > jac_thresh & Fisherp < p_thresh & fdr < fdr_thresh)
  # Compute node information and study information
  ## return(df_output)
  info_mod <- Map(function(cemname, cem){
    # tmpmod <- subset(cem@module, modules != 'Not.Correlated')
    scores <- Map(function(scrname, scr){
      scr <- scr%>%
        gather(key = class, value = value, -pathway)%>%
        mutate(std = cemname)%>%
        mutate(metric = scrname)%>%
        filter(pathway != 'Not.Correlated')%>%
        unite(module, std, pathway, sep = '.')%>%
        select(module, class, value, metric)%>%
        as.data.frame()
      return(scr)
    }, names(cem@enrichment), cem@enrichment)
    scores <- do.call(rbind, c(scores, make.row.names = F))
    # Subset one metric. Loop was kept in case more metrics are needed.
    scores <- subset(scores, metric == gsea_metric)
    scores <- scores[,!colnames(scores) == 'metric']
    cemdf <- cem@module%>%
      filter(modules != 'Not.Correlated')%>%
      group_by(modules)%>%
      summarise(mod_length = n())%>%
      ungroup()%>%
      mutate(std = cemname)%>%
      unite(module, std, modules, sep = '.')
    cemdf <- merge(cemdf, scores, by = 'module')
    return(cemdf)
  }, names(analyses), analyses)
  info_mod <- do.call(rbind, c(info_mod, make.row.names = F))
  info_mod <- info_mod%>%
    spread(key = class, value = value)%>%
    filter(module %in% c(df_output$first_mod, df_output$second_mod))
  info_mod[is.na(info_mod)] <- 0
  # Determine module activity based on limma
  mod_mean <- Map(function(cemname, cem){
    # Group and subject are hardcoded
    contmat <- makeContMatrix(samp_annot  = cem@sample_annotation, 
                              class_column = cem@class_column,
                              comp_group = 'D0', exprs = cem@expression,
                              subject_col = 'subject') 
    toptables <- do.call(makeLimmaComp, contmat)
    toptables <- do.call(rbind, c(toptables, make.row.names = F))
    # Now determine activity of module
    cemsp <- cem@module%>%
      filter(modules != 'Not.Correlated')
    cemsp <- split(cemsp$genes, cemsp$modules)
    cem_actv <- Map(function(modname, mod){
      tmptop <- toptables%>%
        filter(gene %in% mod)%>%
        group_by(comparison)%>%
        summarise(fc_median = median(logFC), p_median = median(P.Value))%>%
        ungroup()%>%
        mutate(module = paste0(cemname, '.', modname))%>%
        select(comparison, module, fc_median, p_median)%>%
        gather(key = parameter, value = value, fc_median, p_median)
      tmptop 
    }, names(cemsp), cemsp)
    cem_actv <- do.call(rbind, c(cem_actv, make.row.names = F))
    cem_actv <- cem_actv %>%
      unite(new_col, comparison, parameter, sep = '.')%>%
      spread(new_col, value)
    cem_actv
  }, names(analyses), analyses)
  mod_mean <- plyr::rbind.fill(mod_mean)
  # Return three dataframes with aggregation of modules.
  meta_info <- list(module_comparison = df_output, metric_df = info_mod, module_info = mod_mean)
  return(meta_info)
}

# Read and process data for dataset 1 (GSE1297)
matriz_1297 <- read.table('1297/matriz_GSE1297.txt', sep = "\t", header = TRUE, dec = ".")
anot_1297 <- read.table('1297/pdata.txt', sep = "\t", header = TRUE, dec = ".")
anot_1297 <- data.frame(SampleName = anot_1297$GEO, Class = anot_1297$Disease)
cem1297 <- cemitool(matriz_1297, anot_1297, filter = TRUE)

# Read and process data for dataset 2 (GSE5281)
matriz_5281 <- read.table('5281/matriz_GSE5281.txt', sep = "\t", header = TRUE, dec = ".")
anot_5281 <- read.table('5281/pdata.txt', sep = "\t", header = TRUE, dec = ".")
anot_5281 <- data.frame(SampleName = anot_5281$GEO, Class = anot_5281$Disease)
cem5281 <- cemitool(matriz_5281, anot_5281, filter = TRUE) # No modules detected in this dataset

# Read and process data for dataset 3 (GSE28146)
matriz_28146 <- read.table('28146/matriz_GSE28146.txt', sep = "\t", header = TRUE, dec = ".")
anot_28146 <- read.table('28146/pdata.txt', sep = "\t", header = TRUE, dec = ".")
anot_28146 <- data.frame(SampleName = anot_28146$GEO, Class = anot_28146$Disease)
cem28146 <- cemitool(matriz_28146, anot_28146, filter = TRUE)

# Read and process data for dataset 4 (GSE36980)
matriz_36980 <- read.table('36980/matriz_GSE36980.txt', sep = "\t", header = TRUE, dec = ".")
anot_36980 <- read.table('36980/pdata.txt', sep = "\t", header = TRUE, dec = ".")
anot_36980 <- data.frame(SampleName = anot_36980$GEO, Class = anot_36980$Disease)
cem36980 <- cemitool(matriz_36980, anot_36980, filter = TRUE)

# Read and process data for dataset 5 (GSE48350)
matriz_48350 <- read.table('48350/matriz_GSE48350.txt', sep = "\t", header = TRUE, dec = ".")
anot_48350 <- read.table('48350/pdata.txt', sep = "\t", header = TRUE, dec = ".")
anot_48350 <- data.frame(SampleName = anot_48350$GEO, Class = anot_48350$Disease)
cem48350 <- cemitool(matriz_48350, anot_48350, filter = TRUE)

# Read and process data for dataset 6 (GSE84422)
matriz_84422 <- read.table('84422/matriz_GSE84422.txt', sep = "\t", header = TRUE, dec = ".", check.names = FALSE)
anot_84422 <- read.table('84422/pdata.txt', sep = "\t", header = TRUE, dec = ".", check.names = FALSE)
anot_84422 <- data.frame(SampleName = anot_84422$GEO, Class = anot_84422$Disease)
cem84422 <- cemitool(matriz_84422, anot_84422, filter = TRUE)

# Calculate overlap between datasets
cemoverlap_df <- cem_overlap(cem1297, cem28146, cem36980, cem48350, cem84422) # Excluding GSE5281

# Filter and export overlap data
data_out <- data.frame(cemoverlap_df$gene1, cemoverlap_df$gene2, cemoverlap_df$edgeCount)

# Write output files
write.table(data_out, file = "cemoverlap_df.txt", sep = "\t", row.names = FALSE)
