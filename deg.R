#' Run MAST wrapper from Seurat 
#' 
#' @param so Seurat object
#' @param min_expr Numeric between 0 and 1. Proportion of cells that must express a gene for the gene
#' to be included in analysis. Default 0.1 (10%)
#' @param genes_remove Optional character vector. Genes to manually remove from DEG analysis. 
#' @param 
run_deg <- function(so, 
                    min_expr = 0.1, 
                    genes_remove = NULL,
                    cells_keep = NULL,
                    meta_cols_keep = c("pbmc_sample_id", "subject.subjectGuid", "Cohort", "Birth.Year", "Sex"),
                    lrt_vars = c("CohortUP2", "SexMale", "Birth.Year"),
                    fdr_method = "BH",
                    zlm_formula = as.formula("~Cohort + Sex + Birth.Year"),
                    zlm_method = "bayesglm",
                    zlm_ebayes = TRUE,
                    zlm_parallel = TRUE, 
                    ...){
    # Filter cells
    normcounts <- so[["RNA"]]@data
    if(!is.null(cells_keep)){
        assertthat::assert_that(all(cells_keep %in% colnames(normcounts)))
        normcounts <- normcounts[, cells_keep]
    }
    selCells <- colnames(normcounts)
    
    # Filter genes
    selGenes <- data.frame(num_cells_expressed = rowSums(normcounts > 0)) %>%
                           rownames_to_column(var = "Gene") %>%
                           filter(num_cells_expressed >= min_expr*length(selCells))
    if(!is.null(genes_remove)){
        selGenes <- setdiff(selGenes$Gene, genes_remove)
    }
    normcounts <- normcounts[selGenes, ]
    
    # Make SCA object
    fdat <- data.frame(rownames(x = normcounts))
    colnames(x = fdat)[1] <- "primerid"
    rownames(x = fdat) <- fdat[, 1]

    ## Cell Data
    cdat <- so@meta.data %>%
            as.data.frame() %>%
            select(all_of(meta_cols_keep))
    cdat <- cdat[colnames(normcounts), , drop = F]
    assertthat::assert_that(all(rownames(cdat) == colnames(normcounts)))

    # make object
    sca <- MAST::FromMatrix(exprsArray = as.matrix(normcounts),
                            cData = cdat,
                            fData = fdat)
                                
    # Run Hurdle Model
    gc()
    zlm_res <- MAST::zlm(formula = zlm_formula,
                 sca = sca,
                 method = zlm_method,
                 ebayes = zlm_ebayes,
                 parallel = zlm_parallel, ...)  
    gc()
                                
    # Extract convergence info
    df_converged <- zlm_res@converged %>%
        as.data.frame() %>%
        tibble::rownames_to_column("primerid") %>%
        dplyr::rename(Converged_C = C, Converged_D = D)
                                
    # perform LRT
    zlm_summary <- MAST::summary(object = zlm_res, doLRT = lrt_vars)
    gc()

    zlm_dt <- zlm_summary$datatable
                                
    # extract p-values and coefficients
    pvaldf <- zlm_dt %>%
          dplyr::filter(component == "H" &
                        contrast %in% lrt_vars) %>%
          dplyr::select(primerid, contrast, `Pr(>Chisq)`)%>%
            as.data.frame()
    logFCdf <- zlm_dt %>%
           dplyr::filter(component == "logFC" &
                         contrast %in% lrt_vars) %>%
           dplyr::rename(logFC = coef) %>%
           dplyr::select(primerid, contrast, logFC) %>%
           as.data.frame()
    de_df <- pvaldf %>%
        dplyr::left_join(logFCdf) %>%
        dplyr::rename(nomP = `Pr(>Chisq)`) %>%
        dplyr::mutate(adjP = p.adjust(nomP, method = fdr_method)) %>%
        dplyr::mutate(adjP_method = fdr_method) %>%
        dplyr::left_join(df_converged, by = "primerid") %>%
        dplyr::arrange(desc(logFC))

    return(de_df)     
}
                

#' Sample Signal Heatmap
#' Heatmap for plotting genes by mean or median value per sample
#' @param so Seurat Object
#' @param select_features The features (ie genes) to plot
#' @param assay Character value. Name of assay in so to plot. Default "RNA"
#' @param ann_color_list Optional list of color definitions for any annotation metadata in meta_cols.
#' 
sample_deg_heatmap <- function(so, 
                               select_features, 
                               assay = "RNA",
                               ann_color_list = NULL,
                               agg_method = "mean",
                               summarize_by = "pbmc_sample_id", 
                               meta_cols=c("Cohort", "Sex", "BirthYearGroup"),
                               cell_column = "barcodes",
                               cluster_rows = TRUE,
                               cluster_columns = FALSE,
                               show_row_names = TRUE, 
                               show_column_names = TRUE){
                               
    agg_method <- match.arg(agg_method, choices = c("mean", "median"))
    
    counts <- so[[assay]]@data
    meta <- so@meta.data
    
    agg_function <- switch(agg_method,
                           mean = mean,
                           median = median)
    
    # plot heatmap of significant genes: taking mean gene expression per sample
    sigDF <- as.data.frame(counts[select_features, , drop= FALSE])
    mat_df <- sigDF %>%
               rownames_to_column(var = "Gene") %>%
               tidyr::gather(key = !!as.name(cell_column), value = "value", -Gene) %>%
               dplyr::left_join(meta, by = cell_column) %>%
               dplyr::select(all_of(c("Gene","value", summarize_by))) %>%
               group_by(across(all_of(c("Gene", summarize_by)))) %>%
               dplyr::summarize(agg_norm_count = agg_function(value, na.rm = T), .groups = "drop") %>%
               tidyr::spread(key = !!as.name(summarize_by), value = agg_norm_count) %>%
               column_to_rownames(var = "Gene")
        
    # scale matrix and define breaks
    mat_scale <- t(scale(t(mat_df))) %>% 
                as.data.frame()
    iNA <- apply(mat_scale, 1, function(x){all(is.na(x))}) 
    mat_scale[iNA,] <- 0 # if using median, all values can end up as 0 which results in NA scaled values. Forcing all NA to 0 for entire row instead of removing.

    # column annotation of heatmap
    meta_df <- meta %>%
            select(all_of(c(meta_cols, summarize_by))) %>%
            distinct() %>%
            arrange(across(all_of(meta_cols))) %>%
            tibble::remove_rownames() %>%
            tibble::column_to_rownames(summarize_by)

    # Sort matrix columns
    mat_scale <- as.matrix(mat_scale)
    rownames(mat_scale) <- rownames(mat_df)
    mat_scale <- mat_scale[, rownames(meta_df), drop = FALSE]
       
    # define custom color palatte
    max_val <- max(abs(mat_scale))
    if(max_val > 3){
        colorPalette <- circlize::colorRamp2(breaks = c(-1*max_val,-3, 0, 3, max_val), colors=c("darkorchid1","purple", "black", "yellow","lawngreen"))
    } else {
        colorPalette <- circlize::colorRamp2(breaks = c(-3, 0, 3), colors=c("purple", "black", "yellow"))
    }
    
    # plot
    ha <- ComplexHeatmap::HeatmapAnnotation(df = meta_df, col = ann_color_list)
    ComplexHeatmap::Heatmap(matrix = mat_scale,
                            col = colorPalette,
                          cluster_columns = cluster_columns,
                          cluster_rows = cluster_rows,
                          show_row_names = show_row_names,
                          show_column_names = show_column_names,
                          top_annotation = ha)
}   




#' Feature Heatmap for Sampled Cells

#' ComplexHeatmap for plotting genes by cells, sampling a set number of cells per sample
#' using an input list of genes and a Seurat object. Will take results from the data slot
#' of the specified assay.

#' @param select_features Character vector of features (ie genes) to plot
#' 

subsample_deg_heatmap <- function(so, 
                               assay = "RNA",
                               select_features, 
                               ann_color_list = NULL,
                               n_cells = 10,
                               summarize_by = "pbmc_sample_id", 
                               meta_cols=c("Cohort", "Sex", "BirthYearGroup"),
                               cell_column = "barcodes",
                               cluster_rows = TRUE,
                               cluster_columns = FALSE,
                               show_row_names = TRUE, 
                               show_column_names = TRUE,
                               seed.use = 3){
    set.seed(seed.use)
    
    # sample cells based on metadata
    meta <- so@meta.data %>%
        mutate(temp_rownames = rownames(.)) %>%
        group_by(across(all_of(summarize_by))) %>%
        sample_n(size = min(n_cells, n()), replace = FALSE) %>%
        ungroup() %>%
        tibble::column_to_rownames("temp_rownames")

    # subset count data by selected cells and genes
    counts <- so[[assay]]@data
    mat_df <- as.data.frame(counts[select_features, rownames(meta), drop=FALSE])

    # scale matrix 
    mat_scale <- t(scale(t(mat_df))) %>% 
                as.data.frame()

    # column annotation of heatmap
    meta_df <- meta %>%
            select(all_of(c(meta_cols, summarize_by,cell_column))) %>%
            arrange(across(all_of(c(meta_cols, summarize_by, cell_column)))) %>%
            select(-!!as.name(cell_column))
    head(meta_df)

    # Sort matrix columns
    mat_scale <- as.matrix(mat_scale)
    rownames(mat_scale) <- rownames(mat_df)
    mat_scale <- mat_scale[, rownames(meta_df), drop = FALSE]

    # define custom color palatte
    max_val <- max(abs(mat_scale))
    if(max_val > 3){
        colorPalette <- circlize::colorRamp2(breaks = c(-1*max_val,-3, 0, 3, max_val), colors=c("darkorchid1","purple", "black", "yellow","lawngreen"))
    } else {
        colorPalette <- circlize::colorRamp2(breaks = c(-3, 0, 3), colors=c("purple", "black", "yellow"))
    }
    
    # plot
    ha <- ComplexHeatmap::HeatmapAnnotation(df = meta_df, col = ann_color_list)
    ComplexHeatmap::Heatmap(matrix = mat_scale, 
                            col = colorPalette,
                          cluster_columns = cluster_columns,
                          cluster_rows = cluster_rows,
                          show_row_names = show_row_names,
                          show_column_names = show_column_names,
                          top_annotation = ha)
}


#' Plot a volcano plot of input deg data frame, with color annotations based on
#' p-value and effect-size cutofffs
plot_volcano <- function(dat, 
                         es_col, 
                         p_col, 
                         p_cutoff = 0.05, 
                         es_cutoff = 0.2, 
                         title = "", 
                         feat_label = NULL, 
                         feat_col = NULL, 
                         labelsize = 10,
                         alpha = 1,
                        max.overlaps=10){
    assign_change <- function(p, p_cutoff, es, es_cutoff){
        change <- rep("", length(p))
        iUp <- which(p < p_cutoff & es > es_cutoff)
        if(length(iUp) > 0){change[iUp] <- "Up"}
        iDown <- which(p < p_cutoff & es < -1*es_cutoff)
        if(length(iDown) > 0){change[iDown] <- "Down"}
        return(change)
    }
    
    temp <- dat %>%
        mutate(Change = assign_change(get(p_col), p_cutoff, get(es_col), es_cutoff))
    
    g_volcano <- ggplot(temp, aes(!!as.name(es_col), -log10(!!as.name(p_col)), col = Change)) +
        geom_point(alpha = alpha) +
        scale_color_manual(values = c("blue", "red", "gray"), breaks = c("Down", "Up","")) +
        geom_hline(yintercept = -log10(p_cutoff), color = "black", linetype = 2) +
        ggtitle(title, subtitle = sprintf("cutoff.p : %s, cutoff.effectsize: %s", p_cutoff, es_cutoff))
    
    if(!is.null(feat_label)){
        if(is.null(feat_col) || !feat_col %in% names(temp)) {
            warning("Could not label features without a valid feat_col value")
        }
        temp <- temp %>%
            mutate(plotlabel = ifelse(get(feat_col) %in% feat_label, get(feat_col), ""))
        g_volcano <- g_volcano +
            ggrepel::geom_text_repel(data = temp, color = "black", size = labelsize, aes(label = plotlabel), max.overlaps = max.overlaps)
        
    }
                           
    return(g_volcano)
}



convert_upset <- function(upset_plot, rel_heights = c(3,1), rel_widths = c(2,3)){
    # Converts upset plot object to a format that can be plotted with cowplot or patchwork multi-plot functions
    # based on this code https://github.com/hms-dbmi/UpSetR/issues/105
    converted <- cowplot::plot_grid(NULL, upset_plot$Main_bar, upset_plot$Sizes, upset_plot$Matrix,
                            nrow=2, align='hv', rel_heights = rel_heights,
                           rel_widths = rel_widths)
    return(converted)
}


#' Make an upset plot from a data frame
#'
#' Make an upset plot from an input deg dataframe and given p-value and effect size cutoff. 
#' 
#' @param deg_df A data.frame containing DEG values. Data should be prefiltered to the relevant data, but the 
# Data frame should be pre-filtered to relevant observations, and the group_col identifies the 
# metadata column across which to compare DEGs.
plot_upset <- function(deg_df, 
                       p_cutoff, 
                       es_cutoff, 
                       p_col = "adjP", 
                       es_col = "logFC", 
                       genecol = "primerid", 
                       group_col = "contrast",
                       convert_grob = TRUE,
                       rel_heights = c(3,1), 
                       rel_widths = c(2,3),
                      ...){
    us_df <- deg_df %>%
        dplyr::filter(get(p_col) <= p_cutoff, (abs(get(es_col)) >= es_cutoff)| is.na(get(es_col)) & es_cutoff == 0) %>%
        dplyr::select(all_of(c(genecol, group_col))) %>%
        dplyr::mutate(has_gene = 1) %>%
        tidyr::spread(key = group_col, value = has_gene) %>%
        tibble::column_to_rownames(genecol)
    us_df[is.na(us_df)] <- 0
    
    comps <- unique(deg_df[[group_col]])
    
    plt <- UpSetR::upset(us_df, sets = intersect(comps, colnames(us_df)), text.scale = 2,...)
    if(convert_grob){
        plt <- convert_upset(plt, rel_heights=rel_heights, rel_widths= rel_widths)
    } 
    return(plt)
}