hm_to_grob <- function(hm){
    grid::grid.grabExpr(ComplexHeatmap::draw(hm))
} 

subsample_deg_heatmap <- function(so, 
                               features, 
                               ann_color_list = NULL,
                               n_cells = 10,
                               summarize_by = "pbmc_sample_id", 
                               meta_cols=c("Cohort", "Sex", "BirthYearGroup"),
                               cell_column = "barcodes",
                               cluster_rows = TRUE,
                               cluster_columns = FALSE,
                               show_row_names = TRUE, 
                               show_column_names = TRUE,
                               assay = "RNA",
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
    mat_df <- as.data.frame(counts[features, rownames(meta), drop=FALSE])

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


make_cluster_hm_meta <- function(meta, 
                                 value = "count",
                                 cluster_col,
                                 meta_col, 
                                 ann_meta=NULL, 
                                 legend_ann_meta = NULL,
                                 color_ann_meta=NULL,
                                 bScale = TRUE, 
                                 plotname = NULL, 
                                 outputGrob=TRUE, 
                                 color_matrix_fun = NULL,...){
    df_counts <- meta %>%
        group_by(across(all_of(c(cluster_col, meta_col))))%>%
        tally() 
    if(value == "count"){
        df_counts<- df_counts %>%
            tidyr::spread(key = all_of(meta_col), value = n) %>%
            tibble::column_to_rownames(cluster_col) %>%
            mutate_all(tidyr::replace_na, 0)
        hm_name <- "Z-score"
    } else if (value == "percent"){
        df_counts<- df_counts %>%
             group_by(across(all_of(c(meta_col)))) %>%
            mutate(pct= n/sum(n)*100) %>%
            select(-n) %>%
            tidyr::spread(key = all_of(meta_col), value = pct) %>%
            tibble::column_to_rownames(cluster_col) %>%
            mutate_all(tidyr::replace_na,0)
        hm_name <- "Percent"
    }
    mat_counts <- as.matrix(df_counts)
    
    if(bScale){
        mat_counts <- t(apply(mat_counts, 1, scale, scale=TRUE, center = TRUE))
        colnames(mat_counts) <- colnames(df_counts)
        rownames(mat_counts) <- rownames(df_counts)
    }
    
    if(is.null(plotname)){
        plotname <- hm_name
    }
    
    if(!is.null(ann_meta)) {
        meta_ann <- meta %>%
            select(all_of(unique(c(meta_col, ann_meta)))) %>%
            distinct() %>%
            remove_rownames() %>%
            arrange(across(all_of(ann_meta))) %>%
            as.data.frame() %>%
            tibble::column_to_rownames(meta_col)
        mat_counts <- mat_counts[ , rownames(meta_ann)]
        
        if(is.null(legend_ann_meta)){
            legend_ann_meta <- rep(TRUE, length(meta_ann))
        } 
        if(is.null(color_ann_meta)){
            color_ann_meta <- lapply(ann_meta, function(x){
                    meta_levels <- unique(unlist(meta_ann[,x, drop = TRUE]))
                    n_levels <- length(meta_levels)
                    cl_vals <- H5weaver::varibow(n_levels)
                    names(cl_vals) <- meta_levels
                    cl_vals
                })
            names(color_ann_meta) <- ann_meta
        }
        
        ha <- ComplexHeatmap::HeatmapAnnotation(df = meta_ann, show_legend = legend_ann_meta, col = color_ann_meta)
        
        if(!is.null(color_matrix_fun)){
#             pt_99 <- quantile(abs(mat), 0.99)
#             color_matrix_fun <- color_matrix_fun(asseq(-1*pt_99, pt_99), length.out = 100)
            hm <- ComplexHeatmap::Heatmap(mat_counts, show_row_names = TRUE, 
                                          name = plotname, top_annotation = ha, 
                                          col = color_matrix_fun(as.vector(mat_counts)),...)
        } else {
            hm <- ComplexHeatmap::Heatmap(mat_counts, show_row_names = TRUE, 
                                          name = plotname, top_annotation = ha, ...)
        }

    } else {
        if(!is.null(color_matrix_fun)){
#             pt_99 <- quantile(abs(mat), 0.99)
#             color_matrix_fun <- color_matrix_fun(seq(-1*pt_99, pt_99), length.out = 100)
            hm <- ComplexHeatmap::Heatmap(mat_counts, show_row_names = TRUE, name = plotname, 
                                          col = color_matrix_fun(as.vector(mat_counts)),...)
        } else {
            hm <- ComplexHeatmap::Heatmap(mat_counts, show_row_names = TRUE, name = plotname, ...)
        }
    }
    
    if(outputGrob){
      return(hm_to_grob(hm))
    } else {
      return(hm)
    }
    
}


# Bar plot within cluster of a metadata variable. 
cluster_bar_plot <- function(metadata, clustercol, clustercol_value=NULL, 
                             xgroup, color_col = NULL, bargroup = NULL, plotvalue = "count", med_line=FALSE, ...){
    sum_table_pct <- metadata %>%
        as.data.frame() %>%
        group_by(across(all_of(c(clustercol, xgroup, bargroup,color_col)))) %>%  
        summarize(N = n(), .groups="drop") %>%
        group_by(across(all_of(xgroup))) %>% 
        mutate(Total_N = sum(N),
               Percent = N/Total_N*100) %>%
        ungroup() 
    
    plotval <- switch(plotvalue,
                  count = "N",
                  percent = "Percent")
    medline_val <- switch(plotvalue,
                  count = "N_med",
                  percent = "Percent_med")
    
    if(!is.null(clustercol_value)) {
        sum_table_pct <- sum_table_pct %>%
            filter(get(clustercol) == clustercol_value) %>%
            arrange(across(all_of(c(xgroup, bargroup))))
    }
    if(is.null(color_col)){
        color_col = xgroup
    }
    sum_table_pct[[color_col]] <- factor(sum_table_pct[[color_col]])
    
    if(!is.null(bargroup)){
        g_cluster_bar <- ggplot(sum_table_pct, aes(!!as.name(xgroup), !!as.name(plotval), fill = !!as.name(color_col), group = !!as.name(bargroup)))
    } else {
        g_cluster_bar <- ggplot(sum_table_pct, aes(!!as.name(xgroup), !!as.name(plotval), fill = !!as.name(color_col)))
    }
    
    g_cluster_bar <- g_cluster_bar +
        geom_col(color = "black") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    if(is.null(clustercol_value)){
        g_cluster_bar <- g_cluster_bar + facet_wrap(as.formula(paste0("~",clustercol)), ...)
    }
    
    g_cluster_bar
}



proportion_cluster_trend <- function(so, 
                                     cluster_col, 
                                     time_col= "Day", 
                                     subject_col = "subject.subjectGuid", 
                                     shape_col = "Sex", 
                                     color_col = "Cohort",
                                     title="Trending Percent Cells per Cluster by Subject over Time",
                                     output_table = FALSE,
                                     ...){
  df_summary <- suppressWarnings(so@meta.data %>%
        group_by(across(all_of(c(subject_col, shape_col, time_col,color_col, cluster_col)))) %>%
        summarize(N = n(), .groups = "drop") %>%
        group_by(across(all_of(c(subject_col, shape_col, time_col,color_col)))) %>%
        mutate(Pct = N/sum(N)*100) %>%
        ungroup() %>%
        group_by(across(all_of(subject_col))) %>%
        mutate(ClusterFoldChange = Pct/Pct[Day==0]) %>%
        ungroup())
    gtrend <-  ggplot(df_summary, aes(factor(!!as.name(time_col)), Pct, group = !!as.name(subject_col), 
                                      color = !!as.name(color_col))) +
            geom_point(aes(shape = !!as.name(shape_col)))+
            geom_line() +
            facet_wrap(as.formula(paste0("~",cluster_col)), scales = "free_y",...) +
            ggtitle(title)
    if(output_table){
      return(list(gtrend,df_summary))
    } else{
      return(gtrend)
    }

}

proportion_cluster_group_box <- function(so, 
                                     cluster_col, 
                                     time_col= "Day", 
                                     subject_col = "subject.subjectGuid", 
                                     group_col = c("Cohort","Sex"),
                                     color_col = "Cohort",
                                     title="Trending Percent Cells per Cluster by Subject over Time",
                                     output_table = FALSE,
                                     ...){
   summary_df <- suppressWarnings(so@meta.data %>%
        group_by(across(all_of(unique(c(subject_col, time_col,group_col,color_col, cluster_col)))))%>%
        summarize(N = n(), .groups = "drop") %>%
        group_by(across(all_of(c(subject_col,group_col)))) %>%
        mutate(Pct = N/sum(N)) %>%
        ungroup() %>%
        group_by(across(all_of(subject_col))) %>%
        mutate(ClusterFoldChange = Pct/Pct[get(time_col)==0]) %>%
        ungroup()%>%
        mutate(Log2FC = log2(ClusterFoldChange))) %>%
        filter(Day !=0) %>%
        unite(col = "group", all_of(group_col), remove = F)
    
    g_foldchange <- ggplot(summary_df, aes(group, Log2FC, color = !!as.name(color_col)), label = !!as.name(subject_col)) +
        geom_point(position=position_jitter(height = 0, width = 0.3, seed = 3)) +
        geom_boxplot(alpha = 0, outlier.alpha = 1, color = "black") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        facet_wrap(...)
    
    if(output_table){
        return(list(g_foldchange, summary_df))
    } else {
        return(g_foldchange)
    }
}  


cluster_marker_hm <- function(so, cluster_col, cluster_df=NULL, genes = NULL, plottitle=NULL,
                              n_cells_cluster=50, n_genes = 10, output_table = FALSE, assay = "RNA"){
    clusters <- unique(cluster_df$cluster)
    if(is.null(genes)){
      genes <- character()
      for(clust in clusters){
          genes_cl <- cluster_df %>%
                  filter(cluster == clust) %>%
                  filter(p_val_adj <0.05) %>%
                  filter(!gene %in% genes) %>%
                  mutate(direction = sign(avg_log2FC)) %>%
                  group_by(direction, cluster) %>%
                  top_n(n_genes, wt = abs(avg_log2FC)) %>%
                  ungroup() %>%
                  distinct(gene) %>%
                  pull(gene)
          genes <- c(genes, genes_cl)
      }
      if(is.null(plottitle)){
          plottitle <- sprintf("Top %s High and Low Significant Markers per Cluster", n_genes)
      }
    } else {
      if(is.null(plottitle)){
          plottitle <- ""
      }
    }
  
    
    so <- Seurat::SetIdent(so, value = so@meta.data[[cluster_col]])
    so <- subset(so, downsample = n_cells_cluster)

    hm_genes <- DoHeatmap(so, features = genes, assay = assay) & 
        ggtitle(plottitle)
    
    if(output_table & !is.null(cluster_df)){
      return(list(hm_genes, cluster_df))
    } else{
      return(hm_genes)
    }
}


# QC Summary Table for Clusters
  summarize_cluster_qc <- function(so, cluster_col){
    assertthat::assert_that("pct.mt" %in% names(so@meta.data),
                            msg = "Cannot find 'pct.mt' column in cell metadata")
    summary_df <- so@meta.data %>%
        group_by(across(all_of(cluster_col))) %>%
        summarize(N_cells = n(),
                  Median_Genes = median(n_genes),
                  Median_UMI = median(n_umis),
                  Median_Pct_MitoUMI = median(pct.mt), .groups = "drop")
    summary_df
  }
  
  # Heatmap by cluster
  subsample_deg_cluster_heatmap <- function(so, 
                               features, 
                               ann_color_list = NULL,
                               n_cells = 10,
                               summarize_by = "pbmc_sample_id", 
                               meta_cols=c("Cluster"),
                               cell_column = "barcodes",
                               cluster_rows = TRUE,
                               cluster_columns = TRUE,
                               show_row_names = TRUE, 
                               show_column_names = TRUE,
                               assay = "RNA",
                               seed.use = 3){
    set.seed(seed.use)
    
    assertthat::assert_that(length(meta_cols) == 1)
    
    # sample cells based on metadata
    meta <- so@meta.data %>%
        mutate(temp_rownames = rownames(.)) %>%
        group_by(across(all_of(summarize_by))) %>%
        sample_n(size = min(n_cells, n()), replace = FALSE) %>%
        ungroup() %>%
        tibble::column_to_rownames("temp_rownames")

    # subset count data by selected cells and genes
    counts <- so[[assay]]@data
    mat_df <- as.data.frame(counts[features, rownames(meta), drop=FALSE])

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
                          column_split = meta_df[[meta_cols]],
                          show_row_names = show_row_names,
                          show_column_names = show_column_names,
                          top_annotation = ha)
}


#' Create a 
qc_plot_seurat_clusters <- function(so, 
                                   reduction_name, 
                                   cluster_col, 
                                   ncol = NULL,
                                   sizes.highlight = 1,
                                   cols.highlight = "#DE2D26"){
    clusters <- sort(as.integer(as.character(unique(so@meta.data[[cluster_col]]))))
    plot_list <- lapply(clusters, function(cl){
        clust_cells <- list(so$barcodes[so[[cluster_col]] == as.character(cl)])
        names(clust_cells)[1] <- paste0("Cluster ",cl)
        DimPlot(so, 
                reduction = reduction_name,
                cells.highlight = clust_cells, 
                ncol = ncol, 
                cols.highlight=cols.highlight, 
                sizes.highlight=sizes.highlight) + 
        ggtitle(paste0("Cluster ",cl))
    })
    patchwork::wrap_plots(plot_list)
}