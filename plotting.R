
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
                         max.overlaps = 10,
                         base_size = 12,
                         themels = list()){
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
        theme_bw(base_size = base_size) +
        ggtitle(title, subtitle = sprintf("cutoff.p : %s, cutoff.effectsize: %s", p_cutoff, es_cutoff)) +
        do.call(theme, themels)
    
    if(!is.null(feat_label)){
        if(is.null(feat_col) || !feat_col %in% names(temp)) {
            warning("Could not label features without a valid feat_col value")
        }
        temp <- temp %>%
            mutate(plotlabel = ifelse(get(feat_col) %in% feat_label, get(feat_col), ""))
        g_volcano <- g_volcano +
            ggrepel::geom_text_repel(data = temp, color = "black", size = labelsize, 
                                     max.overlaps = max.overlaps, aes(label = plotlabel))
        
    }
                           
    return(g_volcano)
}



convert_upset <- function(upset_plot, rel_heights = c(3,1), rel_widths = c(2,3), title = "", titlesize = 10, title_rel_heights = c(1,10), has_meta){
    # Converts upset plot object to a format that can be plotted with cowplot or patchwork multi-plot functions
    # based on this code https://github.com/hms-dbmi/UpSetR/issues/105
    if(has_meta){
        converted <- cowplot::plot_grid(NULL, NULL,upset_plot$Main_bar, upset_plot$Sizes, upset_plot$set.metadata.plots[[1]], upset_plot$Matrix,
                            nrow=2, align='hv', 
                                        rel_heights = rel_heights,
                           rel_widths = rel_widths)
        
    } else {
        converted <- cowplot::plot_grid(NULL,upset_plot$Main_bar, upset_plot$Sizes, upset_plot$Matrix,
                            nrow=2, align='hv', rel_heights = rel_heights,
                           rel_widths = rel_widths)
    }
    
    if(title != ""){
     title_theme <- cowplot::ggdraw() +
          cowplot::draw_label(title, x = 0.05, hjust = 0, size = titlesize)   
     converted <- cowplot::plot_grid(title_theme, converted,
                            nrow=2, align='hv', rel_heights = title_rel_heights)
    }
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
                       title = "", 
                       rel_heights = c(3,1), 
                       rel_widths = c(2,3),
                       titlesize = 20,
                       textscale = 2, 
                       set_colors = NULL,
                       set_meta = NULL,
                       ...){
    us_df <- deg_df %>%
        dplyr::filter(get(p_col) <= p_cutoff, (abs(get(es_col)) >= es_cutoff)| is.na(get(es_col)) & es_cutoff == 0) %>%
        dplyr::select(all_of(c(genecol, group_col))) %>%
        dplyr::mutate(has_gene = 1) %>%
        tidyr::spread(key = group_col, value = has_gene) %>%
        tibble::column_to_rownames(genecol)
    us_df[is.na(us_df)] <- 0
    
    comps <- unique(deg_df[[group_col]])
    sets <- intersect(comps, colnames(us_df))
    if(!is.null(set_meta)){
        set_meta <- set_meta[set_meta[[1]] %in% sets, ]
        set_colors<- set_colors[names(set_colors) %in% set_meta$Group]
        meta_list <- list(data = set_meta, 
                         plots = list(list(type = "heat", 
                                          column = "Group", 
                                           assign = 6,
                                          colors = set_colors)))
#                                            colors = set_colors)))

        plt <- UpSetR::upset(us_df, sets = sets, text.scale = textscale, set.metadata = meta_list,...)
#         set_colors
#         bar_color <- set_colors[sets] 
#         bar_color <- bar_color[order(colSums(us_df), decreasing = T)]
    } else {
#         bar_color <- "gray23"  
        plt <- UpSetR::upset(us_df, sets = sets, text.scale = textscale,...)

    }

#     plt <- UpSetR::upset(us_df, sets = sets, text.scale = textscale, sets.bar.color  = bar_color, ...)
    if(convert_grob){
        plt <- convert_upset(plt, rel_heights, rel_widths, title = title, titlesize = titlesize, has_meta = !is.null(set_meta))
    } 
    return(plt)
}



# this should be used on subsetted count matrix (ie specific genes)--don't use on full gene count matrix
#' Rescale a Normalized count Matrix, excluding 0's
#'
#' @param norm_mat
#' @param margin Numeric value, 1 or 2. Default 1 will scale within row (Seurat matrix = within gene)
#' @param replace_val Value that 0's should be replaced with in result matrix, default NA.
#' @return A Dense matrix of scaled with 0's replaced with specified value.
#' @examples
#' set.seed(3)
#' datamat <- matrix(rpois(n = 25, lambda = 2), nrow = 5)
#' rownames(datamat) <- paste0("gene",LETTERS[1:5])
#' colnames(datamat) <- paste0("cell", 1:5)
#' print(datamat)
#' print(scale_mat_exclude_0)
scale_mat_exclude_0 <- function(norm_mat, margin = 1, replace_val = NA){
    if(nrow(norm_mat) > 100){
        warning("Scaling function intended for use on subsetted matrix. Attempting to create a dense matrix with more than 100 rows.")
    }
    assertthat::assert_that(!is.null(rownames(norm_mat)), msg = "matrix must have row names")
    assertthat::assert_that(!is.null(colnames(norm_mat)), msg = "matrix must have column names")
    
    norm_mat[norm_mat==0] <- NA
    df_scale <- t(
        apply(norm_mat, margin, function(x){
            scale_val <- rep(NA, length(x)); 
            names(scale_val) <- names(x);
            scaled <- scale(x[!is.na(x)]); 
            scale_val[rownames(scaled)] <- scaled; 
            scale_val
       })
    )
    assertthat::assert_that(all(dim(df_scale) == dim(norm_mat)))
    
    return(df_scale)

}


#' Get the clustered row or column order based on a count matrix and clustering on aggregated
#' data within provided metadata level combinations
#' @param features The gene names to cluster on
#' @param group_vars The metadata variables to aggregate marker information by
#' @param so A seurat object, default NULL. If not supplied must supply meta and mat
#' @param meta Optional. A metadata matrix of cells by metadata. Can be supplied instead of a Seurat object
#' @param mat Optional. A count matrix (normalized counts) of genes by cells. Can be supplied instead of a Seurat object
#' @param out_dim Dimension ("markers" or "meta") to output clustered order for
#' @param assay Default "RNA". Assay from Seurat object to pull counts from
#' @param slot Default "data". Slot from Seurat assay to pull counts from
#' @param aggmethod Default mean. currentluy only choice is mean. How to aggregate data within metadata groups prior to clustering
order_dim <- function(features, 
                     group_vars,
                     so = NULL,
                     meta = NULL,
                     mat = NULL,
                     out_dim = c("markers", "meta"),
                     assay="RNA", 
                     slot = "data",
                     aggmethod = "mean"){
    
    out_dim <- match.arg(out_dim,  choices = c("markers", "meta"), several.ok = FALSE)
    aggmethod <- match.arg(aggmethod,  choices = c("mean"), several.ok = FALSE)
     
    if(aggmethod == "mean"){
        # transform scaled data to long then average by metadata groupings
        if(!is.null(so)){
            mat <- Seurat::GetAssayData(so, assay = assay, slot = slot)[features, , drop=FALSE]
            meta <- so@meta.data
        } else {
            assertthat::assert_that(ncol(mat) == nrow(meta))
        }
        df_scaled_no0 <- scale_mat_exclude_0(norm_mat =  mat, margin = 1, replace_val = NA)
        
        avg_mat <- meta %>%
            dplyr::select(all_of(unique(group_vars))) %>%
            cbind(t(as.matrix(df_scaled_no0,nrow = nrow(df_scaled_no0), ncol = ncol(df_scaled_no0)))) %>%
            tidyr::gather(key = "marker", value = "value", all_of(features)) %>%
            group_by(across(all_of(c(group_vars, "marker")))) %>%
            dplyr::summarize(avg = mean(value, na.rm = T), .groups = "drop") %>%
            spread(key = "marker", value = "avg") %>%
            tidyr::unite(col = "newrows",all_of(unique(group_vars))) %>%
            tibble::column_to_rownames("newrows")
        rm(df_scaled_no0)
        mat <- avg_mat[,features, drop = FALSE]
        if(out_dim == "markers"){
            clust_res <- hclust(dist(t(mat)))
        } else if (out_dim == "meta"){
            clust_res <- hclust(dist(mat))
        }
        dimorder <- clust_res$labels[clust_res$order]
    } 
    
    return(dimorder)

}

# order_dim(so, 
#           features=top20_len1, 
#           group_vars = "wsnn_f_res.1", 
#            out_dim = "meta",
#          assay="RNA", 
#          slot = "scale.data",
#          aggmethod = "mean")

#' Plot a dotplot from a seurat object or metadata + counts
#'
#'
myDotPlot <- function(features, 
                      column_var, 
                      so = NULL, 
                      # column_meta=NULL, 
                      column_group_var = NULL,
                      split_col_var=NULL, 
                      split_row_var=NULL, 
                      assay="RNA", 
                      slot = "scale.data",
                      base_size = 14,
                      value_limit = 2.5,
                      clip_pct = 0.99, 
                      color_vals = viridis::turbo(50),  # Seurat::PurpleAndYellow(k = 50)
                      marker_order = NULL,
                      split_col_order = NULL,
                      cluster_markers = TRUE,
                      cluster_split_cols = FALSE){
    
    meta <- so@meta.data %>%
        dplyr::select(all_of(unique(c(column_var, split_col_var, split_row_var, column_group_var))))
       
    # transform count data to long then calc percent expressed by metadata groupings
    df_counts <- Seurat::GetAssayData(so, assay = assay, slot = "counts")[features, , drop=FALSE]
    pct_df <- meta %>%
        cbind(t(as.matrix(df_counts))) %>%
        tidyr::gather(key = "marker", value = "value", all_of(features)) %>%
        dplyr::group_by(across(all_of(c(column_var, split_col_var, column_group_var, "marker")))) %>%
        dplyr::summarize(N = dplyr::n(),
                  pct = sum(value > 0)/N*100, .groups = "drop")
    
    # scale normdata excluding 0's, then average by metadata groupings
    df_norm <- Seurat::GetAssayData(so, assay = assay, slot = "data")[features, , drop=FALSE]
    df_scaled_no0 <- scale_mat_exclude_0(norm_mat =  df_norm, margin = 1, replace_val = NA)
    avg_df <- meta %>%
        cbind(t(as.matrix(df_scaled_no0))) %>%
        tidyr::gather(key = "marker", value = "value", all_of(features)) %>%
        group_by(across(all_of(c(column_var, split_col_var, column_group_var, "marker")))) %>%
        dplyr::summarize(avg = mean(value, na.rm = T), .groups = "drop")
    rm(df_scaled_no0)
    rm(df_norm)
    
    # Determine marker order
    if(!is.null(marker_order)){
        if(!all(features %in% marker_order)){
            warning(sprintf("Missing features in input marker_order [%s]. Reverting to default order.",
                            paste(setdiff(features, marker_order)), collapse=","))
            marker_order <- sort(unique(features))
        }
    } else if(cluster_markers){
        
        marker_order <- order_dim(so,
                             out_dim = "markers",
                             features= features, 
                             group_vars = c(column_var, split_col_var, column_group_var),
                             assay=assay, 
                             slot = "data",
                             aggmethod = "mean")
    } else {
        marker_order <- sort(unique(features)) 
    }
        
    # Facet Order
    if(!is.null(split_col_var) & !is.null(split_col_order)){
        split_lev <- unique(so@meta.data[[split_col_var]])
        if(!all(split_lev %in% split_col_order)){
            warning(sprintf("Missing column values in input split_col_order [%s]. Reverting to default order.",
                            paste(setdiff(split_lev, split_col_order)), collapse=","))
            split_col_order <- sort(split_lev)
        }
    } else if (cluster_split_cols){
        split_col_order <- order_dim(so,
                             features, 
                             group_vars = split_col_var,
                             out_dim = "meta",
                             assay = assay, 
                             slot = "data",
                             aggmethod = "mean")
    } 
        
    # set limit on extreme values to plot color scale nciely--use min range of input percentile or specific value limit
    max.val <- min(quantile(abs(avg_df$avg), clip_pct, na.rm = TRUE), value_limit, na.rm = TRUE)

    temp <- avg_df %>%
        dplyr::mutate(avg= ifelse(avg > max.val, max.val,avg)) %>%  # clip ends
        dplyr::mutate(avg=ifelse(avg < -1*max.val,-1* max.val,avg)) %>% # clip ends
        dplyr::full_join(pct_df, by = intersect(colnames(avg_df),colnames(pct_df))) %>%
        mutate(marker = factor(marker, levels = marker_order))
    
    if(!is.null(split_col_var) & !is.null(split_col_order)){
        temp[[split_col_var]] <- factor(as.character(temp[[split_col_var]]), levels = split_col_order)
    }
    
    g <- ggplot(temp, aes(!!as.name(column_var), marker, color = avg, size = pct)) +
        scale_size_continuous(limits = c(0, 101), range =c(0,7))+
        theme_bw(base_size = base_size)
    if(is.null(column_group_var)){
        g <- g + geom_point()
    } else {
        g <- g + geom_point(aes(group = !!as.name(column_group_var)), position = position_dodge(width = 0.7))
    }

    if(!is.null(split_col_var)){
        g <- g +
            facet_grid(as.formula(paste0("~", split_col_var)), scales = "free_x", space = "free")
    }
    
    g + 
        scale_color_gradientn(
            limits = c(-1*max.val, max.val), 
                       colors = color_vals
        ) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

}

#' Split a Plot List into even Chunks for Paged Plotting
#'
#' @param plot_list A list of ggplots
#' @param plots_per_page Integer value. Number of plots that should be plotted per page.
#' @param dummy_panel Optional default plot that should be used as filler on last
#' page if there are empty spaces. Defaults to a completely blank placeholder.
#' @param ... Additional arguments passed to `patchwork::wrap_plots()` to format the plots on
#' a page as desired
#' @return Plots each 'page' of plots as a separate plot. Intended for use inside of `pdf()`, for example.
plot_list_to_pages <- function(plot_list, 
                               plots_per_page, 
                               dummy_plot = ggplot() + theme_void(),
                               ...){
    n_plots <- length(plot_list)
    n_pgs <- ceiling(n_plots/plots_per_page)
    for(ipg in 1:n_pgs){
        istart <- (ipg-1)*plots_per_page + 1
        iend <- min(ipg*plots_per_page, n_plots)
        plt_subset <- plot_list[istart:iend]
        plt_count <- length(plt_subset)
        if(plt_count < plots_per_page){
            while(plt_count < plots_per_page){
                plt_subset <- c(plt_subset, list(dummy_plot))
                plt_count <- plt_count + 1
            }
        }
        print(patchwork::wrap_plots(plt_subset, ...))
    }
}

#' Save a plotlist to html giving unique dimensions for each plot
#'
#' @param plot_list A list of plots
#' @param output_file Full path to output html file
#' @param heights Numeric. A vector of heights for each plot. Must be same length as plot_list. Default NULL will use RMD default = 5
#' @param widths Numeric. A vector of widths for each plot. Must be same length as plot_list. Default NULL will use RMD default = 7
#' @param titles Character. A vector of titles for each plot. Must be same length as plot_list. Default NULL will use no titles. Titles are
#' recommended to make the output file searchable for easier navigation
plotlist_to_html <- function(plot_list, output_file, heights = NULL, widths = NULL, titles = NULL, in_rmd = "plotlist_to_html.Rmd",...){
  rmarkdown::render(input = in_rmd,
                    params = list(plot_list = plot_list, heights = heights, widths = widths, titles = titles),
                    output_file = output_file,
                   ...)
}
# plotlist_to_html(plot_list,"test.html", c(5,7,10), c(7,12, 20), c("A","B","C"))
