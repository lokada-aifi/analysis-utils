#' Plot Enrichment as in # Same plot as in fgsea::plotEnrichment() but with flexible formatting options
#'
#' @param pathway Gene set to plot. Inherited from fgsea::plotEnrichment() 
#' @param stats Gene-level statistics (weights that were used for pathway enrichment) Inherited from fgsea::plotEnrichment().
#' @param gseaParam GSEA parameter. Inherited from fgsea::plotEnrichment().
#' @param ticksSize Width of vertical line corresponding to a gene (default: 0.2)
#' @param enrichmentcolor Character value, color for the main enrichment plot NES trace, default "green"
#' @param base_size Positive numeric value, global text size parameter
#' @param theme_ls List of named arguments passed to pass to ggplot2::theme()
#' @return A ggplot2 object
plotEnrichment_fmt <- function(pathway, stats, gseaParam = 1, ticksSize = 0.2, enrichmentcolor = "green", base_size = 10, theme_ls = NULL){
    default_theme_ls <- list(
        panel.border = element_blank(), 
        panel.grid.minor = element_blank()
    )
    unspec_param <- setdiff(names(default_theme_ls), names(theme_ls))
    if(length(unspec_param) > 0){
        theme_ls <- c(theme_ls, default_theme_ls[unspec_param])
    }
    
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
        returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    diff <- (max(tops) - min(bottoms))/8
    x = y = NULL
    g <- ggplot(toPlot, aes(x = x, y = y)) + 
        geom_point(color = enrichmentcolor, size = 0.1) +
        geom_hline(yintercept = max(tops), 
                   colour = "red", 
                   linetype = "dashed") + 
        geom_hline(yintercept = min(bottoms), 
                   colour = "red", 
                   linetype = "dashed") + 
        geom_hline(yintercept = 0, colour = "black") + 
        geom_line(color = enrichmentcolor) + 
        theme_bw(base_size = base_size) + 
        geom_segment(data = data.frame(x = pathway), 
                     mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = ticksSize) + 
        do.call(theme, theme_ls) +
        labs(x = "rank", y = "enrichment score")
    g
} 

#' Plot Pathway Enrichment across Comparison Groups
#'
#' Plot enrichment for specified pathways for each set of group gene ranks (ie different celltypes or populations)
#' where each group is a different subplot in the pathway plot
#'
#' @param pathways Character vector of pathway names to plot
#' @param pathway_list List of genes in each pathway, list item names are the same as pathway names
#' @param rank_list List of gene rankings for each analysis group. List item names should be a descriptive group name. Per-gene
#' stats/ranking values should be the names of the associated genes.
#' @param gseaParam Numeric value passed to plotEnrichment(), default 1.
#' @param ticksSize Numeric value for plot tick size passed to plotEnrichment(), default 0.2.
#' @param ... Additional arguments passed to patchwork::wrap_plots() for formatting group panels in each pathway plot
#' @return A list of plots, one per pathway. Each plot is a patchwork panel

plot_enrichment_groups <- function(pathways, pathway_list, rank_list, gseaParam = 1, ticksSize = 0.2, ...){
    group_names <- names(rank_list)
    
    # Iterate over pathways
    pl_list <- lapply(pathways, function(pwy){
        # Get plots for each group
        pw_grp_list <- mapply(SIMPLIFY = FALSE,function(rnk, grp_nm){
            fgsea::plotEnrichment(pathway_list[[pwy]], rnk, gseaParam = gseaParam, ticksSize = ticksSize) +
            ggtitle(pwy, subtitle = grp_nm)
        }, rank_list, group_names)
        
        # Format group plots
        patchwork::wrap_plots(pw_grp_list, ...)
    })   
    
    return(pl_list)
}

#' Plot Pathway Enrichment across Comparison Groups, Coloring by Significance
#'
#' Plot enrichment for specified pathways for each set of group gene ranks (ie different celltypes or populations)
#' where each group is a different subplot in the pathway plot. Color by significance
#' based on gsea results
#'
#' @param gsea_res GSEA result dataframe output of fgsea(), potentially with other metadata columns added.
#' @param pathway_list List of genes in each pathway, list item names are the same as pathway names
#' @param rank_list List of gene rankings for each analysis group. List item names should be a descriptive group name. Per-gene
#' stats/ranking values should be the names of the associated genes.
#' @param group_col Character value. Name of column in gsea_res containing the group names that match the `rank_list` names
#' @param pw_col Character value. Name of column in gsea_res containing the pathway names that match the `pathway_list` names
#' @param pval_col Character value. Name of column in gsea_res containing the pvalue or FDR values
#' @param alpha Numeric value between 0 and 1, default 0.05. Cutoff for p-value/FDR significance.
#' @param gseaParam Numeric value passed to plotEnrichment(), default 1.
#' @param ticksSize Numeric value for plot tick size passed to plotEnrichment(), default 0.2.
#' @param signif_col Character vector of length 2 containing color values for non-significant and significant pathway 
#' enrichment, respectively. Default "blue" and "orange".
#' @param base_size Positive numeric value, global text size parameter
#' @param theme_ls List of named arguments passed to pass to ggplot2::theme()
#' @param ... Additional arguments passed to patchwork::wrap_plots() for formatting group panels in each pathway plot
#' @return A list of plots, one per pathway. Each plot is a patchwork panel if multiple groups were plotted per pathway.

plot_enrichment_groups_signif <- function(gsea_res,
                                          pathway_list, 
                                          rank_list,  
                                          group_col, 
                                          pw_col = "pathway",
                                          pval_col = "padj", 
                                          alpha = 0.05,
                                          gseaParam = 1, 
                                          ticksSize = 0.2,
                                          signif_col = c(`ns` = "blue", `signif`= "orange"),
                                          base_size = 12,
                                          theme_ls = NULL,
                                          ...){
    assertthat::assert_that(all(names(rank_list) %in% unique(gsea_res[[group_col]])))
    group_names <- names(rank_list)
    
    signif_gs <- gsea_res %>%
        filter(get(pval_col) <= alpha) %>%
        pull(all_of(pw_col)) %>% 
        unique()
    
    if(length(signif_gs) > 0){
        pw_plotlist <- lapply(signif_gs, function(pwy){
            pw_grp_list <- mapply(function(rnk, grp_nm){
                p_res <- gsea_res %>%
                    filter(get(pw_col) == pwy, get(group_col) == grp_nm) %>%
                    pull(all_of(pval_col)) 
                is_signif <- p_res <= alpha
                col <- signif_col[[is_signif + 1]]
                plotEnrichment_fmt(pathway = pathway_list[[pwy]], 
                                   stats = rnk, 
                                   gseaParam = gseaParam, 
                                   ticksSize = ticksSize, 
                                   enrichmentcolor = col, 
                                   base_size = base_size, 
                                   theme_ls= theme_ls) +
                    ggtitle(pwy, subtitle = grp_nm)
            }, rank_list, group_names, SIMPLIFY = FALSE)
            pw_panel <- patchwork::wrap_plots(pw_grp_list, ...)
            return(pw_panel)
        })    
        return(pw_plotlist)
    } else {
        cat("No Significant Genesets")
    }

}


#' Get scaled values for custom breaks
#'
#' Get scaled values for custom breaks in a given value vector for use
#' defining breaks in scale_gradientn(), for example
#'
breaks_to_scaledbreaks <- function(breaks, x){
    rescaled_weights <- scales::rescale(x)
    rescaled_breaks <- quantile(rescaled_weights, probs = ecdf(x)(breaks))
    return(rescaled_breaks)
}


#' DotPlot of Pathway Enrichement
#' 
plot_signif_pathways <- function(gsea_res, 
                                 p_col = "padj", 
                                 p_cutoff = 0.05,
                                 x_group = "contrast", 
                                 y_group = "gs", 
                                 split_y_group = NULL,
                                 color_col = "NES",
                                 color_pal = NULL,
                                 color_name = NULL,
                                 base_size = 14){
    gsea_res <- gsea_res %>%
        filter(get(p_col) <= p_cutoff)
    
    g <- ggplot(gsea_res, aes(!!as.name(x_group), !!as.name(y_group))) +
               geom_point(aes(color = !!as.name(color_col), 
                              size = -1*log10(!!as.name(p_col)))) +
        theme_bw(base_size = base_size)
    
    if(!is.null(split_y_group)){
        g <- g +
            facet_grid(paste0("~", split_y_group), scales = "free_x", space = "free")
    }
    
    if(!is.null(color_pal)){
        scaled_breaks <- breaks_to_scaledbreaks(color_pal, gsea_res[[color_col]])
        color_name <- ifelse(is.null(color_name), color_col, color_name)
        g <- g + scale_color_gradientn(colors = names(color_pal), values = scaled_breaks, na.value = "gray",  name = color_name)
    }
    
    return(g)
}