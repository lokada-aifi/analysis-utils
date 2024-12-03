# Helper for converting NULL values in a list to na values (ie to preserve all original indices during downstream operations)

# Fix plotly facets
nullToNA <- function(x) {
            x[sapply(x, is.null)] <- NA
            return(x)
}

equalize_facets <- function(ply_plot, axis_fix = "x", padding = 0.05){
    # calculate facet spacing based on # of unique axises present
    axis_pattern <- paste0(axis_fix, "axis.*")
    axis_names <- grep(axis_pattern, names(ply_plot$x$layout), value = T)
    orig_domains <- lapply(ply_plot$x$layout[axis_names],"[[", "domain")
    orig_domains <- unique(orig_domains)
    orig_domains <- orig_domains[order(sapply(orig_domains,"[[",1))]
    n_domains <- length(orig_domains)
    domain_span <- (1-padding*(n_domains-1))/n_domains
    min_vals <- cumsum(c(0,rep(padding + domain_span, n_domains-1)))
    ranges <- lapply(min_vals, function(x){c(x, x + domain_span)})
    
    # fix the main plotting field facet domains
    for(i in seq_along(axis_names)){
        axis_nm <- axis_names[i]
        orig_domain_axis <- ply_plot$x$layout[[axis_nm]]$domain
        ichange <- which(sapply(orig_domains, function(x){all(x==orig_domain_axis)}))
        ply_plot$x$layout[[axis_nm]]$domain <- ranges[[ichange]]
    }
    
    # fix the facets strips
    var1 <- paste0(axis_fix, 0)
    var2 <- paste0(axis_fix, 1)
    minorig <- sapply(ply_plot$x$layout$shapes,"[[",var1)
    maxorig <- sapply(ply_plot$x$layout$shapes,"[[",var2)
    anchororig <- as.numeric(nullToNA(sapply(ply_plot$x$layout$shapes,"[[","yanchor")))
    for(i in seq_along(orig_domains)){
        oldrange <- orig_domains[[i]] 
        newrange <- ranges[[i]]
        
        ichange <- which(minorig == oldrange[1] & maxorig == oldrange[2])
        ply_plot$x$layout$shapes[ichange] <- lapply(ply_plot$x$layout$shapes[ichange], function(x){
            x[[var1]] <- newrange[1]
            x[[var2]] <- newrange[2]
            x
        }) 
        
        # fix y facet strip anchor
        if(axis_fix == "y"){
            oldanchor <- oldrange[2]
            newanchor <- newrange[2]
            ichange <- which(anchororig == oldanchor)
            
            ply_plot$x$layout$shapes[ichange] <- lapply(ply_plot$x$layout$shapes[ichange], function(x){
                    x$yanchor <- newanchor
                    x
            }) 
        }
    }
    
    
    # Fix the facet text
    if(axis_fix == "x"){
        textx <- sapply(ply_plot$x$layout$annotations,"[[", "x") 
        orig_medians <- sapply(orig_domains, median)
        new_medians <- sapply(ranges, median)
        for(i in seq_along(orig_medians)){
            oldmed <- orig_medians[i] 
            newmed <- new_medians[i]
            ichange <- which(textx == oldmed)

            ply_plot$x$layout$annotations[ichange] <- lapply(ply_plot$x$layout$annotations[ichange], function(x){
                x$x <- newmed
                x
            }) 
        }
    }
    if(axis_fix == "y"){
        texty <- sapply(ply_plot$x$layout$annotations,"[[", "y") 
        for(i in seq_along(orig_domains)){
            oldmax <- orig_domains[[i]][2]
            newmax <- ranges[[i]][2]
            ichange <- which(texty == oldmax)

            ply_plot$x$layout$annotations[ichange] <- lapply(ply_plot$x$layout$annotations[ichange], function(x){
                x$y <- newmax
                x
            }) 

        }
    }
    
    
    return(ply_plot)
}

fix_plotly_facets <- function(ply_plot, padding_x = 0.05, padding_y = 0.05){
    ply_plot <- equalize_facets(ply_plot, axis_fix = "x", padding = padding_x)
    ply_plot <- equalize_facets(ply_plot, axis_fix = "y", padding = padding_y)

    return(ply_plot)
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
#' @param title Character. Single string value to use as HTML title. Default ""

plotlist_to_html <- function(plot_list, output_file, heights = NULL, widths = NULL, titles = NULL, title = '', in_rmd = "plotlist_to_html.Rmd",...){
  rmarkdown::render(input = in_rmd,
                    params = list(plot_list = plot_list, heights = heights, widths = widths, titles = titles, title=title),
                    output_file = output_file,
                   ...)
}
# plotlist_to_html(plot_list,"test.html", c(5,7,10), c(7,12, 20), c("A","B","C"))