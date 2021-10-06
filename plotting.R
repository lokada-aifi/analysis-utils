# Helper for converting NULL values in a list to na values (ie to preserve all original indices during downstream operations)
nullToNA <- function(x) {
            x[sapply(x, is.null)] <- NA
            return(x)
}

# Fix the facet spacing when a ggplot is converted to a plotly object for a single axis
# Makes all facets equally sized and spaced along the specified axis
fix_plotly_facet_axis <- function(ply_plot, axis_fix = "x", padding = 0.05){
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
        texty <- sapply(testplot$x$layout$annotations,"[[", "y") 
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

# wrapper to correct all plotly facet spacing along both x and y axis.
fix_plotly_facets <- function(ply_plot, padding_x = 0.05, padding_y = 0.05){
    ply_plot <- fix_plotly_facet_axis(ply_plot, axis_fix = "x", padding = padding_x)
    ply_plot <- fix_plotly_facet_axis(ply_plot, axis_fix = "y", padding = padding_y)

    return(ply_plot)
}