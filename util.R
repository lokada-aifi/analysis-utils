#' Remove items in a delimited string
#'
#' Remove specific levels from a delimited string. Example application is to
#' select specific directory levels from a file path.
#'
#' @param stringval A character value or vector of values.
#' @param level Single integer, vector of integers, or list of integer vectors of length equal to number of files. 
#' The delimited  levels to remove. 
#' @param delim Default "/". The delimiter to split the string
#' @param keep_last Integer value default NULL. If provided will override level and keep the last N items in the delimited string
#' @return Vector of modified stringval values with specified items removed
remove_delim_level <- function(stringval, level = 1, delim = "/", keep_last = NULL){
    path_ls <- strsplit(stringval, split = delim)
    if(!is.null(keep_last)){
        path_lengths <- sapply(path_ls, length)
        level <- lapply(path_lengths, function(n){1:(n - keep_last)})
    } 
    if(!is.list(level) & length(level)!= 1 | length(level)!= length(path_ls)){
        message(sprintf("Length levels to remove not equal to 1 or length files. Assuming the same multiple items (%s) should be removed from all strings",
               paste(level, collapse = ",")))
        level = list(level)

    }
    
    path_ls <- mapply(function(stringval, lvl){
        stringval[-c(lvl)]
    }, path_ls, level, SIMPLIFY = FALSE)
    
    path_ls <- sapply(path_ls, function(x){paste(x, collapse = delim)})

    return(path_ls)
}

#' Summarize loaded packages in table for portability
#' @return A data.frame of a summary of packages described in \code{sessionInfo()}
summarize_packages <- function(){
    types <- c("basePkgs", "otherPkgs", 'loadedOnly')
    out_df <- purrr::map_df(types,.f= function(x){
        info <- sessionInfo()[[x]]
        if(is.list(info)){
            data.frame(package = names(info), type = x, version = sapply(info, '[[','Version'))
        } else {
            data.frame(package = info, type = x, version = NA)
        }
    })
    row.names(out_df) <- NULL
    out_df
    
}