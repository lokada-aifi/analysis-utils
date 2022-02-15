
#' Summarize proportion of categorical variables within another categorical variable
#'
#' Creates a dataframe with calculted percentages of categorical variable levels within
#' another categorical variable
#'  
#' @param df A data frame
#' @param groupby_col Character vector. One or more column names from df.
#' @param sumcol Character vector. One or more column names from df.
#'
prop_summary <- function(df, groupby_col, sumcol){
  df <- df %>%
    dplyr::group_by(across(all_of(c(bycol, sumcol)))) %>%
    dplyr::summarize(N = n(), .groups = "drop") %>%
    dplyr::group_by(across(all_of(sumcol))) %>%
    dplyr::mutate(N_all_group = sum(N), 
                  Pct = N/sum(N)*100) %>%
    dplyr::ungroup()
  df
}


#' Convert table to matrix
#'
#' Convert output of table() to a matrix
#' @param tab A table 
#' @return A matrix of the input table
#' examples
#' mydf <- data.frame()
tab_to_mat <- function(tab){
    mat <- matrix(tab, ncol = ncol(tab))
    colnames(mat) <- colnames(tab)
    rownames(mat) <- rownames(tab)
    mat
} 


#' Cleanup htest Output
#' 
#' Converts an htest object (such as outputs of wilcox.test, fisher.test, etc) to data frame output
#'
#' @param htest_result A object of class "htest"
#' @param remove_unused_cols Logical value. If TRUE (default), removes the columns that are all NA's from
#' the output data frame. Some htest objects do not have all of these values.
#' @return A one row data.frame containing the htest object values
htest_list_to_df <- function(htest_result, remove_unused_cols = TRUE){
    df <- data.frame(
        data.name = htest_result$data.name,
        class = class(htest_result),
        method = htest_result$method,
        statistic.name = ifelse(!is.null(names(htest_result$statistic)),names(htest_result$statistic),NA),
        statistic = as.numeric(ifelse(!is.null(htest_result$statistic), htest_result$statistic, NA)),
        parameter = as.numeric(ifelse(!is.null(htest_result$parameter), htest_result$parameter,NA)),
        p.value = as.numeric(htest_result$p.value),
        estimate.name = ifelse(!is.null(htest_result$estimate),names(htest_result$estimate), NA),
        estimate = as.numeric(ifelse(!is.null(htest_result$estimate),htest_result$estimate,NA)),
        lower.conf.limit = as.numeric(ifelse(!is.null(htest_result$conf.int), htest_result$conf.int[[1]], NA)),
        upper.conf.limit = as.numeric(ifelse(!is.null(htest_result$conf.int), htest_result$conf.int[[2]], NA)),
        conf.level = as.numeric(ifelse(!is.null(htest_result$conf.int), attributes(htest_result$conf.int)$conf.level, NA)),
        null.value = as.numeric(htest_result$null.value),
        alternative = htest_result$alternative
    )
    
    if(remove_unused_cols){
        iNA <- apply(df, 2, function(x){all(is.na(x))})
        if(sum(iNA) > 0){
            df <- df[, -which(iNA)]
        }
    }
    
    df
}

#' Perform parallelized Wilcoxon test on Long Data
#'
#' Perform a parallelized wilcoxon test. The data will
#' be subsetted by values in supplied columns and tests
#' performed on each subset.
#'
#' @param dat Data frame with one row per assay/result per samples
#' @param subset_cols Categorical variable columns within which to perform the wilcoxon test 
#' @param result_col Numeric result value that will be used for Wilcoxon teset
#' @param contrast_col Categorical binary variable defining the two groups being compared in the 
#' Wilcoxon test
#' @param ncores Integer value. The number of cores to use for parallel processsing. Default NULL
#' @param remove_unused_col Logical value, default TRUE. Passed to \\code{htest_list_to_df()}
#' @param pval_adj_method Character value passed to \\code{p.adjust()}. Default "BY"
do_wilcox <- function(dat, 
                      subset_cols = c("Assay_2", "Day"), 
                      result_col = "NPX_Change", 
                      contrast_col = "Cohort", 
                      ncores = NULL,
                      remove_unused_cols = TRUE,
                      pval_adj_method = "BY"){
    if(is.null(ncores) || ncores < 1){
        ncores <- max(1, parallel::detectCores()-2)
        message(sprintf("N cores not defined. Setting number of cores to %s", ncores))
    }
    
    dat_list <- dat %>%
        unite("wilcox_subset", all_of(subset_cols), remove = FALSE, sep = "__") %>%
        group_by(across(all_of(subset_cols))) %>%
        group_split() 
            
    wilcox_list <- parallel::mclapply(dat_list, function(x){
         wilcox.test(as.formula(sprintf("%s ~ %s", result_col, contrast_col)), data = x) %>%
             htest_list_to_df() %>%
             mutate(group = unique(x$wilcox_subset)) %>%
             relocate(group, data.name) %>%
             separate(group, into = subset_cols, remove = FALSE, sep = "__")
    }, mc.cores = ncores,mc.cleanup = TRUE)
    
    medians_df <- dat %>%
        group_by(across(all_of(c(subset_cols, contrast_col)))) %>%
        summarize(Median = median(get(result_col), na.rm = TRUE),
                  N = sum(!is.na(get(contrast_col))), .groups = "drop") %>%
        pivot_wider(names_from = all_of(contrast_col), values_from = c(N, Median),names_sep = "_") 
        
    wilcox_df <- bind_rows(wilcox_list) %>%
        mutate(p.adj = p.adjust(p.value, method = pval_adj_method)) %>%
        mutate(p.adj.method = pval_adj_method) %>%
        left_join(medians_df, by = subset_cols) %>%
        arrange(p.value, p.adj)
    
    return(wilcox_df)
}