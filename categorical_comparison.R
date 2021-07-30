
#' Summarize proportion of categorical variables within another categorical variable
#'
#' Creates a dataframe with calculted percentages of categorical variable levels within
#' another categorical variable
#'  
#' @param 
#' 
#'
#'
#'
prop_summary <- function(df, bycol, sumcol){
  df <- df %>%
    dplyr::group_by(!!as.name(bycol), !!as.name(sumcol)) %>%
    dplyr::summarize(N = n(), .groups = "drop") %>%
    dplyr::group_by(!!as.name(bycol)) %>%
    dplyr::mutate(N_all = sum(N), 
                  Pct = N/sum(N)*100) %>%
    dplyr::ungroup()
  df
}
