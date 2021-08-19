
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