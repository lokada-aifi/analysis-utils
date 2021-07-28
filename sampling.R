#' Randomly Split Observations into Groups with Stratification
#'
#' Randomly splits a set of observations into groups. If a dataframe of 
#' factor-like metadata is supplied, observations will be balanced across 
#' stratification groups. Number of observations will be balanced across 
#' groups by default unless group proportions are supplied. 
#'
#' @param obs Vector of unique observations, ie SampleID's.
#' @param n_groups Integer. The number of groups to split the data into.
#' @param grp_proportions Vector of proportions. Default (NULL) assumes groups
#' have even numbers of observations. If values supplied, must be length equal to 
#' the number of groups and add up to 1.
#' @param strat_levels Data frame. Default (NULL) assumes there are no statification
#' levels. If supplied, should be an observation x metadata data frame with rows equal
#' to the number of observations and one column per stratification variable. Row order
#' must correspond to observation order.
#' @param seed_use An integer. Randomization seed for reproducibility. Default is 3.
#' @verbose Logical value. If TRUE (default) will print the final observation groups and 
#' tables of observations by stratification variables per group, if applicable.
#' @return A list of length n_groups, with all values of obs assigned to one list element.
#' @examples
#' @export
#' my_df <- data.frame(Sample = paste0("Sample", 1:100),
#'                     Cohort = c(rep("CohortA", 30), rep("CohortB", 70)),
#'                     Diagnosis = c(rep("Healthy",10), rep("Disease",20),
#'                                   rep("Healthy",20), rep("Disease",50))
#'                                          
#' my_obs <- my_df$Sample
#' my_variables <- my_df[,c("Cohort","Diagnosis")]
#' my_grp_prop <- c(0.7, 0.15, 0.15)
#' grp_assignments <- randomize_samples(obs = my_obs, 
#'                                 n_groups = 3, 
#'                                 grp_proportions=my_grp_prop, 
#'                                 strat_levels = my_variables,
#'                                 seed_use = 1)
#' grp_assignments
randomize_samples <- function(obs, 
                              n_groups,
                              grp_proportions = NULL,
                              strat_levels = NULL, 
                              seed_use = 3, 
                              verbose = TRUE){
  
  # Check input
  assertthat::assert_that(all(!duplicated(obs)),
                          msg = sprintf("Duplicate values detected in obs ['%s']. Please supply a unique vector",
                                        paste0(obs[duplicated(obs)], collapse = "', ''")
                          ))
  if(!is.null(strat_levels)){
    assertthat::assert_that(nrow(strat_levels) == length(obs),
                            msg = "Number of rows in strat_levels must be equal to length of obs")
    assertthat::assert_that("data.frame" %in% class(strat_levels),
                            msg = "strat_levels should be a data.frame")                        
  }
  
  n_obs <- length(obs)
  group_names <- paste0("Group", 1:n_groups)
  cat(sprintf("seed is %s", seed_use), sep = "/n")
  set.seed(seed_use)
  
  # Determine group proportions
  if(!is.null(grp_proportions)){
    assertthat::assert_that(length(grp_proportions) == n_groups,
                            msg = "grp_proportions should be length n_groups.")
    assertthat::assert_that(sum(grp_proportions) >= 0.99 & sum(grp_proportions) <= 1.01,
                            msg = "grp_proportions should add up to 1")
    balance_total_grp_size <- FALSE
  } else {
    # Proportions default to be equal across groups
    grp_proportions <- rep(1/n_groups, n_groups)
    balance_total_grp_size <- TRUE
  }
  
  # sampling workflow for stratified and non-stratified    
  if(!is.null(strat_levels)){  # stratified workflow
    # Get all stratification groups
    level_combos <- interaction(strat_levels)  # only works for data.frame, not matrix
    unique_levels <- unique(level_combos)
    
    # Run the non-stratified sample randomization workflow per stratification group
    # strata_list is a list of length equal to stratification levels
    # with each list item being a list of length n_groups containing sampled observations
    strata_list <- lapply(unique_levels, function(x){
      i_temp <- which(level_combos == x)
      obs_temp <- obs[i_temp]
      grps_temp <- randomize_samples(obs=obs_temp, 
                                     n_groups=n_groups, 
                                     grp_proportions=grp_proportions, 
                                     strat_levels=NULL,
                                     seed_use=seed_use + runif(1,1,1000), 
                                     verbose = FALSE)
      names(grps_temp) <- x
      print(x)
      grps_temp
    })
    
    # Combine per-strata groups to make final groups, aggregating one group per strata.
    # For equal group proportions, select combinations of strata groups 
    # that will balance total N per group. For nonequal grouos, just combine strata based 
    # on their existing order.
    if(balance_total_grp_size){
      obs_grouped_strat <- strata_list[[1]]  # start with first strata
      for(i in 2:length(strata_list)){ # for each additional strata, merge in each strata group to one of the final groups
        # get next strata list, order groups by decreasing number of samples
        strata_grp <- strata_list[[i]]
        lengths_new <- sapply(strata_grp, length)
        order_new <- order(lengths_new, decreasing = TRUE) 
        strata_grp <- strata_grp[order_new]
        
        # get current merged list, order groups by increasing number of samples
        lengths_current <- sapply(obs_grouped_strat, length)
        order_current <- order(lengths_current)
        obs_grouped_strat <- obs_grouped_strat[order_current]
        
        # merge the merged list index-wise with the ordered strata list
        obs_grouped_strat <- Map(c, obs_grouped_strat, strata_grp)
      }
    } else {
      obs_grouped_strat <- Reduce(function(x,y){Map(c, x, y)}, strata_list)
    }
    
    names(obs_grouped_strat) <- group_names[1:length(obs_grouped_strat)]
    
    # Output to console 
    if(verbose){
      print(obs_grouped_strat)
      sapply(obs_grouped_strat, length)
      tab_list <- lapply(obs_grouped_strat, function(x){
        temp_df <- cbind(obscol = x, strat_levels[match(x, obs), ])
        temp_df <- temp_df %>%
          dplyr::group_by(across(all_of(colnames(strat_levels)))) %>%
          dplyr::summarize(N = n(), .groups = "drop")
        temp_df
      })
      print(tab_list)
    }
    
    return(obs_grouped_strat)
    
  } else {  # non-stratified workflow
    # Determine the number of samples per group
    grp_size <- round(n_obs*grp_proportions)  # general starting point
    extra_n <- n_obs - sum(grp_size)  # see if any counts should be adjusted
    if(abs(extra_n) > 0){  
      n_correction <- c(rep(sign(extra_n)*1, abs(extra_n)), 
                        rep(0, n_groups-abs(extra_n)))
      
      # randomizes which group will get the sample number adjustment
      # n_correction <- sample(n_correction, length(n_correction), replace = FALSE)  
      
      # favors adding extra samples to groups with fewer samples or removing samples from groups with more
      if(extra_n > 0){
        n_correction <- n_correction[rank(grp_size, ties.method = "random")]
      } else {
        n_correction <- n_correction[rank(grp_size*-1, ties.method = "random")]
      }
      grp_size <- grp_size + n_correction
    }
    
    # Create a sequence of group assignments based on n observations per group
    grp_seq <- sapply(1:n_groups, function(i){
      rep(i, times = grp_size[i])
    })
    grp_seq <- unlist(grp_seq)
    
    # Randomly shuffle the observations, then split into list by the group sequence 
    obs_random <- sample(obs, n_obs, replace = FALSE)
    obs_grouped <- split(obs_random, grp_seq)
    names(obs_grouped) <- group_names[unique(grp_seq)]
    
    # Add dummy list elements if there were no samples for a given group
    if(length(obs_grouped) < n_groups){
      n_null <- n_groups - length(obs_grouped)
      dummy_list = rep(list(NULL),n_null)
      names(dummy_list) <- setdiff(group_names, names(obs_grouped))
      obs_grouped <- c(obs_grouped, dummy_list)
      obs_grouped <- obs_grouped[order(names(obs_grouped))]
    } 
    
    # Output to console 
    if(verbose){
      print(obs_grouped)
      sapply(obs_grouped, length)
    }
    
    return(obs_grouped)
    
  }
}