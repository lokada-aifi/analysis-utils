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
#' the number of groups
#' @param strat_levels Data frame. Default (NULL) assumes there are no statification
#' levels. If supplied, should be an observation x metadata data frame with rows equal
#' to the number of observations and one column per stratification variable. Row order
#' must correspond to observation order.
#' @param seed_use An integer. Randomization seed for reproducibility. Default is 3.
#' @return A list of length n_groups, with all values of obs assigned to one list element.
#' @examples
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
#' sapply(grp_assignments, length)
#' lapply(grp_assignments, function(x){
#'    with(my_df %>% filter(Sample %in% x), 
#'         table(Cohort, Diagnosis))
#' })
randomize_samples <- function(obs, 
                              n_groups,
                              grp_proportions = NULL,
                              strat_levels = NULL, 
                              seed_use = 3, 
                              verbose = TRUE){
    
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
    set.seed(seed_use)
    
    if(!is.null(grp_proportions)){
        assertthat::assert_that(length(grp_proportions) == n_groups,
                                msg = "grp_proportions should be length n_groups.")
        assertthat::assert_that(sum(grp_proportions) >= 0.99 & sum(grp_proportions) <= 1.01,
                                msg = "grp_proportions should add up to 1")
        balance_total_grp_size <- FALSE
    } else {
        grp_proportions <- rep(1/n_groups, n_groups)
        balance_total_grp_size <- TRUE
    }
    
    if(!is.null(strat_levels)){ 
        level_combos <- interaction(strat_levels)  # only works for data.frame, not matrix
        unique_levels <- unique(level_combos)
        
        sampled_levels <- lapply(unique_levels, function(x){
            i_temp <- which(level_combos == x)
            obs_temp <- obs[i_temp]
            grps_temp <- randomize_samples(obs=obs_temp, 
                                      n_groups=n_groups, 
                                      grp_proportions=grp_proportions, 
                                      strat_levels=NULL,
                                      seed_use=seed_use, 
                                      verbose = FALSE)
            names(grps_temp) <- x
            cat(x)
            grps_temp
        })
        
        if(balance_total_grp_size){
            obs_grouped <- sampled_levels[[1]]
            for(i in 2:length(sampled_levels)){
                temp <- sampled_levels[[i]]
                lengths_current <- sapply(obs_grouped, length)
                order_current <- order(lengths_current)
                lengths_new <- sapply(temp, length)
                order_new <- order(lengths_new, decreasing = TRUE)
                
                # Reorder groups by number of samples
                obs_grouped <- obs_grouped[order_current]
                temp <- temp[order_new]
                
                obs_grouped <- Map(c, obs_grouped, sampled_levels[[i]])
            }
        } else {
            obs_grouped <- Reduce(function(x,y){Map(c, x, y)}, sampled_levels)
        }
        
        names(obs_grouped) <- paste0("Group", 1:n_groups)
        
        if(verbose){
            print(obs_grouped)
            sapply(obs_grouped, length)
            tab_list <- lapply(obs_grouped, function(x){
                temp_df <- cbind(obscol = x, strat_levels[match(x, obs), ])
                temp_df <- temp_df %>%
                    dplyr::group_by(across(all_of(colnames(strat_levels)))) %>%
                    dplyr::summarize(N = n(), .groups = "drop")
               temp_df
            })
            print(tab_list)
        }
        
        return(obs_grouped)
        
    } else {
        grp_size <- round(n_obs*grp_proportions)
        extra_n <- n_obs - sum(grp_size)
        
        if(extra_n > 0){  # the group that gets more n will be random.
            n_correction <- c(rep(sign(extra_n)*1, extra_n), 
                              rep(0, n_groups-extra_n))
            n_correction <- sample(n_correction, length(n_correction), replace = FALSE)
            grp_size <- grp_size + n_correction
        }
        
        grp_seq <- sapply(1:n_groups, function(i){
            rep(i, each = grp_size[i])
        })
        grp_seq <- unlist(grp_seq)
        
        obs_random <- sample(obs, n_obs, replace = FALSE)
        
        obs_grouped <- split(obs_random, grp_seq)
        names(obs_grouped) <- paste0("Group", 1:n_groups)
        
        if(verbose){
            print(obs_grouped)
            sapply(obs_grouped, length)
        }
    
        return(obs_grouped)
        
    }
}

