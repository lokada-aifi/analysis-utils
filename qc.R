# Add qc flags to dataset based on existing numeric variables and specific thresholds
# Currently a little slow--can speed up
# Examples:
# testdf <- data.frame(pct.mito = c(0,20,49,50,51,80,NA,10,10),
#                     n_genes = c(NA, 50,0,199,200,201, 500,1000,1000),
#                     n_umis = c(299,300,301, 5000, 24999,25000,25001,10000,NA))
# flaglist <- data.frame(variable = c("pct.mito", "n_genes","n_umis", "n_umis"),
#                        direction_pass = c("lt", "gte", "gte", "lt"),
#                        threshold = c(50, 200, 500, 25000))
# flaglist2 <- data.frame(variable = c("pct.mito", "n_genes","n_umis", "n_umis","adfas"),
#                        direction_pass = c("lt", "gte", "gte", "lt", "gt"),
#                        threshold = c(50, 200, 500, 25000, 300))
# flaglist3 <- data.frame(variable = c("pct.mito", "n_genes","n_umis", "n_umis","n_genes"),
#                        direction_pass = c("lt", "gte", "gte", "lt", "gtes"),
#                        threshold = c(50, 200, 500, 25000, 300))
# test1 <- add_lowquality_flags(df = testdf, flaglist = flaglist)
# test1
# test2 <- add_lowquality_flags(df = testdf, flaglist = flaglist2)
# test2
# test3 <- add_lowquality_flags(df = testdf, flaglist = flaglist3)
# test3
add_qc_flags <- function(df, flaglist, newcolprefix = "QC_", add_summary_col = TRUE){
    good_ops <- c("lt", "gt", "lte", "gte", "eq")
    assertthat::assert_that(all(flaglist$direction_pass %in% good_ops),
                           msg = sprintf("All direction_pass values must be one of ['%s']. Needs correction: ['%s']. Should describe passing values relative to the threshold.",
                           paste(good_ops, collapse = "','"), 
                           paste(unique(setdiff(flaglist$direction_pass, good_ops)), collapse = "','")))
    newvarnames <- character()

    # loop through each flag
    for(i in 1:nrow(flaglist)){
        flag_vals <- flaglist[i, ]
        if(!flag_vals$variable %in% colnames(df)){
            warning(sprintf("Flag variable %s not found in input data.frame. Skipping flag.", flag_vals$variable))
            next()
        }
        compfun <- switch(flag_vals$direction_pass,
                         lt = `<`,
                         gt = `>`,
                         lte = `<=`,
                         gte = `>=`,
                         eq = `==`)
        newvar <- sprintf("%s%s_%s%s", newcolprefix,flag_vals$variable, flag_vals$direction_pass, flag_vals$threshold)
        newvarnames <- c(newvarnames, newvar)
        flag <- compfun(df[,flag_vals$variable], flag_vals$threshold)
        flag <- factor(flag, levels = c(FALSE, TRUE), labels = c("FAIL","PASS"))
        df[,newvar] <- flag
    }
    
    if(add_summary_col){
        summaryflag <- apply(df[,newvarnames, drop = FALSE], 1, function(x){ifelse(any(is.na(x)), FALSE, all(x == "PASS"))})
        summaryflag <- factor(summaryflag, levels = c(FALSE, TRUE), labels = c("FAIL","PASS"))
        summary_colname <- sprintf("%sFlagSummary", newcolprefix)
        assertthat::assert_that(length(summaryflag) == nrow(df))

        df[,summary_colname] <- summaryflag
    }
    
    return(df)
}