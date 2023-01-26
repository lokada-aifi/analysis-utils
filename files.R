#' Move files matching pattern from one directory to another
#' @param from String. From directory path
#' @param to String. To directory path
#' @param pattern. String value. String or regex for pattern matching files in 
#' 'from' to move to 'to'
#' @param recursive Logical, default TRUE. Whether or not to search for files in
#' 'from' recursively
#' @param verbose Logical, default TRUE. Whether or not to print out progress/
#' status messages
#' @param overwrite Logical, default FALSE. Whether files in 'to' are allowed to 
#' be overwritten if an attempt is made to copy into a file that already exists
move_files <- function(from, to, pattern, recursive = TRUE, verbose = TRUE, overwrite = FALSE){
    # create output directory
    if(!dir.exists(to)){
        if(verbose){
            message(sprintf("Creating output directory '%s'", to))
        }
        dir.create(to, recursive = TRUE)
    }
    
    # identify downloaded files
    all_files <- list.files(path = from, pattern = pattern, recursive = recursive, full.names = TRUE)
    if (verbose){
        message(sprintf("%s files found for pattern [%s]: \n\t%s", 
                        length(all_files), pattern, paste(all_files, collapse = "\n\t")))
    }
    
    if(length(all_files) > 0){
        for (fname in all_files){
            to_file <- file.path(to, basename(fname))
            if (file.exists(to_file) & overwrite == FALSE){
                message(sprintf("Destination file %s already exists. File not copied. Resolve copies or set overwrite = TRUE.", to_file))
            } else {
                file.rename(from = fname, to = to_file)
                message(sprintf("Moved file %s", to_file))

                if(length(list.files(dirname(from), recursive = TRUE)) == 0){
                    system(paste("rm -r", dirname(from)))
                    message(sprintf("Removed empty parent directory %s", dirname(from)))
                }
            }

        }
    }
    return("Finished")
} 