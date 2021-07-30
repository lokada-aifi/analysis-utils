# Set output plot dimensions
plotdim <- function(w, h){
    options(repr.plot.width = w, repr.plot.height = h)
}


#' Make a logfile function
#' 
#' Creates a custom function that can be used to log comments/text to a 
#' logfile within an analysis. log is output in format:
#' [timestamp]text with the bracketed timestamp and text separated
#' by delimter of choice. The timestamp timezone can be set when
#' creating the custom fucntion.
#'
#' @param outfile A file path (character value). Path of the output
#' logfile to create
#' @param delim Character value. A delimiter to use to separate tiime and text
#' @param tzone A lubridate timezone value. If supplied will format the 
#' time for your time zone of choice. Default is null, will use default Sys.time() o
#' output
#' @return A function to create a log file entry. The function has the following arguments:
#' * text (the text to log)
#' * append Logical value, default TRUE. Whether or not the log comment should be appended to
#'   the existing file or overwrite it.
makelogfun <- function(outfile, delim = " ", tzone = NULL){
    newfun <- function(text, 
                       bappend = TRUE){
        text = as.character(text)
        if(!is.null(tzone)){
            time_text <- sprintf("[%s]%s%s", 
                             lubridate::with_tz(Sys.time(), tzone = tzone),
                             delim, 
                             text)
        } else {
            time_text <- sprintf("[%s]%s%s", Sys.time(), delim, text)
        }
        if(!file.exists(outfile) | !bappend){
            write(time_text, file = outfile, append = FALSE)
        } else {
            write(time_text, file = outfile, append = bappend)
        }
    }
    return(newfun)
}