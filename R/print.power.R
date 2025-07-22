print.power <- function(x, digits = 3, latex.output = FALSE, template = 1, summaries = TRUE, ...) {

    class(x) <- paste("power", template, sep = "")
    
    print(x, digits, latex.output, summaries, ...)
    
    
}



