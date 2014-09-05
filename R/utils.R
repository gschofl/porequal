#' @importFrom assertthat assert_that is.readable has_attr
#' @importFrom rbamtools cigarData
#' @importFrom Biostrings readDNAStringSet
NULL

trim <- function (x, trim = '\\s+') {
  assert_that(is.vector(x))
  gsub(paste0("^", trim, "|", trim, "$"), '', x)
}

#' Return the CIGAR string from a \code{bamAlign} object.
#' 
#' @param align A \code{\linkS4class{bamAlign}} object
#' @return A CIGAR string.
#' @export
cigar_string <- function(align) {
  assert_that(is(align, "bamAlign"))
  cig <- cigarData(align)
  paste0(paste0(cig[, "Length"], cig[, "Type"]), collapse = "")
}

#' Read a fasta file.
#'
#' @param path Path to fasta file.
#' @param as Return sequence as a character vector or as a
#' \code{\linkS4class{DNAStringSet}} object.
#' @return A named character vector or a \code{\linkS4class{DNAStringSet}} object.
#' @export
read_fasta <- function(path, as = c("character", "DNAStringSet")) {
  as <- match.arg(as, c("character", "DNAStringSet"))
  if (as == "DNAStringSet") {
    return(readDNAStringSet(path))
  } 
  con <- file(path, open = "r")
  on.exit(close(con))
  lines <- readLines(con)
  headers <- trim(sub(">", "", lines[grepl("^>", lines)]))
  sequences <- trim(lines[grepl("^[^>]", lines)])
  assert_that(length(headers) == length(sequences))
  names(sequences) <- headers
  sequences
}

#' Multiple plot function
#' 
#' @details
#' ggplot objects can be passed in \code{...}, or to \code{plotlist} (as a list of ggplot objects)
#' - cols:   Number of columns in layout
#' - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#'
#' If the layout is something like \code{matrix(c(1,2,3,3), nrow=2, byrow=TRUE)},
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#' 
#' @param ... ggplot objects.
#' @param plotlist A list of ggplot objects.
#' @param cols Number of columns in layout.
#' @param layout  A matrix specifying the layout. If present, 'cols' is ignored.
#' @author Winston Chang
#' @references \url{http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/}
#' @export
multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots <- length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }  
  if (numPlots == 1) {
    print(plots[[1]]) 
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

