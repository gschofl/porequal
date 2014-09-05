#' @include utils.R
#' @importFrom rbamtools bamReader bamRange bamClose alignDepth getPos getDepth
NULL

globalVariables("Align.Depth")

#' Plot alignment depth.
#' 
#' @param bamfile Path to the bam file.
#' @param reffile Path to the reference sequence (Fasta format).
#' @return Side effects. Graphs alignment depth as a line-plot and summarises the alignment depth in a
#' boxplot.
#' @export
alignment_depth <- function(bamfile, reffile) {
  reader <- bamReader(bamfile, idx = TRUE)
  reflen <- nchar(read_fasta(reffile))
  coords <- c(0, 0, reflen)
  range <- bamRange(reader, coords)
  bamClose(reader)
  ad <- alignDepth(range)
  df <- data.frame(Position = getPos(ad), Align.Depth = getDepth(ad))
  refname <- paste("Refname:", ad@refname)
  
  adp <- ggplot(df, aes(x = Position, y = Align.Depth)) +
    geom_line(size = 0.25, colour = "grey50") +
    ggtitle(bquote(atop("Alignment Depth", atop(italic(.(refname)), "")))) +
    ylab("Alignment Depth") +
    theme_bw()
  
  adbox <- ggplot(df, aes(x = 1, y = Align.Depth)) + geom_boxplot() +
    scale_color_continuous(breaks = NULL) +
    ggtitle(bquote(atop("Alignment Depth", atop(italic(.(refname)), "")))) +
    ylab("Alignment Depth") +
    theme(axis.title.x = element_blank()) +
    theme_bw()
  
  lay <- matrix(c(1,1,1,2), nrow = 1, byrow = TRUE)
  multiplot(adp, adbox, layout = lay)
  
  invisible(NULL)
}
