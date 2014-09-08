#' @include utils.R
#' @importFrom GenomicAlignments readGAlignmentsFromBam sequenceLayer cigar 
#' @importFrom Biostrings consensusMatrix start 
#' @importFrom GenomicRanges seqlengths seqnames
#' @importFrom Rsamtools ScanBamParam
#' @importFrom IRanges mcols
NULL

#' Generate a consensus table between a BAM alignment and a reference sequence.
#' 
#' @details
#' \code{consensus_table} generates a \code{\link{data.table}} with 5 columns:
#' The \emph{Position} along the reference, the \emph{Nucleotide} (A, C, G, T, or Gap [-]),
#' the number of reads at a Position for each Nucleotide, the 
#' 
#' @param bamfile Path to the BAM file.
#' @param reffile Path to the reference sequence (Fasta format).
#' @return A data.table. See \code{Details}.
#' @export
consensus_table <- function(bamfile, reffile) {
  assert_that(is.readable(bamfile), is.readable(reffile))
  ## make sure to import the SEQ field from the BAM file
  param <- ScanBamParam(what = "seq")
  ## other fields omitted by default are the query quality (QUAL) field,
  ## and the mapping quality (MAPQ) fields.
  aln <- readGAlignmentsFromBam(bamfile, index = bamfile, param = param)
  refname <- levels(seqnames(aln))
  qseq <- mcols(aln)$seq
  rseq <- read_fasta(reffile, as = "DNAStringSet")
  qseq_on_ref <- sequenceLayer(qseq, cigar(aln))
  
  ## start() returns the 1-based leftmost position of the clipped query. It should
  ## be zero-based (see package rbamtools), so we have to subtract 1
  qmat <- consensusMatrix(qseq_on_ref, as.prob = FALSE, shift = start(aln) - 1L, 
                          width = seqlengths(aln), baseOnly = TRUE)
  rmat <- consensusMatrix(rseq, baseOnly = TRUE)
  assert_that(all(dim(qmat) == dim(rmat)))
  
  ## report probabilities
  prob.qmat <- sweep(qmat, 2, colSums(qmat), `/`)

  ## create data.tables
  qtbl <- as.data.table(qmat)
  setnames(qtbl, names(qtbl), as.character(seq_len(ncol(qtbl))))
  qtbl[, Nucleotide := c("A", "C", "G", "T", "-")]
  qtbl <- melt(qtbl, id.vars = "Nucleotide", value.name = "Nreads", variable.name = "Position")
  setkeyv(qtbl, c("Position", "Nucleotide"))
  
  prob.qtbl <- as.data.table(prob.qmat)
  setnames(prob.qtbl, names(prob.qtbl), as.character(seq_len(ncol(prob.qtbl))))
  prob.qtbl[, Nucleotide := c("A", "C", "G", "T", "-")]
  prob.qtbl <- melt(prob.qtbl, id.vars = "Nucleotide", value.name = "Probability.query", variable.name = "Position")
  setkeyv(prob.qtbl, c("Position", "Nucleotide"))
  
  rtbl <- as.data.table(rmat)
  setnames(rtbl, names(rtbl), as.character(seq_len(ncol(rtbl))))
  rtbl[, Nucleotide := c("A", "C", "G", "T", "-")]
  rtbl <- melt(rtbl, id.vars = "Nucleotide", value.name = "Probability.ref", variable.name = "Position")
  setkeyv(rtbl, c("Position", "Nucleotide"))
  
  ## join
  res <- qtbl[prob.qtbl][rtbl]
  attr(res, "refname") <- refname
  res
}


plot_quality_line <- function(prob.tbl, match.mismatch.tbl, smooth = TRUE,
                              annotate = TRUE, annot.size = 6, refname = "") {
  p <- ggplot(prob.tbl, aes(x = as.numeric(Position), y = Probability.query)) +
    geom_line(alpha = 0.75, size = 0.25, colour = "gray30") +    
    xlab("Nucleotide position") +
    ylab("Proportion query reads matching reference") +
    ggtitle(bquote(atop("Alignment Quality", atop(italic(.(refname)), "")))) +
    theme_bw()
  if (annotate) {
    p <- p + 
      geom_point(data = match.mismatch.tbl, aes(x = as.numeric(Position), y = Probability.query),
                            shape = 16, colour = "#f5410b", size = 2) +
      geom_text(data = match.mismatch.tbl, aes(y = Probability.query - .02, label = Mismatch),
                       size = annot.size)
  }
  if (smooth) {
    p <- p + 
      stat_smooth(method = "gam", formula = y ~ s(x))
  }
  p
}


plot_quality_box <- function(prob.tbl) {
  p <- ggplot(prob.tbl, aes(x = 1, y = Probability.query)) + 
    geom_boxplot() +
    scale_color_continuous(breaks = NULL) +
    ylab("Proportion query reads matching reference") +
    theme(axis.title.x = element_blank()) +
    theme_bw()
  p
}


plot_depth_line <- function(prob.tbl, refname = "") {
  p <- ggplot(prob.tbl, aes(x = as.numeric(Position), y = total_reads)) +
    geom_line(alpha = 0.75, size = 0.25, colour = "gray30") +
    xlab("Nucleotide position") +
    ylab("Alignment Depth") +
    ggtitle(bquote(atop("Alignment Depth (Number of reads/bp)", atop(italic(.(refname)), "")))) +
    theme_bw()
  p
}


plot_depth_box <- function(prob.tbl) {
  p <- ggplot(prob.tbl, aes(x = 1, y = total_reads)) + 
    geom_boxplot() +
    scale_color_continuous(breaks = NULL) +
    ylab("Alignment Depth") +
    theme(axis.title.x = element_blank()) +
    theme_bw()
  p
}


#' Plot alignment quality.
#' 
#' This function graphs the probability to recover the reference sequence from
#' an multiple alignment of nanopore reads.
#' 
#' @param x A \code{\link{consensus_table}}.
#' @param annotate Annotate mismatch positions.
#' @param annot.size If \code{annotate = TRUE}, size of the annotation text. 
#' @param smooth Plot a smoothing function (gam: \code{y ~ s(x)}).
#' @return Plots the probability to recover the reference sequence and summarises it in a
#' boxplot.
#' @export
plot_alignment_quality <- function(x, smooth = TRUE, annotate = TRUE, annot.size = 4) {
  assert_that(is(x, "data.table"), has_attr(x, "refname"))
  refname <- paste("Refname:", attr(x, "refname"))
  ## extract the proportion of query reads matching the reference sequence for
  ## each base position.
  prob.tbl <- x[Probability.ref == 1L, list(Position, Nucleotide, Probability.query)]
  ## aggregate the total number of reads for each positions.
  total_reads <- x[, list(total_reads = sum(Nreads)), by = "Position"]
  prob.tbl <- prob.tbl[total_reads]
  ## A Match-Mismatch table comparing the "majority rule" consensus sequence
  ## to the reference
  mismatch.tbl <- x[, .SD[which.max(Probability.query)], by = "Position"][Probability.ref == 0]
  setkey(mismatch.tbl, "Position")
  match.mismatch.tbl <- 
    x[Probability.ref == 1][mismatch.tbl][, list(
      Position,
      Probability.ref = Probability.query,
      Probability.query = Probability.query.1,
      Mismatch = paste(Nucleotide, Nucleotide.1, sep = "/")
    )]
  ## Plots
  pdl <- plot_depth_line(prob.tbl, refname)
  pdb <- plot_depth_box(prob.tbl)
  pql <- plot_quality_line(prob.tbl, match.mismatch.tbl, smooth,
                           annotate, annot.size, refname)
  pqb <- plot_quality_box(prob.tbl)
  ## Multiplot
  lay <- matrix(c(1,1,1,2,3,3,3,4), nrow = 2, byrow = TRUE)
  multiplot(pdl, pdb, pql, pqb, layout = lay)
  
  invisible(list(pql = pql, pdl = pdl, match.mismatch.tbl = match.mismatch.tbl))
}
