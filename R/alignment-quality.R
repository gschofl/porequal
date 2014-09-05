#' @include utils.R
#' @importFrom GenomicAlignments readGAlignmentsFromBam sequenceLayer cigar 
#' @importFrom Biostrings consensusMatrix start 
#' @importFrom GenomicRanges seqlengths seqnames
#' @importFrom Rsamtools ScanBamParam
#' @importFrom IRanges mcols
NULL

#' Generate a consensus table between a BAM alignment and a reference sequence.
#' 
#' @param bamfile Path to the BAM file.
#' @param reffile Path to the reference sequence (Fasta format).
#' @param offset Using \code{\link[GenomicAlignments]{sequenceLayer}}, the refseq and
#' the query reads seem to be one off. Use this offset to cut off bases from the 
#' beginning of the query reads and the end of the refseq.
#' @return A data.table
#' @export
consensus_table <- function(bamfile, reffile, offset = 1) {
  assert_that(is.readable(bamfile), is.readable(reffile))
  param <- ScanBamParam(what = "seq")
  aln <- readGAlignmentsFromBam(bamfile, param = param)
  refname <- levels(seqnames(aln))
  
  qseq <- mcols(aln)$seq
  rseq <- read_fasta(reffile, "DNAStringSet")
  qseq_on_ref <- sequenceLayer(qseq, cigar(aln))
  qmat <- consensusMatrix(qseq_on_ref, as.prob = FALSE, shift = start(aln), 
                          width = seqlengths(aln))  
  ## UGLY HACK
  ## The refseq and the query reads seem to be one off
  ## So I take out the first column from the queries and the last column from refseq (offset).
  if (offset > 0) {
    qmat <- qmat[c(1:4, 16), -(1:offset)]
    rmat <- consensusMatrix(rseq)
    rmat <- rmat[c(1:4, 16), -(ncol(rmat) + 1 - 1:offset)]
  }
  assert_that(all(dim(qmat) == dim(rmat)))

  ## report probabilities
  prob.qmat <- sweep(qmat, 2, colSums(qmat), `/`)
  
  ## create data.tables
  qtbl <- as.data.table(qmat)
  setnames(qtbl, names(qtbl), as.character(seq_len(ncol(qtbl))))
  qtbl[, Nucleotide := rownames(qmat)]
  qtbl <- melt(qtbl, id.vars = "Nucleotide", value.name = "Nreads", variable.name = "Position")
  setkeyv(qtbl, c("Position", "Nucleotide"))
  
  prob.qtbl <- as.data.table(prob.qmat)
  setnames(prob.qtbl, names(prob.qtbl), as.character(seq_len(ncol(prob.qtbl))))
  prob.qtbl[, Nucleotide := rownames(qmat)]
  prob.qtbl <- melt(prob.qtbl, id.vars = "Nucleotide", value.name = "Probability.query", variable.name = "Position")
  setkeyv(prob.qtbl, c("Position", "Nucleotide"))
  
  rtbl <- as.data.table(rmat)
  setnames(rtbl, names(rtbl), as.character(seq_len(ncol(rtbl))))
  rtbl[, Nucleotide := rownames(rmat)]
  rtbl <- melt(rtbl, id.vars = "Nucleotide", value.name = "Probability.ref", variable.name = "Position")
  setkeyv(rtbl, c("Position", "Nucleotide"))
  
  ## join
  res <- qtbl[prob.qtbl][rtbl]
  attr(res, "refname") <- refname
  res
}

#' Plot alignment quality.
#' 
#' This function graphs the probability to recover the reference sequence from
#' an multiple alignment of nanopore reads.
#' 
#' @param x A \code{\link{consensus_table}}.
#' @param smooth Plot a smoothing function (gam: \code{y ~ s(x)}).
#' @return Plots the probability to recover the reference sequence and summarises it in a
#' boxplot.
#' @export
plot_alignment_quality <- function(x, smooth = TRUE) {
  assert_that(is(x, "data.table"), has_attr(x, "refname"))
  refname <- paste("Refname:", attr(x, "refname"))
  prob.tbl <- x[Probability.ref == 1, list(Position, Nucleotide, Probability.query)]
  nreads <- x[, list(Nreads = sum(Nreads)), by = "Position"]
  prob.tbl <- prob.tbl[nreads]
  
  qual.line <- ggplot(prob.tbl, aes(x = as.numeric(Position), y = Probability.query)) +
    geom_line(alpha = 0.75, size = 0.25, colour = "gray30") +
    xlab("Nucleotide position") +
    ylab("Fraction of queries matching reference") +
    ggtitle(bquote(atop("Alignment Quality", atop(italic(.(refname)), "")))) +
    theme_bw()  
  if (smooth) {
    qual.line <- qual.line + stat_smooth(method = "gam", formula = y ~ s(x))
  }
  
  qual.box <- ggplot(prob.tbl, aes(x = 1, y = Probability.query)) + 
    geom_boxplot() +
    scale_color_continuous(breaks = NULL) +
    ylab("Fraction of queries matching reference") +
    theme(axis.title.x = element_blank()) +
    theme_bw()
  
  depth.line <- ggplot(prob.tbl, aes(x = as.numeric(Position), y = Nreads)) +
    geom_line(alpha = 0.75, size = 0.25, colour = "gray30") +
    xlab("Nucleotide position") +
    ylab("Alignment Depth") +
    ggtitle(bquote(atop("Alignment Depth (Number of reads/bp)", atop(italic(.(refname)), "")))) +
    theme_bw()
  
  depth.box <- ggplot(prob.tbl, aes(x = 1, y = Nreads)) + 
    geom_boxplot() +
    scale_color_continuous(breaks = NULL) +
    ylab("Alignment Depth") +
    theme(axis.title.x = element_blank()) +
    theme_bw()  
  
  lay <- matrix(c(1,1,1,2,3,3,3,4), nrow = 2, byrow = TRUE)
  multiplot(depth.line, depth.box, qual.line, qual.box, layout = lay)
  
  invisible(NULL)
}
