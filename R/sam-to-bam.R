#' @include utils.R
#' @importFrom gsmisc has_command replace_ext
NULL

#' Convert a samfile to a bamfile
#' 
#' @param samfile Path to the sam file.
#' @param reffile Path to the reference sequence (Fasta format)
#' @return A list containing the paths to a bam file (.bam) and
#'   an index file (.bai)
#' @export
sam2bam <- function(samfile, reffile) {
  assert_that(file.exists(samfile))
  assert_that(file.exists(reffile))
  assert_that(has_command("samtools"))
  bamfile <- replace_ext(samfile, replacement = "bam")
  cmd <- paste("samtools faidx", reffile)
  system(cmd, intern = TRUE)
  faifile <- paste0(reffile, ".fai")
  assert_that(file.exists(faifile))
  cmd <- paste("samtools view -bt", faifile, samfile, ">", bamfile)
  system(cmd, intern = TRUE)
  assert_that(file.exists(bamfile))
  ## generate index file
  baifile <- paste0(bamfile, ".bai")
  cmd <- paste("samtools index", bamfile)
  system(cmd, intern = TRUE)
  assert_that(file.exists(baifile))
  invisible(list(bam = bamfile, bai = baifile))
}
