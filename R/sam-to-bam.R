#' @include utils.R
#' @importFrom gsmisc has_command replace_ext
NULL

#' Convert a SAM file to a BAM file.
#' 
#' @param samfile Path to the SAM file.
#' @param reffile Path to the reference sequence (Fasta format)
#' @return The paths to the BAM file.
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
  bamfile
}
