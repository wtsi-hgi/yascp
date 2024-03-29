#!/usr/bin/env Rscript
library(Matrix)

args = commandArgs(trailingOnly=TRUE)
cs_in_dir = args[1] #'samplename'
out_dir =args[2] # '../donor_ids.tsv'
filtered_vcf_file = args[3] 

print('cs_in_dir')
print(cs_in_dir)

print('out_dir')
print(out_dir)

print('filtered_vcf_file')
print(filtered_vcf_file)
#' Update three sparse matrices
#'
#' Update three sparse matrices `AD`, `DP`, and `OTH` given a VCF file storing
#' a list of SNPs that passed filtering.
#'
#' @param in_dir The input cellSNP dir.
#' @param out_dir The output dir, should not be the same with `in_dir`.
#' @param filtered_vcf_file The VCF file storing SNPs that passed filtering.
#' @param is_vcf_gzipped A bool. Whether the `base.vcf` in `in_dir` is gzipped.
#' @return Void.
#' @examples

update_cellsnp_matrices <- function(
  in_dir,
  out_dir,
  filtered_vcf_file,
  is_vcf_gzipped = TRUE)
{
  if (in_dir == out_dir)
    stop("Error: output dir should not be the input dir!")

  if (! dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)

  suffix <- ifelse(is_vcf_gzipped, ".gz", "")
  raw_vcf_file <- paste0(in_dir, "/cellSNP.base.vcf", suffix)
  assert_e(raw_vcf_file)

  for (mtx_name in c("AD", "DP", "OTH")) {
    mtx_file <- sprintf("%s/cellSNP.tag.%s.mtx", in_dir, mtx_name)
    assert_e(mtx_file)

    print(sprintf("Update %s matrix ...", mtx_name))
    update_matrix(
      matrix_file = mtx_file,
      new_matrix_file = sprintf("%s/cellSNP.tag.%s.mtx", out_dir, mtx_name),
      raw_vcf_file = raw_vcf_file,
      filtered_vcf_file = filtered_vcf_file 
    )
  }
}


assert_e <- function(path) {
  if (is.null(path) || (! file.exists(path)))
    stop(sprintf("Error: '%s' does not exist!", path))
}


update_matrix <- function(
  matrix_file,
  new_matrix_file,
  raw_vcf_file,
  filtered_vcf_file,
  verbose = FALSE
)
{
  mtx <- Matrix::readMM(file = matrix_file)

  raw_vcf <- read.delim(raw_vcf_file, header = FALSE, comment.char = "#",
                        stringsAsFactors = FALSE)
  if (verbose)
    str(raw_vcf)

  raw_features <- paste0(raw_vcf$V1, "_", raw_vcf$V2)     # chrom_pos
  rownames(mtx) <- raw_features
  
  flt_vcf <- read.delim(filtered_vcf_file, header = FALSE, comment.char = "#",
                        stringsAsFactors = FALSE)
  if (verbose)
    str(flt_vcf)

  flt_features <- paste0(flt_vcf$V1, "_", flt_vcf$V2)     # chrom_pos
  flt_mtx <- mtx[raw_features %in% flt_features, ]

  Matrix::writeMM(flt_mtx, new_matrix_file) 
}


update_cellsnp_matrices(
   in_dir = cs_in_dir,
   out_dir = out_dir,
   filtered_vcf_file = filtered_vcf_file
 )
