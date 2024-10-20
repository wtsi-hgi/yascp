#!/usr/bin/env Rscript

# __date__ <- '2020-06-26'
# __version__ <- '0.0.1'

library(data.table)

# Function to parse command line arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  arg_list <- list()

  for (i in seq(1, length(args), by = 2)) {
    arg_name <- sub("--", "", args[i])
    arg_value <- args[i + 1]
    arg_list[[arg_name]] <- arg_value
  }

  # if ("version" %in% names(arg_list)) {
  #   cat("Version", __version__, "\n")
  #   quit(status = 0)
  # }

  if (!all(c("az_file", "mapping", "out_file") %in% names(arg_list))) {
    stop("Missing required argument(s). Required: --az_file, --mapping, --out_file", call. = FALSE)
  }

  return(arg_list)
}

# Parse arguments
args <- parse_args()

# Read data with error handling for compressed files

tryCatch({
  All_Data <- fread(args$az_file, sep = "\t", showProgress = FALSE)
  compressed <- TRUE
}, error = function(e) {
  All_Data <- fread(args$az_file, sep = "\t", showProgress = FALSE)
  compressed <- FALSE
})
# print(All_Data)
# print('loaded')
# Read mapping file
Mappings <- fread(args$mapping, sep = "\t", data.table = FALSE)
# print(Mappings)
# Data manipulation
D1 <- All_Data
D1$idx1 <- seq_len(nrow(D1))
setDT(D1)

# Perform the replacement
# We will create a new column for each of the columns in Mappings
print(D1)
for (col in colnames(Mappings)) {
  if (col == "L2") {
    next  # Skip this iteration if col is L2
  }
  D1[, paste0(col, '_predicted.celltype.l2') := Mappings[[col]][match(predicted.celltype.l2, Mappings$L2)]]
}
print(D1)
# Write to output file
fwrite(D1, file = args$out_file, sep = "\t", quote = FALSE)

cat("Done\n")
