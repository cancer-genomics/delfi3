


if (!require(data.table)){
  install.packages('data.table')
  library(data.table)
}

#setDTthreads(percent = 100)

if (!exists("bedops_path")) stop("Bedops is not available!")

unstarch_path <- paste0(bedops_path, "/unstarch")

#' Check starch
#'
#' Check if file is in BEDOPS starch format
#'
#' @param file Path to starch file
#'
#' @return Boolean (T/F)
#' @export
#'
#' @examples
#' check_starch(file = "/path/to/file.starch")
check_starch <- function(file) {
  if (missing(file)) {
    stop("No file argument provided")
  }
  # Check if file exists
  if (!file.exists(file)) {
    stop(paste("File does not exist:", file))
  }

  # Check file is in starch format
  is_starch <-
    system(paste(unstarch_path, " --is-starch", file, "2>/dev/null"),
           intern = T)
  if (is_starch == "0") {
    return(FALSE)
  } else{
    return(TRUE)
  }
}


#' List chromosomes
#'
#' Lists all chromosomes in BEDOPS starch file
#'
#' @param file Path to starch file
#'
#' @return Character vector of chromosome names
#' @export
#'
#' @examples
#' list_chromosomes(file = "/path/to/file.starch")
list_chromosomes <- function(file) {
  # Check that starch file exists and is properly formatted
  if (!check_starch(file)) {
    stop(paste("File is not in starch format:", file))
  }

  chromosomes <-
    system(paste(unstarch_path, "--list-chromosomes", file),
           intern = T)
  return(chromosomes)
}


#' Starch length
#'
#' Returns the number of rows in the starch file (very fast).
#'
#' @param file Path to starch file
#' @param chromosome Chromosome to extract (very fast).
#'
#' @return Integer
#' @export
#'
#' @examples
#' starch_length(file = "/path/to/file.starch")
#' starch_length(file = "/path/to/file.starch", chromosome = "chrX")
starch_length <- function(file, chromosome) {
  # Check that starch file exists and is properly formatted
  if (!check_starch(file)) {
    stop(paste("File is not in starch format:", file))
  }


  # Unstarch file
  if (missing(chromosome)) {
    n_elements <-
      system(paste(unstarch_path, "--elements", file),
             intern = T)
  } else{
    # If chromosome argument has been provided, check if chromosome exists within file
    if (!chromosome %in% list_chromosomes(file)) {
      stop(paste(chromosome, "not in", file))
    }

    n_elements <-
      system(paste(unstarch_path, as.character(chromosome), "--elements", file),
             intern = T)
  }


  return(as.integer(n_elements))
}


#' Unstarch
#'
#' Decompress file in BEDOPS starch format. Wrapper for data.table::fread(). Check data.table documentation for additional arguments.
#'
#' @param file Path to starch file
#' @param chromosome Chromosome to extract (very fast).
#' @param ...
#'
#' @return data.table in bed format
#' @export
#'
#' @examples
#' unstarch(file = "/path/to/file.starch")
#' unstarch(file = "/path/to/file.starch", chromosome = "chrX")
unstarch <- function(file, chromosome, tmpdir = Sys.getenv("MYSCRATCH"), ...) {
  ### Check that starch file exists and is properly formatted
  if (!check_starch(file)) {
    stop(paste("File is not in starch format:", file))
  }

  ### Unstarch file

  # Create tmpdir if it does not already exist
  if(is.null(tmpdir) | tmpdir==""){
    tmpdir <- "tmp_bedops"
  }
  dir.create(tmpdir,showWarnings = F, recursive = T)

  if (missing(chromosome)) {
    bed <- fread(cmd = paste(unstarch_path, file), sep = "\t",
                 fill = F, tmpdir = tmpdir,
                 ...)

    if(nrow(bed)!=starch_length(file)){
      stop(paste("Unstarch failed:", file))
    }
  } else{
    ### If chromosome argument has been provided, check if chromosome exists within file
    if (!chromosome %in% list_chromosomes(file)) {
      stop(paste(chromosome, "not in", file))
    }

    bed <-
      fread(cmd = paste(unstarch_path, as.character(chromosome), file), sep = "\t",
            fill = F, tmpdir = tmpdir,
            ...)

    if(nrow(bed)!=starch_length(file,as.character(chromosome))){
      stop(paste("Unstarch failed:", fragfile))
    }
  }

  # Format header
  # if("#chrom" %in% colnames(bed)){
  #   setnames(bed, "#chrom", "chrom")
  # }
  colnames(bed)[1:3]<-c("chrom","chromStart","chromEnd")

  # Return data.table in BED format
  return(bed)
}

#' Unstarch to GRanges
#'
#' Decompress file in BEDOPS starch format, convert to GRanges. Wrapper for data.table::fread() and GenomicRanges::makeGRangesFromDataFrame(). Check data.table documentation for additional arguments.
#'
#' @param file Path to starch file
#' @param chromosome Chromosome to extract (very fast).
#' @param ...
#'
#' @return GRanges object
#' @export
#'
#' @examples
#' unstarch_to_granges(file = "/path/to/file.starch")
#' unstarch_to_granges(file = "/path/to/file.starch", chromosome = "chrX")
unstarch_to_granges <- function(file, chromosome, tmpdir = Sys.getenv("MYSCRATCH"), ...){
  if(! "GenomicRanges" %in% installed.packages()){
    print("Installing package GenomicRanges")
    if (!require("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
    }
    BiocManager::install("GenomicRanges")
  }
  bed <- unstarch(file, chromosome, tmpdir, ...)
  gr <- GenomicRanges::makeGRangesFromDataFrame(df = bed,
                                                starts.in.df.are.0based = T, keep.extra.columns = T)
  return(gr)
}
