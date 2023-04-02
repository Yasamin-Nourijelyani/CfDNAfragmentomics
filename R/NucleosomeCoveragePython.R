#Add citation to Jonathan

#' Plot coverage from pythpn file
#' sample_bed_file <- system.file("extdata", "p1.bed", package = "CfDNAfragmentomics")
#' sample_bed <- read.table(sample_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
#'
#' library(readr)
#' TFBS_loc <- read_csv("inst/extdata/TFBS_loc.csv")
#'
#' pip install tqdm
#'
#' @export
nucleosomeCoverage1 <- function(bed_file, tfbs, window, name) {

  library(reticulate)

  #import from Jonathan's coverage plot
  mymodule <- source_python("./R/nucleosome_occupancy.py")

  use_python("myenv", required=FALSE)
  conda_create("myenv", python = "3.10")


  coverage = get_coverage(bed_file)

  aggregated_coverage = aggregate_coverage(coverage, tfbs, window)

  plot_coverage(aggregated_coverage, name.split('.')[1])


}




