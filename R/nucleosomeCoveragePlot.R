# Purpose: Plotting of cfDNA coverage data of patient data
# Author: Yasamin Nouri Jelyani
# Date: 2022-12-05
# Version: 0.1.0
# Bugs and Issues: Does not produce plots because chromosome length is too long to put on the x-axis

#------------------------Exported Function------------------------


#' Plot of the coverage for cfDNA fragments
#'
#' A visualization function that generates a coverage plot showing
#' nucleosome nucleosome patient data to visually
#' detect TFBS.
#'
#'
#' @param sample_bed A dataframe for the a sample data
#' (unknown if cancerous or not individual)
#' with at least 3 columns: first column is a string representing the chromosome
#'  number for the cfDNA for instance: "chr1",
#'  the second column are integer values indicating the start location of
#'  the read and third column are integer values indicating end location of
#'  the read.
#' These are the patient who we want to check if they have cancer or not.
#' The fragment lengths
#' of patients will be compared to these healthy data similar to what is done in
#' the paper by Katsman et. al (1).
#' These samples are compared to the controls using a t-test.
#'
#' @param chromosome_length The integer length of the chromosome being analysed
#' Lengths can also be referenced from http://www.insilicase.com/Web/Chromlen.aspx
#' given the chromosome number.
#'
#'
#'
#' @return Returns a plot of cfDNA coverage
#'
#' @references
#' 1- Katsman, E., Orlanski, S., Martignano, F., Fox-Fisher, I., Shemer,
#'  R., Dor, Y., Zick, A.,
#' Eden, A., Petrini, I., Conticello, S. G., & Berman, B. P. (2022).
#' Detecting cell-of-origin and cancer-specific methylation features of
#' cell-free DNA from Nanopore sequencing.
#' Genome biology, 23(1), 158.
#' https://doi.org.myaccess.library.utoronto.ca/10.1186/s13059-022-02710-1
#'
#' 2- Illumina, https://support.illumina.com/downloads/nextera-flex-for-enrichment-BED-files.html
#'
#' 3- “The T-Test.” JMP,
#' https://www.jmp.com/en_ca/statistics-knowledge-portal/t-test.html#:~
#' :text=t%2DTest%20assumptions&amp;text=The%20data%20are%20continuous.,
#' The%20distribution%20is%20approximately%20normal.
#'
#'@examples
#' \dontrun{
#'
#' # Example 1:
#'
#' sample_bed_file <- system.file("extdata", "d1.bed", package = "CfDNAfragmentomics")
#' sample_bed <- read.table(sample_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
#'
#' ## basic usage
#'
#' library(ggplot2)
#' plot <- nucleosomeCoveragePlot(sample_bed = sample_bed, 237995947)
#'
#'
#'
#'
#'}
#'
#'
#' @export
#' @import ggplot2
#' @importFrom gridExtra grid.arrange

nucleosomeCoveragePlot <- function(coverage_df) {
  # coverage_df is made from aggregation in aggregate_coverage function

  # plot histogram of cfDNA coverage

  #library(ggplot2)
  plot_coverage<-ggplot(coverage_df, aes(x = nucleotide, y = coverage)) +
    geom_point() +
    ggtitle("Patients")



}

# [END]
