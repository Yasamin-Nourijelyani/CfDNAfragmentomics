# Purpose: Plotting of cfDNA fragment lengths density of patient data as
# compared with control data
# Author: Yasamin Nouri Jelyani
# Date: 2022-12-05
# Version: 0.1.0
# Bugs and Issues:

#------------------------Exported Function------------------------


#' Plot of the density for nucleosome fragment lengths
#'
#' A visualization function that generates a density plot showing
#' nucleosome fragment lengths for control and patient data to visually
#' compare the mono-nucleosome and di-nucleosome fragment length densities.
#'
#' @param controls_bed A dataframe for the a control data (healthy individual)
#' with at least 3 columns: first column is a string representing the chromosome
#'  number for the cfDNA for instance: "chr1",
#'  the second column are integer values indicating the start location of
#'  the read and third column are integer values indicating end location of
#'  the read.
#' These are the healthy individuals who do not have cancer.The fragment lengths
#' of patients will be compared to these healthy data similar to what is done in
#' the paper by Katsman et. al (1).
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
#' @return Returns a plot of cfDNA fragment size densities for control
#' and patient data.
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
#' controls_bed_file <- system.file("extdata", "d2.bed", package = "CfDNAfragmentomics")
#' controls_bed <- read.table(controls_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
#'
#' sample_bed_file <- system.file("extdata", "p1.bed", package = "CfDNAfragmentomics")
#' sample_bed <- read.table(sample_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
#'
#' ## basic usage
#'
#' library(ggplot2)
#' plot <- nucleosomeDensityPlot(controls_bed = controls_bed, sample_bed = sample_bed)
#'
#'
#' # Example 2:
#'
#' controls_bed_file <- system.file("extdata", "p1.bed", package = "CfDNAfragmentomics")
#' controls_bed <- read.table(controls_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
#'
#' sample_bed_file <- system.file("extdata", "d1.bed", package = "CfDNAfragmentomics")
#' sample_bed <- read.table(sample_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
#'
#' ## basic usage
#' plot <- nucleosomeDensityPlot(controls_bed = controls_bed, sample_bed = sample_bed)
#'
#'
#'}
#'
#'
#' @export
#' @import ggplot2
#' @importFrom gridExtra grid.arrange

nucleosomeDensityPlot <- function(controls_bed, sample_bed) {

  #check if the input values are correct.
  if (! is.data.frame(controls_bed) | ! is.data.frame(sample_bed)) {
    stop("Please put valid inputs for nucleosomeRatio function")
  }


  #assume the columns 2 and 3 are numeric values are mentioned in
  #function description

  #find fragment lengths of control data
  controls_bed_lengths <- data.frame(controls_bed[3] - controls_bed[2])
  colnames(controls_bed_lengths)[1] = "length"

  #find fragment lengths of patient data
  patient_bed_lengths <- data.frame(sample_bed[3] - sample_bed[2])
  colnames(patient_bed_lengths)[1] = "length"


  # plot histogram of fragment sizes for the patient and control data to
  #visualize fragment sizes differences.

  #library(ggplot2)
  plot_control <- ggplot2::ggplot(controls_bed_lengths, mapping=aes(x=length)) +
         geom_density() +
    ggtitle("Controls")

  plot_sample <- ggplot2::ggplot(patient_bed_lengths, mapping=aes(x=length)) +
    geom_density(color="red") +
    ggtitle("Patients")

  gridExtra::grid.arrange(plot_control, plot_sample, nrow=2)
}

# [END]
