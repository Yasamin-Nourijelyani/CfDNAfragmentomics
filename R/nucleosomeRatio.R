# Purpose: Analysis of cfDNA fragment lengths of patient data as compared with
# control data using t-test.
# Author: Yasamin Nouri Jelyani
# Date: 2022-12-05
# Version: 0.1.0
# Bugs and Issues:


#------------------------Helper Function------------------------


#' Helper function: Detecting dinucleosome lengths of controls
#'
#' @param control_frag_sizes A dataframe for the a control data
#' (healthy individual) cfDNA fragment lengths.
#' This data frame has 1 column, showing the cfDNA fragment lengths for
#' one control individual.
#'
#' @return Return the control length dataframe with only dinucleosome sizes
#' as defined in the research by Katsman et.al (1).
#'
#' @references
#' 1- Katsman, E., Orlanski, S., Martignano, F., Fox-Fisher, I., Shemer,
#'  R., Dor, Y., Zick, A., Eden, A., Petrini, I., Conticello, S. G.,
#' Berman, B. P. (2022). Detecting cell-of-origin and cancer-specific
#' methylation features of cell-free DNA from Nanopore sequencing.
#' Genome biology, 23(1), 158.
#' https://doi.org.myaccess.library.utoronto.ca/10.1186/s13059-022-02710-1
#'
get_controles_di <- function(control_frag_sizes){

  #dinucleosome ratio: (275-400 base pairs (Katsman et. al.)):
  #filter the dataframe rows for dinucleosome lengths

  subject_1 <- dplyr::filter(control_frag_sizes, 275 <= length & length<= 400)

  return(subject_1)
}



#------------------------Helper Function------------------------


#' Helper function: Detecting mononucleosome lengths of controls
#'
#' @param control_frag_sizes A dataframe for the a control data
#' (healthy individual) cfDNA fragment lengths.
#' This data frame has 1 column, showing the cfDNA fragment lengths for
#' one control individual.
#'
#'@return Return the control length data frame with only mononucleosome sizes
#'as defined in the research by Katsman et.al (1).
#'
#' @references
#' 1- Katsman, E., Orlanski, S., Martignano, F., Fox-Fisher, I., Shemer,
#'  R., Dor, Y., Zick, A., Eden, A., Petrini, I., Conticello, S. G.,
#' Berman, B. P. (2022). Detecting cell-of-origin and cancer-specific
#' methylation features of cell-free DNA from Nanopore sequencing.
#' Genome biology, 23(1), 158.
#' https://doi.org.myaccess.library.utoronto.ca/10.1186/s13059-022-02710-1
#'
get_controles_mono <- function(control_frag_sizes){

  #mononucleosome ratio: (<150 base pairs (Katsman et. al.)):
  #filter the dataframe rows for mononucleosome lengths

  subject_1 <- dplyr::filter(control_frag_sizes, length < 150)

  return(subject_1)
}

#------------------------Helper Function------------------------


#'Helper function: Detecting dinucleosome lengths of patient
#'
#' @param sample_frag_sizes A dataframe for a patient data cfDNA
#' fragment lengths. This data frame has 1 column, showing the cfDNA fragment
#' lengths for one patient.
#'
#' @return Return the patient length data frame with only dinucleosome sizes
#' as defined in the research by Katsman et.al (1).
#'
#' @references
#' 1- Katsman, E., Orlanski, S., Martignano, F., Fox-Fisher, I., Shemer,
#'  R., Dor, Y., Zick, A., Eden, A., Petrini, I., Conticello, S. G.,
#' Berman, B. P. (2022). Detecting cell-of-origin and cancer-specific
#' methylation features of cell-free DNA from Nanopore sequencing.
#' Genome biology, 23(1), 158.
#' https://doi.org.myaccess.library.utoronto.ca/10.1186/s13059-022-02710-1
#'
#'
get_sample_di <- function(sample_frag_sizes){


  #dinucleosome ratio: (275-400 base pairs (Katsman et. al.)):
  #filter the dataframe rows for dinucleosome lengths

  subject_1 <- dplyr::filter(sample_frag_sizes, 275 <= length & length<= 400)

  return(subject_1)
}


#------------------------Helper Function------------------------


#'Helper function: Detecting mononucleosome lengths of patient
#'
#' @param sample_frag_sizes A dataframe for a patient data
#' fragment lengths. This data frame has 1 column, showing the
#' cfDNA fragment lengths for one patient.
#'
#'@return Return the control length data frame with only mononucleosome sizes
#'as defined in the research by Katsman et.al (1).
#'
#' @references
#' 1- Katsman, E., Orlanski, S., Martignano, F., Fox-Fisher, I., Shemer,
#'  R., Dor, Y., Zick, A., Eden, A., Petrini, I., Conticello, S. G.,
#' Berman, B. P. (2022). Detecting cell-of-origin and cancer-specific
#' methylation features of cell-free DNA from Nanopore sequencing.
#' Genome biology, 23(1), 158.
#' https://doi.org.myaccess.library.utoronto.ca/10.1186/s13059-022-02710-1
#'
get_sample_mono <- function(sample_frag_sizes){

  #mononucleosome:(< 150 base pairs (Katsman et. al.))
  #filter the dataframe rows for mononucleosome lengths

  subject_1 <- dplyr::filter(sample_frag_sizes, length < 150)

  return(subject_1)
}


#------------------------Exported Function------------------------



#' Wilcox-Test of nucleosome ratios for patient and control cfDNA fragment
#' reads data.
#'
#' Fragmentation sizes of cfDNA (cell free DNA) molecules are potential
#' cancer biomarkers (1).
#'  Hence, to find if a patient data file contains cancer information,
#' it can be compared with a control dataset with known healthy cfDNA fragments.
#' The Wilcox test is a non-parametric test used to determine
#' significance of the difference in fragment
#' size of the control and patient cfDNA lengths to determine if the patient
#' data contains the cancer biomarker.
#' Here, we are comparing two independent samples, which are the
#' control and patient data, and we are making inference about the state
#' (cancer positive or negative) of the
#' population of cfDNA molecules in the patient as being potentially cancer
#'  positive or negative. Both control and patient
#' data inputted to the functions are assumed to be reads for the specific loci
#' corresponding to the cancer type of interest. Hence, both the control
#' and patient data are required in the analysis to ensure the loci
#' being analysed are corresponding to the specific cancer type of interest.
#' Additionally, performing a non-parametric wilcox-test does not assume that the data is
#' normally distributed (5).
#'  It is also assumed that the data is real human data and contains
#'  both mono-nucleosome and di-nucleosome length data.
#'
#' Note 1: The data to generate examples are from the illumina website (4).
#'  This example data is not representative of real patient data due to
#'  limitations to access of real patient data.
#'  However, users can use their own real patient datasets for
#'  analysis.
#'  More examples can be found at https://github.com/Yasamin-Nourijelyani/CfDNAfragmentomics/tree/master/inst/extdata
#'
#' Note 2:
#' BED files can be generate from BAM files using bedtools bamtobed function
#' in the command line (3).
#'
#' Since BED files are commonly
#' generated by computational biologists and are easily convertible, from
#' BAM files, they can be used as input. Otherwise, tab seperated
#' .txt files that generate data in the same format of bed files
#' with the same columns as described in the parameters are also
#' excepted.
#'
#' @param controls_bed A dataframe for the a control data (healthy individual)
#' with at least 3 columns: first column is a string representing the chromosome
#'  number for the cfDNA for instance: "chr1",
#'  the second column includes integer values indicating the start location of
#'  the read and third column includes integer values indicating end location of
#'  the read.
#' These are the healthy individuals who do not have cancer.The fragment lengths
#' of patients will be compared to these healthy data similar to what is done in
#' the paper by Katsman et. al (1).
#' @param sample_bed A dataframe for the a sample data
#' (unknown if cancerous or not individual)
#'  with at least 3 columns: first column is a string representing the chromosome
#'  number for the cfDNA for instance: "chr1",
#'  the second column are integer values indicating the start location of
#'  the read and third column are integer values indicating end location of
#'  the read.
#' These are the patient who we want to check if they have cancer or not.
#' The fragment lengths
#' of patients will be compared to these healthy data similar to what is done in
#' the paper by Katsman et. al (1).
#' These sample data are compared to the controls using a t-test.
#'
#' @param di_nucleosome_p_value A numeric value between 0 and 1
#' indicating the p-value for the dinuclesome ratio that is considered
#' significant. Default value is 0.05.
#' @param mono_nucleosome_p_value A numeric value between 0 and 1
#' indicating the p-value for the dinuclesome ratio that is considered
#' significant. Default value is 0.05
#'
#' @return Return a list of outputs:
#' The function returns a list of outputs.
#' This list contains a real number: ttest_mono_pvalue, a real number: ttest_di_pvalue,
#' and a boolean value: cancerous.
#' ttest_mono_pvalue is the p-value calculated for the difference in mononucleosome
#' lengths between patient and controls, ttest_di_pvalue is the p-value calculated
#' for the difference in dinucleosome lengths, and if the patient input files shows significant shorter
#' mononucleosome and dinucleosome sizes compared to the control with respect to the
#' inputted p-values that are considered significant, cancerous will be TRUE.
#' and cancerous will be FALSE otherwise. Note that even if the p-values show significance,
#' if the patient cfDNA lengths are not shorter, then the patient is not considered as being cancerous.
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
#' 2- “Mann Whitney U Test in R Programming.” GeeksforGeeks, 28 Dec. 2021,
#' https://www.geeksforgeeks.org/mann-whitney-u-test-in-r-programming/.
#'
#' 3- Quinlan, Aaron R., Kindlon, Neil. "bamtobed" readthedocs, May 26, 2022,
#'  https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html
#'
#' 4- Illumina, https://support.illumina.com/downloads/nextera-flex-for-enrichment-BED-files.html
#'
#'
#'@examples
#' \dontrun{
#'
#' controls_bed_file <- system.file("extdata", "p1.bed", package = "CfDNAfragmentomics")
#' controls_bed <- read.table(controls_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
#'
#' sample_bed_file <- system.file("extdata", "d1.bed", package = "CfDNAfragmentomics")
#' sample_bed <- read.table(sample_bed_file, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
#'
#' ## basic usage
#'
#' # Example 1:
#'
#' ratio <- nucleosomeRatio(controls_bed = controls_bed, sample_bed = sample_bed, 0.3, 0.3)
#' ratio$cancerous
#'
#' # Example 2:
#'
#' ratio <- nucleosomeRatio(controls_bed = controls_bed, sample_bed = sample_bed)
#' ratio$cancerous
#'
#' }
#'
#'
#' @export
#' @importFrom stats t.test
#' @importFrom dplyr between

nucleosomeRatio <- function(controls_bed, sample_bed, di_nucleosome_p_value=0.05,
                            mono_nucleosome_p_value=0.05){

  #check if the input values are correct.
  if (! is.data.frame(controls_bed) | ! is.data.frame(sample_bed) |
      ! dplyr::between(di_nucleosome_p_value, 0, 1) |
      ! dplyr::between(mono_nucleosome_p_value, 0, 1) | ncol(controls_bed) < 3 |
      ncol(sample_bed) < 3 | (! is.numeric(controls_bed[, 2])) |
      (! is.numeric(controls_bed[, 3])) | (! is.numeric(sample_bed[, 2])) |
      (! is.numeric(sample_bed[, 3])) | ! is.numeric(di_nucleosome_p_value) |
      ! is.numeric(mono_nucleosome_p_value)) {
    stop("Please put valid inputs for nucleosomeRatio function")
  }




  #assume the columns 2 and 3 are numeric values are mentioned in
  #function description

  #find fragment lengths of control data
  controls_bed_lengths <- data.frame(controls_bed[, 3] - controls_bed[, 2])
  colnames(controls_bed_lengths)[1] = "length"

  #find fragment lengths of patient data
  patient_bed_lengths <- data.frame(sample_bed[, 3] - sample_bed[, 2])
  colnames(patient_bed_lengths)[1] = "length"


  #--------------------Controls fragment sizes--------------------------------

  #----------------di---------------

  #df for ratio of dinucleosomes for each control subject
  controls_di_df <- get_controles_di(controls_bed_lengths)
  #----------------MONO------------

  #df for ratio of mononucleosomes for each control subject
  controls_mono_df <- get_controles_mono(controls_bed_lengths)


  #--------------------sample fragment sizes------------------------------------

  #----------------di--------------

  #df for ratio of dinucleosomes for each sample subject
  patient_di_df <- get_sample_di(patient_bed_lengths)

  #----------------MONO------------

  #df for ratio of mononucleosomes for each sample subject
  patient_mono_df <- get_sample_mono(patient_bed_lengths)



  #--------------------stat test for mono- and di- nucleosome ratios----------

  # perform mann-whitney-wilcoxen rank sum test (non-parametric)
  #to check if there is significant difference between
  # patient and control mono and di nucleosome fragment sizes.
  # Fragment lengths are useful biomarkers in the detection of cancer.

  if (nrow(controls_mono_df) != 0 & nrow(patient_mono_df) != 0){
    w_test_mono <- stats::wilcox.test(as.numeric(controls_mono_df$length),
                                      as.numeric(patient_mono_df$length))

  }
  else{
    stop("There are no mono-nucleosome data in this file. Please
         input valid patient cfDNA data")


  }

  if (nrow(controls_di_df) != 0 & nrow(patient_di_df) != 0){
    w_test_di <- stats::wilcox.test(as.numeric(controls_di_df$length),
                                    as.numeric(patient_di_df$length))


  }
  else{
    stop("There are no di-nucleosome data in this file. Please
         input valid patient cfDNA data")


  }




  # To determine if the patient data is potentially cancerous, given
  # their cfDNA fragment lengths, we will check if the
  # fragment sizes are significantly shorter than control.

  if (w_test_mono$p.value < mono_nucleosome_p_value &
      w_test_di$p.value < di_nucleosome_p_value){

    # there is significant difference in the sizes of
    # patient and control mono- and di- nucleosomes

    # check if patient fragment lengths for mono and di nucleosomes
    # are shorter, which is a biomarker for cancer.

    if ((mean(controls_mono_df$length) >
         mean(patient_mono_df$length)) &
        (mean(controls_di_df$length) >
         mean(patient_di_df$length))){

      #if the average sizes of the patient are shorter (significantly)
      #compared the control,then the sample is likely cancerous,
      #so return true
      Results <- list(wtest_mono_pvalue = w_test_mono$p.value,
                      wtest_di_pvalue = w_test_di$p.value,
                      cancerous = TRUE)

      return(Results)

    }else{ #the sample lengths are not significantly shorter
      Results <- list(wtest_mono_pvalue = w_test_mono$p.value,
                      wtest_di_pvalue = w_test_di$p.value,
                      cancerous = FALSE)

      return(Results)
    }

  }else { #the difference between sample and control is not significant

    Results <- list(wtest_mono_pvalue = w_test_mono$p.value,
                    wtest_di_pvalue = w_test_di$p.value,
                    cancerous = FALSE)

    return(Results)  }
}




# [END]

